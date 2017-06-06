// This code is part of the project "Ligra: A Lightweight Graph Processing
// Framework for Shared Memory", presented at Principles and Practice of
// Parallel Programming, 2013.
// Copyright (c) 2013 Julian Shun and Guy Blelloch
//
// This code has been extended in the project "GraphGrind: Addressing Load
// Imbalance of Graph Partitioning", presented at International Symposium
// on Supercomputing, 2017.
// Copyright (c) 2017 Jiawen Sun, Hans Vandierendonck, Dimitrios S. Nikolopoulos
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights (to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject to
// the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
// LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
// OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
// WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
//#define PUSH 1
#if NUMA
#include "ligra-numa.h"
#else
#include "ligra.h"
#endif
#include "math.h"
static intT CAScounter;
template <class vertex>
struct PR_F
{
    vertex* V;
    double* Delta, *nghSum;
    static const bool use_cache = false;
    struct cache_t
    {
        double nghSum;
    };
    PR_F(vertex* _V, double* _Delta, double* _nghSum) :
        V(_V), Delta(_Delta), nghSum(_nghSum) {}
    inline bool update(intT s, intT d)
    {
        nghSum[d] += Delta[s]/V[s].getOutDegree();
        return 1;
    }
    inline bool updateAtomic (intT s, intT d)
    {
        //__sync_fetch_and_add(&CAScounter,1);
        writeAdd(&nghSum[d],Delta[s]/V[s].getOutDegree());
        return 1;
    }

    inline void create_cache(cache_t &cache, intT d)
    {
        cache.nghSum = nghSum[d];
    }
    inline bool update(cache_t &cache, intT s)
    {
        cache.nghSum = Delta[s]/V[s].getOutDegree();
        return 1;
    }

    inline void commit_cache(cache_t &cache, intT d)
    {
        writeAdd(&nghSum[d],cache.nghSum);
    }

    inline bool cond (intT d)
    {
        return cond_true(d);
    }
};

struct PR_Vertex_F_FirstRound
{
    double damping, addedConstant, one_over_n, epsilon2;
    double* p, *Delta, *nghSum;
    PR_Vertex_F_FirstRound(double* _p, double* _Delta, double* _nghSum, double _damping, double _one_over_n,double _epsilon2) :
        p(_p),
        damping(_damping), Delta(_Delta), nghSum(_nghSum), one_over_n(_one_over_n),
        addedConstant((1-_damping)*_one_over_n),
        epsilon2(_epsilon2) {}
    inline bool operator () (intT i)
    {
        Delta[i] = damping*(p[i]+nghSum[i])+addedConstant-p[i];
        p[i] += Delta[i];
        Delta[i]-=one_over_n; //subtract off delta from initialization
        return (fabs(Delta[i]) > epsilon2 * p[i]);
    }
};

struct PR_Vertex_F
{
    double damping, epsilon2;
    double* p, *Delta, *nghSum;
    PR_Vertex_F(double* _p, double* _Delta, double* _nghSum, double _damping, double _epsilon2) :
        p(_p),
        damping(_damping), Delta(_Delta), nghSum(_nghSum),
        epsilon2(_epsilon2) {}
    inline bool operator () (intT i)
    {
        Delta[i] = nghSum[i]*damping;
        p[i]+=Delta[i];
        return (fabs(Delta[i]) > epsilon2*p[i]);
    }
};

struct PR_Vertex_Reset
{
    double* nghSum;
    PR_Vertex_Reset(double* _nghSum) :
        nghSum(_nghSum) {}
    inline bool operator () (intT i)
    {
        nghSum[i] = 0.0;
        return 1;
    }
};

template <class GraphType>
void Compute(GraphType &GA, long start)
{
    typedef typename GraphType::vertex_type vertex; // Is determined by GraphType
    const partitioner &part = GA.get_partitioner();
    const int perNode = part.get_num_per_node_partitions();
    const double damping = 0.85;
    const double epsilon = 0.0000001;
    const double epsilon2 = 0.01;
    intT n = GA.n;
    intT m = GA.m;

    double one_over_n = 1/(double)n;
    
    mmap_ptr<double> p;
    p.part_allocate (part);
    mmap_ptr<double> nghSum;
    nghSum.part_allocate (part);
    mmap_ptr<double> Delta;
    Delta.part_allocate (part);
    
    loop(j,part,perNode, p[j]=0.0);
    loop(j,part,perNode, Delta[j]=one_over_n);
    loop(j,part,perNode, nghSum[j]=0.0);

    partitioned_vertices Frontier = partitioned_vertices::bits(part,n, m);
    partitioned_vertices All = partitioned_vertices::bits(part,n,m);
 //   CAScounter=0; 
    intT round = 0;
    while(1)
    {
        round++;
        partitioned_vertices output = edgeMap(GA,Frontier,PR_F<vertex>(GA.getWholeGraph().V,Delta,nghSum),m/20,DENSE_FORWARD);
        output.del();
        //vertexSubset active
        partitioned_vertices active
            = (round == 1) ?
              vertexFilter(GA,All,PR_Vertex_F_FirstRound(p,Delta,nghSum,damping,one_over_n,epsilon2)) :
              vertexFilter(GA,All,PR_Vertex_F(p,Delta,nghSum,damping,epsilon2));
        //compute L1-norm (use nghSum as temp array)
        {
            loop(j,part,perNode,nghSum[j] = fabs(Delta[j]));
        }
        double L1_norm = sequence::plusReduce(nghSum.get(),n);
        if(L1_norm < epsilon) break;
        //reset
        vertexMap(part,All,PR_Vertex_Reset(nghSum));
        Frontier.del();
        Frontier = active;
    }
 //   cerr<<"CAScounter"<<CAScounter<<endl;
    Frontier.del();
    p.del();
    Delta.del();
    nghSum.del();
    All.del();
}
