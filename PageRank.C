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
#if NUMA
#include "ligra-numa.h"
#else
#include "ligra.h"
#endif
#include "math.h"
int MaxIter=10;
template <class vertex>
struct PR_F
{
    double* p_curr, *p_next;
    vertex* V;
    static const bool use_cache = true;
    struct cache_t
    {
        double p_next;
    };
    PR_F(double* _p_curr, double* _p_next, vertex* _V) :
        p_curr(_p_curr), p_next(_p_next), V(_V) {}
    inline bool update(intT s, intT d)  //update function applies PageRank equation
    {
        p_next[d] += p_curr[s]/V[s].getOutDegree();
        return 1;
    }
    inline bool updateAtomic (intT s, intT d)   //atomic Update
    {
        writeAdd(&p_next[d],p_curr[s]/V[s].getOutDegree());
        return 1;
    }

    inline void create_cache(cache_t &cache, intT d)
    {
        cache.p_next = p_next[d];
    }
    inline bool update(cache_t &cache, intT s)
    {
        cache.p_next += p_curr[s]/V[s].getOutDegree();
        return 1;
    }

    inline void commit_cache(cache_t &cache, intT d)
    {
	// Cache is used only in sequential mode
	p_next[d] = cache.p_next;
    }

    inline bool cond (intT d)
    {
        return cond_true(d);    //does nothing
    }
};

//vertex map function to update its p value according to PageRank equation
struct PR_Vertex_F
{
    double damping;
    double addedConstant;
    double* p_curr;
    double* p_next;
    PR_Vertex_F(double* _p_curr, double* _p_next, double _damping, intT n) :
        damping(_damping), addedConstant((1-_damping)*(1/(double)n)),
        p_curr(_p_curr), p_next(_p_next) {}
    inline bool operator () (intT i)
    {
        p_next[i] = damping*p_next[i] + addedConstant;
        return 1;
    }
};

//resets p
struct PR_Vertex_Reset
{
    double* p_curr;
    PR_Vertex_Reset(double* _p_curr) :
        p_curr(_p_curr) {}
    inline bool operator () (intT i)
    {
        p_curr[i] = 0.0;
        return 1;
    }
};

template <class GraphType>
void Compute(GraphType &GA, long start)
{
    typedef typename GraphType::vertex_type vertex; // Is determined by GraphTyp
    const partitioner &part = GA.get_partitioner();
    const int perNode = part.get_num_per_node_partitions();
    intT n = GA.n;
    intT m = GA.m;
    const double damping = 0.85;
    const double epsilon = 0.0000001;
    //Data Array
    //p_curr and p_next to do special allocation
    //frontier also need special node allocation
    //blocksize equal to the szie of each partitioned
    double one_over_n = 1/(double)n;
    
    mmap_ptr<double> p_curr;
    p_curr.part_allocate (part);
    mmap_ptr<double> p_next;
    p_next.part_allocate (part);
    loop(j,part,perNode,p_curr[j]=one_over_n);
    loop(j,part,perNode,p_next[j]=0);

    int count=0;
    partitioned_vertices Frontier = partitioned_vertices::bits(part,n, m);
    while(1 && count<MaxIter)
    {
        count++;
        partitioned_vertices output = edgeMap(GA, Frontier, PR_F<vertex>(p_curr,p_next,GA.getWholeGraph().V),m/20);
        vertexMap(part,Frontier,PR_Vertex_F(p_curr,p_next,damping,n));
        //compute L1-norm between p_curr and p_next
        {
            loop(j,part,perNode,p_curr[j] = fabs(p_curr[j]-p_next[j]));
        }
        double L1_norm = sequence::plusReduce(p_curr.get(),n);
        if(L1_norm < epsilon) break;
        //reset p_curr
        vertexMap(part, Frontier,PR_Vertex_Reset(p_curr));
        swap(p_curr,p_next);
        Frontier.del();
        Frontier = output;
        Frontier.bit=true;
    }
    Frontier.del();
    p_curr.del();
    p_next.del();
}

