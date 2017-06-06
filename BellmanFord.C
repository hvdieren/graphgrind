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
#define WEIGHTED 1
#define VERTEX 1
#if NUMA
#include "ligra-numa.h"
#else
#include "ligra.h"
#endif

struct BF_F
{
    intT* ShortestPathLen;
    intT* Visited;
    struct cache_t
    {
        intT shotestPathLen;
        intT visited;
    };
    static const bool use_cache = false;
    BF_F(intT* _ShortestPathLen, intT* _Visited) :
        ShortestPathLen(_ShortestPathLen), Visited(_Visited) {}
    inline bool update (intT s, intT d, intE edgeLen)   //Update ShortestPathLen if found a shorter path
    {
        intT newDist = ShortestPathLen[s] + edgeLen;
        if(ShortestPathLen[d] > newDist)
        {
            ShortestPathLen[d] = newDist;
            if(Visited[d] == 0)
            {
                Visited[d] = 1 ;
                return 1;
            }
        }
        return 0;
    }
    inline bool updateAtomic (intT s, intT d, intE edgeLen)  //atomic Update
    {
        intT newDist = ShortestPathLen[s] + edgeLen;
        return (writeMin(&ShortestPathLen[d],newDist) &&
                CAS(&Visited[d],intT(0),intT(1)));
    }

    inline void create_cache(cache_t &cache, intT d)
    {
        cache.shotestPathLen = ShortestPathLen[d];
        cache.visited = Visited[d];
    }
    inline bool update(cache_t &cache, intT s, intE edgeLen)
    {
        cache.shotestPathLen = ShortestPathLen[s] + edgeLen;
        if(cache.visited == 0)
        {
            cache.visited = 1 ;
            return 1;
        }
        return 0;
    }

    inline intT commit_cache(cache_t &cache, intT d)
    {
        writeMin(&ShortestPathLen[d],cache.shotestPathLen);
        return Visited[d]=cache.visited;
    }

    inline bool cond (intT d)
    {
        return cond_true(d);    //does nothing
    }
};

//reset visited vertices
struct BF_Vertex_F
{
    intT* Visited;
    BF_Vertex_F(intT* _Visited) : Visited(_Visited) {}
    inline bool operator() (intT i)
    {
        Visited[i] = 0;
        return 1;
    }
};

template <class GraphType>
void Compute (GraphType &GA, long start)
{
    intT n = GA.n;
    intT m = GA.m;
    const partitioner &part = GA.get_partitioner();
    //initialize ShortestPathLen to "infinity"
    const int perNode = part.get_num_per_node_partitions();

    mmap_ptr<intT> ShortestPathLen;
    ShortestPathLen.part_allocate (part);
    
    loop(j,part,perNode,ShortestPathLen[j] = INT_MAX/2);
    ShortestPathLen[start] = 0;
  
    mmap_ptr<intT> Visited;
    Visited.part_allocate (part);
    loop(j,part,perNode,Visited[j] = 0);
    //vertexSubset Frontier(n,start); //initial frontier
    partitioned_vertices Frontier=partitioned_vertices::create(n,start,GA.getWholeGraph().V[start].getOutDegree());

    intT round = 0;
    while(!Frontier.isEmpty())
    {
        if(round == n)
        {
            //negative weight cycle
            {
                loop(j,part,perNode,ShortestPathLen[j] = -(INT_MAX/2));
            }
            break;
        }
        partitioned_vertices output = edgeMap(GA,Frontier,BF_F(ShortestPathLen,Visited), m/20,DENSE_FORWARD);
        vertexMap(part, output,BF_Vertex_F(Visited));
        Frontier.del();
        Frontier = output;
        round++;
    }
    Frontier.del();
    Visited.del();
    ShortestPathLen.del();
}
