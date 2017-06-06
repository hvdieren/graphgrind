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
#define VERTEX 1
#if NUMA
#include "ligra-numa.h"
#else
#include "ligra.h"
#endif
//static intT CAScounter;
struct BFS_F
{
    struct cache_t
    {
        intT parents;
    };
    static const bool use_cache = false;
    intT* Parents;
    BFS_F(intT* _Parents) : Parents(_Parents) {}
    inline bool update (intT s, intT d)   //Update
    {
        if(Parents[d] == -1)
        {
            Parents[d] = s;
            return 1;
        }
        else return 0;
    }
    inline bool updateAtomic (intT s, intT d)  //atomic version of Update
    {
        //__sync_fetch_and_add(&CAScounter,1);
        return (CAS(&Parents[d],(intT)-1,s));
    }
    //For partitioned graph used
    inline void create_cache(cache_t &cache, intT d)
    {
        cache.parents=Parents[d];
    }

    inline bool update(cache_t &cache, intT s)
    {
        if (cache.parents==-1)
        {
            cache.parents=s;
            return 1;
        }
        else return 0;
    }

    inline void commit_cache(cache_t &cache, intT d)
    {
        Parents[d]=cache.parents;
    }
    //cond function checks if vertex has been visited yet
    inline bool cond (intT d)
    {
        return (Parents[d] == -1);
    }
};

template <class GraphType>
void Compute(GraphType &GA, long start)
{
    typedef typename GraphType::vertex_type vertex; // Is determined by GraphType
    const intT n = GA.n;
    const intT m = GA.m;
    const partitioner &part = GA.get_partitioner();
    const int perNode = part.get_num_per_node_partitions();
//creates Parents array, initialized to all -1, except for start

    mmap_ptr<intT> Parents;
    Parents.part_allocate(part);
    loop(j,part,perNode,Parents[j]=-1);
    Parents[start] = start;
// creates initial frontier
    partitioned_vertices Frontier=partitioned_vertices::create(n,start, GA.getWholeGraph().V[start].getOutDegree());
    while(!Frontier.isEmpty())  //loop until frontier is empty
    {
        partitioned_vertices output=edgeMap(GA,Frontier,BFS_F(Parents), m/20);
        Frontier.del();
        Frontier = output; //set new frontier
    }
    Frontier.del();
    Parents.del();
}
