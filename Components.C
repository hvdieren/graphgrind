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

struct CC_F
{
    intT* IDs;
    intT* prevIDs;
    static const bool use_cache = true;
    struct cache_t
    {
        intT ID_d, prevID_d;
    };
    CC_F(intT* _IDs, intT* _prevIDs) :
        IDs(_IDs), prevIDs(_prevIDs) {}
    inline bool update(intT s, intT d)  //Update function writes min ID
    {
        intT origID = IDs[d];
        // cerr<<"the origID["<<d<<"] is "<<IDs[d]<<" the IDs["<<s<<"] is "<<IDs[s]<<endl;
        if(IDs[s] < origID)
        {
            IDs[d] = min(origID,IDs[s]);
            if(origID == prevIDs[d]) return 1;
        }
        return 0;
    }

    inline bool updateAtomic (intT s, intT d)   //atomic Update
    {
        intT origID = IDs[d];
        return (writeMin(&IDs[d],IDs[s])
                && origID == prevIDs[d]);
    }
    inline void create_cache(cache_t &cache, intT d)
    {
        cache.ID_d = IDs[d];
        cache.prevID_d = prevIDs[d];
    }
    inline bool update(cache_t &cache, intT s)
    {
        intT origID = cache.ID_d;
        if(IDs[s]<origID)
        {
            cache.ID_d = min(origID, IDs[s]);
            if (origID==cache.prevID_d) return 1;
        }
        return 0;
    }

    inline void commit_cache(cache_t &cache, intT d)
    {
	// The cache is used only when d is accessed by a single thread
	IDs[d] = cache.ID_d;
    }

    inline bool cond (intT d)
    {
        return cond_true(d);    //does nothing
    }
};

//function used by vertex map to sync prevIDs with IDs
struct CC_Vertex_F
{
    intT* IDs;
    intT* prevIDs;
    CC_Vertex_F(intT* _IDs, intT* _prevIDs) :
        IDs(_IDs), prevIDs(_prevIDs) {}
    inline bool operator () (intT i)
    {
        prevIDs[i] = IDs[i];
        return 1;
    }
};

template <class GraphType>
void Compute(GraphType &GA, long start)
{
    typedef typename GraphType::vertex_type vertex; // Is determined by GraphType
    const partitioner &part = GA.get_partitioner();
    const int perNode = part.get_num_per_node_partitions();
    intT n = GA.n;
    intT m = GA.m;

    mmap_ptr<intT> IDs;
    IDs.part_allocate (part);
    mmap_ptr<intT> prevIDs;
    prevIDs.part_allocate (part);
    loop(j,part,perNode,IDs[j]=j);

    partitioned_vertices Frontier = partitioned_vertices::bits(part,n,m);

    while(!Frontier.isEmpty())  //iterate until IDS converge
    {
        vertexMap(part,Frontier,CC_Vertex_F(IDs,prevIDs));
        partitioned_vertices output=edgeMap(GA, Frontier, CC_F(IDs,prevIDs), m/20);
        Frontier.del();
        Frontier = output;
    }
    Frontier.del();
    IDs.del();
    prevIDs.del();
}

