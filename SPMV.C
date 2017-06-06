#define WEIGHTED 1
#if NUMA
#include "ligra-numa.h"
#else
#include "ligra.h"
#endif
#include "math.h"
//static intT CAScounter; 
int maxIter=1;
template <class vertex>
struct SPMV_F
{
    double* p_curr, *p_next;
    static const bool use_cache = false;
    struct cache_t
    {
        double p_curr;
    };
    SPMV_F(double* _p_curr, double* _p_next) :
        p_curr(_p_curr), p_next(_p_next) {}
    inline bool update(intT s, intT d, int edgeLen)
    {
        p_next[d] += p_curr[s] * edgeLen;
        return 1;
    }
    inline bool updateAtomic (intT s, intT d,int edgeLen)
    {
        writeAdd(&p_next[d], p_curr[s]*edgeLen);
      //  __sync_fetch_and_add(&CAScounter,1);
        return 1;
    }

    inline void create_cache(cache_t &cache, intT d)
    {
        cache.p_curr = p_next[d];
    }
    inline bool update(cache_t &cache, intT s,int edgeLen)
    {
        cache.p_curr += p_curr[s]*edgeLen;
        return 1;
    }

    inline void commit_cache(cache_t &cache, intT d)
    {
        writeAdd(&p_next[d],cache.p_curr);
    }

    inline bool cond (intT d)
    {
        return cond_true(d);
    }
};


//resets p
struct SPMV_Vertex_Reset
{
    double* p_curr;
    SPMV_Vertex_Reset(double* _p_curr) :
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
    typedef typename GraphType::vertex_type vertex; // Is determined by GraphType
    const partitioner &part = GA.get_partitioner();
    const int perNode = part.get_num_per_node_partitions();
    intT n = GA.n;
    intT m = GA.m;

    double one_over_n = 1/(double)n;

    mmap_ptr<double> p_curr;
    p_curr.part_allocate(part);
    mmap_ptr<double> p_next;
    p_next.part_allocate(part);
   
    loop(j,part,perNode,p_curr[j] = one_over_n);
    loop(j,part,perNode,p_next[j] = 0);
    //CAScounter=0;
    partitioned_vertices Frontier = partitioned_vertices::bits(part,n, m);
    partitioned_vertices output = edgeMap(GA, Frontier, SPMV_F<vertex>(p_curr,p_next),m/20,DENSE_FORWARD);
    //cerr<<"CAScounter:"<<CAScounter<<endl;
    vertexMap(part,Frontier,SPMV_Vertex_Reset(p_curr));
    swap(p_curr,p_next);
    Frontier.del();
    output.del();
    p_curr.del();
    p_next.del();
}

