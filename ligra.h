// -*- C++ -*-
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

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <cstring>
#include <string>
#include <assert.h>
#include <algorithm>
#include <sys/mman.h>

#include "parallel.h"
#include "gettime.h"
#include "utils.h"
#include "graph.h"
#include "IO.h"
#include "parseCommandLine.h"
#include "gettime.h"
#include "papi_code.h"

#if REUSE_DIST
#include "reuse.h"
RDLogger<intT> * rd_bwd_src = 0;  // old frontier
RDLogger<intT> * rd_bwd_data = 0;  // old data
RDLogger<intT> * rd_bwd_next = 0;  // next frontier
RDLogger<intT> * rd_fwd_dest = 0; // next frontier
RDLogger<intT> * rd_fwd_new = 0; // new data
#endif


//Timer design
using namespace std;
static double tm_edgemap_ = 0;
static double tm_edgemap_tosparse_alloc_ = 0;
static double tm_edgemap_tosparse_collect_ = 0;
static double tm_edgemap_tosparse_itself_ = 0;
static double tm_edgemap_tosparse_reduce_ = 0;
static double tm_edgemap_todense_ = 0;
static double tm_edgemap_dense_fwd_ = 0;
static double tm_edgemap_dense_ = 0;
static double tm_edgemap_dense_bwd_ = 0;
static double tm_edgemap_dense_sync_ = 0;
static double tm_edgemap_sparse_ = 0;
static double tm_edgemap_sparse_sync_ = 0;
static double dm_dense_each=0;
static double sm_sparse_each=0;
static double fm_forward_each=0;

static double tm_edgemap_setup_ = 0;
static double tm_edgemap_map_ = 0;
static double tm_edgemap_collect_ = 0;

double tmlog( timer & tm, double & cnt )
{
    double val = tm.next();
    cnt += val;
    return val;
}

void timePrint()
{
    std::cerr << "\ntotal edgemap setup: " << tm_edgemap_setup_
              << "\ntotal edgemap map(dense): " << tm_edgemap_dense_
              << "\ntotal edgemap map(sparse): " << tm_edgemap_sparse_
              << "\ntotal edgemap collect: " << tm_edgemap_collect_
              << std::endl;
}

void timeprint()
{
    std::cerr << "edgemap: " << tm_edgemap_
              << " edgemap_tosparse_alloc: " << tm_edgemap_tosparse_alloc_ <<"\n"
              << " edgemap_tosparse_collect: " << tm_edgemap_tosparse_collect_ <<"\n"
              << " edgemap_tosparse_itself: " << tm_edgemap_tosparse_itself_ <<"\n"
              << " edgemap_tosparse_reduce: " << tm_edgemap_tosparse_reduce_ <<"\n"
              << " edgemap_todense: " << tm_edgemap_todense_ <<"\n"
              << " edgemap_dense_fwd: " << tm_edgemap_dense_fwd_ <<"\n"
              << " edgemap_dense_bwd: " << tm_edgemap_dense_bwd_ <<"\n"
              << " edgemap_dense_sync: " << tm_edgemap_dense_sync_ <<"\n"
              << " edgemap_sparse: " << tm_edgemap_sparse_ <<"\n"
              << " edgemap_sparse_sync_: " << tm_edgemap_sparse_sync_ <<"\n"
              <<std::endl;
}

typedef pair<intT, intT> intTpair;
//#define STEP (2048)
#define STEP (4096)
//*****START FRAMEWORK*****

//Design for partitioned_graph (dense, sparse)
//Dense contains a whole boolean array for each partitioned graph
//to share and update
//Sparse uses the vertices structure, m is the active vertex number
// in each partitioned
//intT * s is the intT array for update the out-degree array

class partitioned_vertices
{
public:
    intT numVertices;
    intT d_m;           // active vertices number for dense part
    bool *d;            // dense part, boolean array
    bool has_sparse;   // flag if sparse representation is present
    bool has_dense;    // flag if dense representation is present
    intT *s;            //sparse active vertices array
    intT num_out_edges; //acitve vertices's out-degree for sparse/dense selection

    //Create the initial partitioned frontier with one start vertex.
    //Used for traversal algorithm, such as BFS,BC and BellmanFord
    static partitioned_vertices create(intT n, intT v, intT initialOutdegree)
    {

        partitioned_vertices pv;
        pv.numVertices=n;
        pv.d = 0;
        pv.d_m = 1;
        pv.has_dense = false;
        pv.has_sparse = false;

        pv.s = new intT [1];
        pv.s[0] = v;
        pv.num_out_edges = initialOutdegree;
        return pv;
    }

    //Empty frontier to stop the programming
    static partitioned_vertices empty()
    {
        partitioned_vertices pv;
        pv.d = 0;
        pv.d_m = 0;
        pv.numVertices=0;
        pv.has_dense = false;
        pv.has_sparse = false;

        pv.s=NULL;
        pv.num_out_edges = 0;
        return pv;
    }

    //Create dense frontier for the dense part
    static partitioned_vertices dense(intT n, const partitioner & part)
    {
        partitioned_vertices pv;
        pv.numVertices=n;
        pv.d = new bool [n];
        pv.d_m = 0;
        pv.has_dense = true;
        pv.has_sparse = false;
        pv.s = NULL;
        pv.num_out_edges = 0;
        return pv;
    }
    //Create sparse frontier for the sparse part
    static partitioned_vertices sparse(intT n)
    {
        partitioned_vertices pv;
        pv.numVertices=n;
        pv.d = 0;
        pv.d_m = 0;
        pv.has_dense = false;
        pv.has_sparse = true;
        pv.s = NULL;
        pv.num_out_edges = 0;
        return pv;
    }


    //create of frontier with n from the boolean array as d[i]
    //vertices(intT _n,intT _m, bool* bits)
    //:n(_n),m(_m),s(NULL),isDense(1){}
    //and create the frontier boolean array inside the bits()
    //not to copy in each algorithm defined
    //Use for Components, BC, BellmanFord and PR, PRDelta

   static partitioned_vertices bits(const partitioner & part, intT n, intT initialOutdegree)
    {
        partitioned_vertices pv;
        pv.num_out_edges = initialOutdegree;
        pv.s = NULL;
        pv.numVertices=n;
        pv.d = new bool [n];
        parallel_for(intT i=0; i<n; i++) pv.d[i]=1;

        pv.d_m = n;
        pv.has_dense = true;
        pv.has_sparse = false;
        return pv;
    }
    //Used for vertexFilter to create the array with number of vertices and boolean array
    //This is using for the radii algorithm, and VertexFilter
    static partitioned_vertices boolean(intT n, bool* bits,intT activeM, intT out_edges)
    {
        partitioned_vertices pv;
        pv.numVertices=n;
        pv.num_out_edges = out_edges;
        pv.d_m = activeM;
        pv.s = NULL;
        pv.d = bits;
        pv.has_dense = true;
        pv.has_sparse = false;
        return pv;
    }


    void del()
    {
        if(s!=NULL) delete [] s;
        //if(s!=NULL) s.del();
        if( d )
            delete [] d;
            //d.del();
    }
    //get the intT array for sparse iteration
    intT* get_partition( unsigned p )
    {
        return s;
    }

    //check the empty boolean function for while loop
    bool isEmpty()
    {
        return d_m == 0;
    }

    //Represent the frontier as the dense format (boolean dense array [0,1,0,1,0,1,1,1,.....]), used for the edgeMapDense function
    void toDense()
    {
        // The dense pull operation requires to know for every vertex whether
        // it is active. So we convert all sparse partitions. The dense push
        // only needs to know for vertices within partitions traversed in
        // dense format. In those cases we could avoid creating the dense
        // representation, e.g., if the flag is set to use DENSE_FORWARD.
        if (!has_dense)
        {
            intT n = numVertices;
            if (!d) d = new bool [n];
            {
                parallel_for(intT i=0; i<n; i++) d[i] = 0;
            }
            {
                parallel_for(intT i=0; i<d_m; i++) d[s[i]] = 1;
            }
            has_dense = true;
        }
    }
    // Represent the frontier as the sparse format (only active vertex [1,5,8,..]), used for the edgeMapSparse Function
    void toSparse()
    {
        // Convert the partition p to sparse. Only one of them!
        if( has_dense && !has_sparse )
        {
            _seq<intT> R = sequence::packIndex(d,numVertices);
            if (d_m != R.n)
            {
                cout<<"bad stored value of m"<<endl;
                abort();
            }
            s=R.A;
            d_m = R.n;
            has_sparse = true;
            has_dense = false;
        }
    }

    // parttioned_vertices: numRows(); return n (all vertices number of frontier);
    intT numRows()
    {
        return numVertices;
    }
    // partitioned_vertices: return non active vertices
    intT numNonzeros()
    {
        return d_m;
    }
};

//To collect all outdegree for partitioned graph
template<class vertex>
class GoutDegree
{
public:
    partitioned_graph<vertex> &PG;
    bool* dense;

    GoutDegree(partitioned_graph<vertex>& pg, bool* pdense):PG(pg),dense(pdense) {}
    pair<intT, intT> operator()(intT i)
    {
        return make_pair((intT)dense[i], dense[i]? PG.getFullOutDegree(i):(intT)0);
    }

};
//To collect part outdegree for partitioned graph
template<class vertex>
class GoutDegreePart
{
public:
    partitioned_graph<vertex> &PG;
    bool* dense;
    intT p;

    GoutDegreePart(partitioned_graph<vertex>& pg, bool* pdense, intT p_)
        : PG(pg), dense(pdense), p(p_) {}
    pair<intT, intT> operator()(intT i)
    {
        bool inc = dense[i]; // && PG.get_partitioner().partition_of(i) == p;
        return make_pair((intT)inc, inc ? PG.get_partition(p).V[i].getOutDegree():(intT)0);
    }
};
//Get part out-degree of partitioned graph
template<class vertex>
class GoutDegreeV
{
public:
    graph<vertex> G;
    intT *s;

    GoutDegreeV(graph<vertex> G_, intT* s_) : G(G_), s(s_) { }
    intT operator()( intT i )
    {
        return G.V[s[i]].getOutDegree();
    }
};
//Get whole out-degree of partitioned graph
template<class vertex>
class GoutDegreeVFull
{
public:
    partitioned_graph<vertex> G;
    intT *s;

    GoutDegreeVFull(partitioned_graph<vertex> G_, intT* s_) : G(G_), s(s_) { }
    intT operator()( intT i )
    {
        return G.getFullOutDegree(s[i]);
    }
};

//Use for filter function for set the out-degree array without -1 for sparse
struct nonNegF
{
    bool operator() (intT a)
    {
        return (a>=0);
    }
};

//options to edgeMap for different versions of dense edgeMap (default is DENSE)
enum options { DENSE, DENSE_FORWARD};

//remove duplicate integers in [0,...,n-1]
void remDuplicates(intT* indices, intT* flags, intT m, intT n)
{
    //make flags for first time
    if(flags == NULL)
    {
        flags = new intT [n];
        {
            parallel_for(intT i=0; i<n; i++) flags[i]=-1;
        }
    }
    {
        parallel_for(intT i=0; i<m; i++)
        {
            if(indices[i] != -1 && flags[indices[i]] == -1)
                CAS(&flags[indices[i]],(intT)-1,i);
        }
    }
    //reset flags
    {
        parallel_for(intT i=0; i<m; i++)
        {
            if(indices[i] != -1)
            {
                if(flags[indices[i]] == i)  //win
                {
                    flags[indices[i]] = -1; //reset
                }
                else indices[i] = -1; //lost
            }
        }
    }
}

//*****EDGE FUNCTIONS*****
// Apply an edgeMapDense(Backward) to all vertices in the frontier.
// By construction of the graph, only vertices within the partition
// corresponding to GA will be set. Pull operator, walk all vertices,
// obtain values from active vertices following in-edges, update value
//and state (active) of vertex for next round. A pull operator writes
//the value of the active vertex by updating it with the neighbours' values.
//This is out-degree id check, shuttle different start vertex

// Note: pos could be a generic "edge index" for BP
//edgeOpFwd: Forward operator, traversal active source vertices (u) src->id
//update destination dst->ngh (V.getOutNeighbor(j))
template<class F>
inline bool edgeOpFwd( intT src, intT pos, intT dst, intE w, F f )
{
#ifdef MORE_ARG
    return f.cond(dst) && f.updateAtomic(src,dst,pos);
#else // MORE_ARG
#ifndef WEIGHTED
    return f.cond(dst) && f.updateAtomic(src,dst);
#else
    return f.cond(dst) && f.updateAtomic(src,dst,w);
#endif // WEIGHTED
#endif // MORE_ARG
}
//edgeOpBwd: Backward operator, src->ngh (destination) , dst->id (source)
template<class F>
inline bool edgeOpBwd( intT src, intT pos, intT dst, intE w, F f,
                       bool *frontier )
{
#ifdef MORE_ARG
    return frontier[src] && f.update(src,dst,pos);
#else // MORE_ARG
#ifndef WEIGHTED
    return frontier[src] && f.update(src,dst);
#else
    return frontier[src] && f.update(src,dst,w);
#endif // WEIGHTED
#endif // MORE_ARG
}

template<class F>
inline bool edgeOpBwdAtomic( intT src, intT pos, intT dst, intE w, F f,
                             bool *frontier )
{
#ifdef MORE_ARG
    return frontier[src] && f.updateAtomic(src,dst,pos);
#else // MORE_ARG
#ifndef WEIGHTED
    return frontier[src] && f.updateAtomic(src,dst);
#else
    return frontier[src] && f.updateAtomic(src,dst,w);
#endif // WEIGHTED
#endif // MORE_ARG
}

template<class F>
inline bool edgeOpIn( intT src, intT pos, intT dst, intE w, F f,
                      bool *frontier, bool *next )
{
#if REUSE_DIST
    rd_bwd_src->access(src);
    if( frontier[src] )
        rd_bwd_data->access(src);
#endif
    if( edgeOpBwd( src, pos, dst, w, f, frontier ) )
    {
#if REUSE_DIST
        rd_bwd_next->access(dst);
#endif
        next[dst]=1;
        return f.cond(dst);
    }
    return true;
}

template<class F>
inline bool edgeOpInAtomic( intT src, intT pos, intT dst, intE w, F f,
                            bool *frontier, bool *next )
{
    // std::cerr << src << " -> " << dst << '\n';

#if REUSE_DIST
    rd_bwd_src->access(src);
    if( frontier[src] )
        rd_bwd_data->access(src);
#endif
    if( edgeOpBwdAtomic( src, pos, dst, w, f, frontier ) )
    {
#if REUSE_DIST
        rd_bwd_next->access(dst);
#endif
        next[dst]=1;
        // return f.cond(dst);
    }
    return true;
}

template<class F>
inline bool edgeOpIn( intT src, typename F::cache_t & cache,
                      intT pos, intT dst, intE w, F f,
                      bool *frontier, bool *next )
{
#ifdef MORE_ARG
    if (frontier[src] && f.update(cache,src,pos))
    {
        next[dst]=1;
        return f.cond(dst);
    }
#else
#ifndef WEIGHTED
    if( frontier[src] && f.update(cache,src) )
    {
        next[dst]=1;
        return f.cond(dst);
    }
#else
    if( frontier[src] && f.update(cache,src,w) )
    {
        next[dst]=1;
        return f.cond(dst);
    }
#endif
#endif
    return true;
}

#if PULL_EDGES
#if HILBERT_EDGES || DEST_EDGES || 1
template<class F, class Edge>
bool* edgeMapDense(const EdgeList<Edge> & EL, bool* vertices, F f, bool *next,
                   bool parallel = false)
{
    typename EdgeList<Edge>::const_iterator E=EL.cend();
    parallel_for( typename EdgeList<Edge>::const_iterator
                  I=EL.cbegin(); I != E; ++I )
    {
        const Edge &eref = *I;

        intT src = eref.getSource();
        intT dst = eref.getDestination();
        intE wgh = eref.getWeight();
        /* When operating on partitions that have been created such that
         * all incoming edges to a node are located in the current partition,
         * then we can execute a non-atomic edge-op. */
        if( f.cond(dst) )
        {
            edgeOpInAtomic( src, -1, dst, wgh, f, vertices, next );
            //edgeOpIn( src, /*unused*/-1, dst, wgh, f, vertices, next );
        }
    }
    return next;
}
#else
template<class F, class Edge>
bool* edgeMapDense(const EdgeList<Edge> & EL, bool* vertices, F f, bool *next,
                   bool parallel = false)
{
    // This version of the code is valid when edgelist is sorted by destination
    static const size_t step_size = 2048;
    static const size_t m = EL.get_num_edges();
    parallel_for( size_t si=0; si < m; si += step_size )
    {
        size_t e = std::min( si + step_size, m );
        intT prevv = si == 0 ? ~intT(0) : EL[si-1].getDestination();
        intT nextv = e == m ? EL.get_num_vertices() : EL[e].getDestination();
        size_t i;
        for( i=si; i < e; ++i )
        {
            const Edge &eref = EL[i];

            intT src = eref.getSource();
            intT dst = eref.getDestination();
            intE wgh = eref.getWeight();

            if( f.cond(dst) )
            {
                if( dst != prevv )
                {
                    edgeOpIn( src, /*unused*/-1, dst, wgh, f, vertices, next );
                    ++i;
                    break;
                }
                else
                    edgeOpInAtomic( src, -1, dst, wgh, f, vertices, next );
            }
        }
        for( ; i < e; ++i )
        {
            const Edge &eref = EL[i];

            intT src = eref.getSource();
            intT dst = eref.getDestination();
            intE wgh = eref.getWeight();

            if( f.cond(dst) )
            {
                if( dst == nextv )
                {
                    edgeOpInAtomic( src, -1, dst, wgh, f, vertices, next );
                    ++i;
                    break;
                }
                else
                    edgeOpIn( src, /*unused*/-1, dst, wgh, f, vertices, next );
            }
        }
        for( ; i < e; ++i )
        {
            const Edge &eref = EL[i];

            intT src = eref.getSource();
            intT dst = eref.getDestination();
            intE wgh = eref.getWeight();

            if( f.cond(dst) )
                edgeOpInAtomic( src, -1, dst, wgh, f, vertices, next );
        }
    }
    return next;
}
#endif
#endif // PULL_EDGES

//If partitioning by source vertices, avoiding data race
//using the atomic update function
#if COMPRESSED_CSCCSR
template<class F, class vertex>
bool* edgeMapDenseCSC(graph<vertex> GA, bool * vertices,
                      F f, bool *next, bool source, bool parallel = 0)
{
    // timer dm;
    // dm.start();
    pair<intT,vertex> *G=GA.CSCV;
    intT CSCVn = GA.CSCVn;
    parallel_for (intT i=0; i<CSCVn; i++)
    {
        intT id = G[i].first;
        if (f.cond(id))
        {
            vertex V = G[i].second;
            intT d = V.getInDegree();
            if (source)
            {
                parallel_for(intT j=0; j<d; j++)
                {
                    intT ngh = V.getInNeighbor(j);
                    edgeOpInAtomic( ngh, j, id, V.getInWeight(j), f, vertices, next );
                }
            }
            else
            {
                if(!parallel || d < 1000)
                {
                    if(f.use_cache)
                    {
                        typename F::cache_t cache;
                        f.create_cache(cache,id);

                        for(intT j=0; j<d; j++)
                        {

                            intT ngh = V.getInNeighbor(j);
                            if( !edgeOpIn( ngh, cache, j, id, V.getInWeight(j), f, vertices, next ) )
                                break;
                        }

                        f.commit_cache(cache,id);
                    }
                    else     //use_cache ==false
                    {
                        for(intT j=0; j<d; j++)
                        {
                            intT ngh = V.getInNeighbor(j);
                            if( !edgeOpIn( ngh, j, id, V.getInWeight(j), f, vertices, next ) )
                                break;
                        }
                    }
                }
                else     //parallel&&d>1000
                {
                    parallel_for(intT j=0; j<d; j++)
                    {
                        intT ngh = V.getInNeighbor(j);
                        edgeOpInAtomic( ngh, j, id, V.getInWeight(j), f, vertices, next );
                    }
                }
            }
        }
    }
    // dm_dense_each=dm.next();
    // cerr << "MapDensePullCompressed: "<< dm_dense_each << endl;
    return next;
}
template <class F, class vertex>
bool* edgeMapDenseForwardCSR(graph<vertex> GA, bool *vertices, F f, bool* next)   // partitioned_vertices instead of vertices
{
    // timer fm;
    // fm.start();
    // intT end = GA.n_src_end;
    // intT start = GA.n_src_start;
    pair<intT,vertex> *G = GA.CSRV;
    intT CSRVn=GA.CSRVn;
// #define pragma cilk grainsize = 1
    parallel_for (intT i=0; i<CSRVn; i++)
    {
        intT id = G[i].first;
        if (vertices[id])
        {
            vertex V = G[i].second;
            intT d = V.getOutDegree();
            if(d < 1000)
            {
                for(intT j=0; j<d; j++)
                {
                    uintT ngh = V.getOutNeighbor(j);
#if REUSE_DIST
                    rd_fwd_dest->access(ngh);
#endif
                    if( edgeOpFwd( id, j, ngh, V.getOutWeight(j), f ) )
                    {
#if REUSE_DIST
                        rd_fwd_new->access(ngh);
#endif
                        next[ngh] = 1;
                    }
                }
            }
            else //d>1000
            {
                parallel_for(intT j=0; j<d; j++)
                {
                    uintT ngh = V.getOutNeighbor(j);
#if REUSE_DIST
                    rd_fwd_dest->access(ngh);
#endif
                    if( edgeOpFwd( id, j, ngh, V.getOutWeight(j), f ) )
                    {
#if REUSE_DIST
                        rd_fwd_new->access(ngh);
#endif
                        next[ngh] = 1;
                    }
                }
            }
        }
    }
    // fm_forward_each=fm.next();
    return next;
}
#else //COMPRESSED_CSCCSR
#if !COMPRESSED_VERTICES
template<class F, class vertex>
bool* edgeMapDense(graph<vertex> GA, bool * vertices,
                   F f, bool *next, bool source, bool parallel = 0)
{
//    timer dm;
//    dm.start();
    vertex *G=GA.V;
    intT start = GA.n_dst_start;
    intT end = GA.n_dst_end;
    // if ( source ) { traverse atomically only } else { traverse with non-atomic and cache opt; }
    parallel_for (intT id=start; id<end; id++)
    {
        vertex & V = G[id];
        if (f.cond(id))
        {
            intT d = V.getInDegree();
            if (source)
            {
                parallel_for(intT j=0; j<d; j++)
                {
                    intT ngh = V.getInNeighbor(j);
                    edgeOpInAtomic( ngh, j, id, V.getInWeight(j), f, vertices, next );
                }
            }
            else
            {
                if(!parallel || d < 1000)
                {
                    if(f.use_cache && d > 1)
                    {
                        typename F::cache_t cache;
                        f.create_cache(cache,id);

                        for(intT j=0; j<d; j++)
                        {
                            intT ngh = V.getInNeighbor(j);
                            if( !edgeOpIn( ngh, cache, j, id, V.getInWeight(j), f, vertices, next ) )
                                break;
                        }
                        f.commit_cache(cache,id);
                    }
                    else //use_cache ==false
                    {
                        for(intT j=0; j<d; j++)
                        {
                            intT ngh = V.getInNeighbor(j);
                            if( !edgeOpIn( ngh, j, id, V.getInWeight(j), f, vertices, next ) )
                                break;
                        }
                    }
                }
                else //parallel&&d>1000
                {
                    parallel_for(intT j=0; j<d; j++)
                    {
                        intT ngh = V.getInNeighbor(j);
                        edgeOpInAtomic( ngh, j, id, V.getInWeight(j), f, vertices, next );
                    }
                }
            }
        }
    }
    return next;
}

template<class F, class vertex>
bool* edgeMapDenseWithSkip(graph<vertex> GA, bool * vertices,
                           F f, bool *next, bool source, bool parallel = 0)
{
    timer dm;
    dm.start();
    vertex *G=GA.V;
    intT start = GA.n_dst_start;
    intT end = GA.n_dst_end;

    parallel_for (intT big_id=start; big_id<end; big_id+=STEP)
    {
        intT end_id = std::min(big_id+STEP, end);
        intT id = big_id;
        while( id < end_id )
        {
            vertex & V = G[id];
            if (f.cond(id))
            {
                intT d = V.getInDegree();
                if (source)
                {
                    parallel_for(intT j=0; j<d; j++)
                    {
                        intT ngh = V.getInNeighbor(j);
                        edgeOpInAtomic( ngh, j, id, V.getInWeight(j), f, vertices, next );
                    }
                }
                else
                {
                    if(!parallel || d < 1000)
                    {
                        if(f.use_cache && d > 1)
                        {
                            typename F::cache_t cache;
                            f.create_cache(cache,id);

                            for(intT j=0; j<d; j++)
                            {
                                intT ngh = V.getInNeighbor(j);
                                if( !edgeOpIn( ngh, cache, j, id, V.getInWeight(j), f, vertices, next ) )
                                    break;
                            }
                            f.commit_cache(cache,id);
                        }
                        else //use_cache ==false
                        {
                            for(intT j=0; j<d; j++)
                            {
                                intT ngh = V.getInNeighbor(j);
                                if( !edgeOpIn( ngh, j, id, V.getInWeight(j), f, vertices, next ) )
                                    break;
                            }
                        }
                    }
                    else //parallel&&d>1000
                    {
                        parallel_for(intT j=0; j<d; j++)
                        {
                            intT ngh = V.getInNeighbor(j);
                            edgeOpInAtomic( ngh, j, id, V.getInWeight(j), f, vertices, next );

                        }
                    }
                }
            }
            id = V.getInDegId();
        }
    }
    return next;

}

//Apply an edgemapDense(forward) to all vertices in the frontier.
//by construction of the graph, only vertices within the partition
//corresponding to GA will be set. Push operator, walk
//active vertex array, push updated value to all target vertices
//by following out-edges, and enables those vertices for the next round.
//Push operator reads the value of the active vertex and updates its neighbours.

template <class F, class vertex>
bool* edgeMapDenseForward(graph<vertex> GA, bool *vertices, F f, bool* next)
{
    // timer fm;
    // fm.start();
    vertex *G = GA.V;
    intT end = GA.n_src_end;
    intT start = GA.n_src_start;
    parallel_for (intT big_id=start; big_id<end; big_id+=STEP)
    {
        intT id=big_id;//initial the non-zero-out-degree vertex id
        intT max_id = min(big_id+STEP,end);
        while (id<max_id)
        {
            vertex V = G[id];
            if (vertices[id])
            {
                intT d = V.getOutDegree();
                if(d < 1000)
                {
                    for(intT j=0; j<d; j++)
                    {
                        uintT ngh = V.getOutNeighbor(j);
#if REUSE_DIST
                        rd_fwd_dest->access(ngh);
#endif
                        if( edgeOpFwd( id, j, ngh, V.getOutWeight(j), f ) )
                        {
#if REUSE_DIST
                            rd_fwd_new->access(ngh);
#endif
                            next[ngh] = 1;
                        }
                    }
                }
                else //d>1000
                {
                    parallel_for(intT j=0; j<d; j++)
                    {
                        uintT ngh = V.getOutNeighbor(j);
#if REUSE_DIST
                        rd_fwd_dest->access(ngh);
#endif
                        if( edgeOpFwd( id, j, ngh, V.getOutWeight(j), f ) )
                        {
#if REUSE_DIST
                            rd_fwd_new->access(ngh);
#endif
                            next[ngh] = 1;
                        }
                    }
                }
            }
            id = V.getOutDegId();
        }
    }

// fm_forward_each=fm.next();
    return next;
}

#else  //COMPRESSED_VERTICES
template<class F, class vertex>
bool* edgeMapDenseCompressed(graph<vertex> GA, bool * vertices,
                             F f, bool *next, bool source, bool parallel = 0)
{
    // timer dm;
    // dm.start();
    pair<intT,vertex> *G=GA.CV;
    intT CVn = GA.CVn;
    parallel_for (intT i=0; i<CVn; i++)
    {
        intT id = G[i].first;
        if (f.cond(id))
        {
            vertex V = G[i].second;
            intT d = V.getInDegree();
            if (source)
            {
                parallel_for(intT j=0; j<d; j++)
                {
                    intT ngh = V.getInNeighbor(j);
                    edgeOpInAtomic( ngh, j, id, V.getInWeight(j), f, vertices, next );
                }
            }
            else
            {
                if(!parallel || d < 1000)
                {
                    if(f.use_cache)
                    {
                        typename F::cache_t cache;
                        f.create_cache(cache,id);

                        for(intT j=0; j<d; j++)
                        {

                            intT ngh = V.getInNeighbor(j);
                            if( !edgeOpIn( ngh, cache, j, id, V.getInWeight(j), f, vertices, next ) )
                                break;
                        }

                        f.commit_cache(cache,id);
                    }
                    else     //use_cache ==false
                    {
                        for(intT j=0; j<d; j++)
                        {
                            intT ngh = V.getInNeighbor(j);
                            if( !edgeOpIn( ngh, j, id, V.getInWeight(j), f, vertices, next ) )
                                break;
                        }
                    }
                }
                else     //parallel&&d>1000
                {
                    parallel_for(intT j=0; j<d; j++)
                    {
                        intT ngh = V.getInNeighbor(j);
                        edgeOpInAtomic( ngh, j, id, V.getInWeight(j), f, vertices, next );
                    }
                }
            }
        }
    }
    // dm_dense_each=dm.next();
    // cerr << "MapDensePullCompressed: "<< dm_dense_each << endl;
    return next;
}

template <class F, class vertex>
bool* edgeMapDenseForwardCompressed(graph<vertex> GA, bool *vertices, F f, bool* next)   // partitioned_vertices instead of vertices
{
    // timer fm;
    // fm.start();
    // intT end = GA.n_src_end;
    // intT start = GA.n_src_start;
    pair<intT,vertex> *G = GA.CV;
    intT CVn=GA.CVn;
// #define pragma cilk grainsize = 1
    parallel_for (intT i=0; i<CVn; i++)
    {
        intT id = G[i].first;
        if (vertices[id])
        {
            vertex V = G[i].second;
            intT d = V.getOutDegree();
            if(d < 1000)
            {
                for(intT j=0; j<d; j++)
                {
                    uintT ngh = V.getOutNeighbor(j);
#if REUSE_DIST
                    rd_fwd_dest->access(ngh);
#endif
                    if( edgeOpFwd( id, j, ngh, V.getOutWeight(j), f ) )
                    {
#if REUSE_DIST
                        rd_fwd_new->access(ngh);
#endif
                        next[ngh] = 1;
                    }
                }
            }
            else //d>1000
            {
                parallel_for(intT j=0; j<d; j++)
                {
                    uintT ngh = V.getOutNeighbor(j);
#if REUSE_DIST
                    rd_fwd_dest->access(ngh);
#endif
                    if( edgeOpFwd( id, j, ngh, V.getOutWeight(j), f ) )
                    {
#if REUSE_DIST
                        rd_fwd_new->access(ngh);
#endif
                        next[ngh] = 1;
                    }
                }
            }
        }
    }
    // fm_forward_each=fm.next();
    return next;
}
#endif //COMPRESSED_VERTICES
#endif // COMPRESSED_CSCCSR
template <class F, class vertex>
pair<uintT,intT*> edgeMapSparseWithG(wholeGraph<vertex> GA, partitioned_vertices frontier, uintT Totalm, F f, intT remDups=0, intT* flags=NULL)
{
    //timer sm;
    // sm.start();
    vertex *V=GA.V;
    uintT *degrees = new uintT [Totalm];
    parallel_for (intT i = 0; i < Totalm; i++)
    degrees[i] = V[frontier.s[i]].getOutDegree();


    uintT* offsets = degrees;
    uintT outEdgeCount = sequence::plusScan(offsets, degrees, Totalm);
    intT* outEdges = new intT [outEdgeCount];

    parallel_for (intT k = 0; k < Totalm; k++)
    {
        intT v = frontier.s[k];
        intT o = offsets[k];
        vertex vert = V[v];
        intT d = vert.getOutDegree();
        if(d < 1000)
        {
            for (intT j=0; j < d; j++)
            {
                intT ngh = vert.getOutNeighbor(j);
                outEdges[o+j]
                    = edgeOpFwd( v, j, ngh, vert.getInWeight(j), f ) ? ngh : -1;
            }
        }
        else
        {
            parallel_for (intT j=0; j < d; j++)
            {
                intT ngh = vert.getOutNeighbor(j);
                outEdges[o+j]
                    = edgeOpFwd( v, j, ngh, vert.getInWeight(j), f ) ? ngh : -1;
            }
        }
    }
//Collect the active m of the localfrontier
    intT* nextIndices;
    uintT nextM;
    if(outEdgeCount==0)
    {
        nextIndices=0;
        nextM=0;
    }
    else
    {
        nextIndices = new intT [outEdgeCount];
        if(remDups) remDuplicates(outEdges,flags,outEdgeCount,remDups);
        // Filter out the empty slots (marked with -1)
        nextM = sequence::filter(outEdges,nextIndices,outEdgeCount,nonNegF());
    }
    delete [] outEdges;
    delete [] degrees;
    return pair<uintT,intT*>(nextM, nextIndices);
}


static int edgesTraversed = 0;

// decides on sparse or dense base on number of nonzeros in the active vertices
template <class F, class vertex>
partitioned_vertices edgeMap(wholeGraph<vertex> WG, partitioned_graph<vertex> GA, partitioned_vertices Localfrontier, F f, intT threshold = -1,
                             char option=DENSE, bool remDups=false)
{
    // timer tm;
    //tm.start();
    timer tm_setup;
    tm_setup.start();
    const partitioner &part = GA.get_partitioner();
    int partitions=part.get_num_partitions();
    // const partitioner &PVpart = Localfrontier.get_partitioner();
    // this should use data for all partitions
    intT numVertices = GA.n;
    if(threshold == -1) threshold = GA.m/20; //default threshold

    intT m = Localfrontier.numNonzeros();
    if (numVertices != Localfrontier.numRows())
    {
        cerr << "edgeMap: Sizes Don't match" << endl;
        abort();
    }

    intT TotalOutDegrees = Localfrontier.num_out_edges;
    edgesTraversed += TotalOutDegrees;
    //cerr << "m=" << m << " outdegs=" << TotalOutDegrees << endl;
    if(TotalOutDegrees == 0) return partitioned_vertices::empty();

    // m = NNZ, outDegrees is #out-edges for active vertices, threshold=edges in graph/20

#if REUSE_DIST
    rd_bwd_src->step();
    rd_bwd_data->step();
    rd_bwd_next->step();
    rd_fwd_dest->step();
    rd_fwd_new->step();
#endif

    // bool do_sparse[partitions];
    bool any_dense = false;
    bool any_sparse = false;
    partitioned_vertices v1;

    if(  m+TotalOutDegrees < threshold)
        // if( !Localfrontier.all_set && m+TotalOutDegrees < threshold)
    {
        any_sparse=true;
        Localfrontier.toSparse();
        //tm_edgemap_todense_ += tm.next();
    }
    else
    {
        any_dense = true;
        Localfrontier.toDense();
        //tm_edgemap_tosparse_itself_ += tm.next();
    }

    tmlog( tm_setup, tm_edgemap_setup_ );
    // cerr << "Setup: "<< tmlog( tm_setup, tm_edgemap_setup_ ) << endl;
    // Here try to remodify the order of graph traversal
    // The order is: if any_dense = false && any_sparse = true-> sparse only
    // if any_dense = true && any_sparse= false -> dense only
    // if any_dense = ture && any_sparse = true -> dense and sparse with boolean
    if (any_dense)
    {
        v1 = partitioned_vertices::dense(numVertices,part);
        parallel_for(intT i=0; i <numVertices; ++i)
        v1.d[i]=false;

        parallel_for( int p=0; p < partitions; ++p )
                {
                    if( GA.get_partition(p).m == 0 ) continue;

                    if( option == DENSE_FORWARD)
                    {
#if COMPRESSED_CSCCSR
                        if(GA.get_partition(p).CSRV)
                            edgeMapDenseForwardCSR(GA.get_partition(p), Localfrontier.d, f, v1.d);
#else
#if COMPRESSED_VERTICES
                        if(GA.get_partition(p).CV)
                            edgeMapDenseForwardCompressed(GA.get_partition(p), Localfrontier.d, f, v1.d);
#else
                        edgeMapDenseForward(GA.get_partition(p), Localfrontier.d, f, v1.d);
#endif // COMPRESSED_VERTICES
#endif // COMPRESSED_CSCCSR

                    }

                    else
                    {
#if COMPRESSED_CSCCSR
                        if(GA.get_partition(p).CSCV)
                            edgeMapDenseCSC(GA.get_partition(p),Localfrontier.d,f, v1.d, GA.source, (bool)option);
#else
#if COMPRESSED_VERTICES
                        if(GA.get_partition(p).CV)
                            edgeMapDenseCompressed(GA.get_partition(p),Localfrontier.d,f, v1.d, GA.source, (bool)option);
#else // COMPRESSED_VERTICES
#if PULL_EDGES
                        edgeMapDense(GA.get_edge_list_partition(p),Localfrontier.d, f, v1.d, (bool)option);
#else// PULL_EDGES
                        edgeMapDense/*WithSkip*/(GA.get_partition(p),Localfrontier.d,f, v1.d, GA.source, (bool)option);
#endif // PULL_EDGES
#endif // COMPRESSED_VERTICES
#endif // COMPRESSED_CSCCSR
                    }
                }

        tmlog( tm_setup, tm_edgemap_dense_ );
    }
    else    //sparse with sparse output
    {
        v1 = partitioned_vertices::sparse(numVertices);
        if( remDups )
        {
            pair<uintT,intT*> R
                = edgeMapSparseWithG( WG, Localfrontier, m, f,
                                      remDups, WG.flags );
            v1.s = R.second;
            v1.d_m = R.first;
        }
        else
        {
            pair<uintT,intT*> R
                = edgeMapSparseWithG( WG, Localfrontier, m, f );
            v1.s = R.second;
            v1.d_m = R.first;
        }

        tmlog( tm_setup, tm_edgemap_sparse_ );
        //cerr << "Map: "<< tmlog( tm_setup, tm_edgemap_map_ ) << endl;
    }

    // Heuristic: if any partition is dense, we do everything dense.
    // Collect the outdegree,active vertices set for next iteration
    if( any_dense )
    {
        // Calculate statistics on active vertices and their out-degree
        // Modify reduce function for partitioned not whole, but value is not correct,
        intTpair p = sequence::reduce<intT>((intT)0, (intT)GA.n,
                                            GoutDegree<vertex>(GA,v1.d));
        v1.d_m=p.first;
        v1.num_out_edges = p.second;
        v1.has_sparse = false;
    }
    else     // all partitions done sparsely
    {
        v1.has_dense = false;
        if (v1.d_m==0)
            v1.num_out_edges=0;
        else
            v1.num_out_edges = sequence::reduce<intT>((intT)0, v1.d_m, addF<intT>(), GoutDegreeVFull<vertex>(GA, v1.s));

    }
    //cerr << "Collect: "<< tmlog( tm_setup, tm_edgemap_collect_ ) << ' ' << ( any_dense ? "dense" : "sparse" ) << endl;
    tmlog( tm_setup, tm_edgemap_collect_ );
    return v1;
}

//*****VERTEX FUNCTIONS*****
//Note: this is the optimized version of vertexMap which does not
//perform a filter
template <class F>
void vertexMap(partitioned_vertices V, F add)
{

    // the vertexMap() is called before we call toDense/toSparse, so we
    // should only have one or the other representation, not both.
    if(V.has_dense)
    {
        parallel_for(intT i=0; i<V.numVertices; i++)
        if(V.d[i]) add(i);

    }
    else
    {
        parallel_for(intT i=0; i<V.d_m; i++)
        {
            add(V.s[i]);
        }
    }
}

//Note: this is the version of vertexMap in which only a subset of the
//input partitioned_vertices is returned
template <typename vertex, class F>
partitioned_vertices vertexFilter(partitioned_graph<vertex> GA, partitioned_vertices V, F filter)
{
//   const partitioner &part = V.get_partitioner();
    intT n = V.numRows();
    uintT m = V.numNonzeros();
    V.toDense();
    bool* d_out = new bool [n];
    {
        parallel_for(intT i=0; i<n; i++) d_out[i]=0;
    }
    {
        parallel_for(intT i=0; i<n; i++)
        if(V.d[i]) d_out[i] = filter(i);
    }
    intTpair p = sequence::reduce<intT>((intT)0, n, GoutDegree<vertex>(GA,d_out));
    intT activeM=p.first;
    intT out_edges=p.second;
    return partitioned_vertices::boolean(n,d_out,activeM,out_edges);
}


//cond function that always returns true
inline bool cond_true (intT d)
{
    return 1;
}

template<class GType, class GraphType>
void Compute( GType&, GraphType&, long);

//driver
int parallel_main(int argc, char* argv[])
{
    commandLine P(argc,argv," [-s] <inFile>");
    char* iFile = P.getArgument(0);
    bool symmetric = P.getOptionValue("-s");
    bool binary = P.getOptionValue("-b");
    long start = P.getOptionLongValue("-r",100);
    long rounds = P.getOptionLongValue("-rounds",1);
    int numOfNode = P.getOptionLongValue("-p", 4);
    char *part_how = P.getOptionValue("-P");
    bool part_src = true;
    if( !part_how || !strcmp( part_how, "dest" ) )
        part_src = false;
    else if( !strcmp( part_how, "source" ) )
        part_src = true;
    else
    {
        std::cerr << "Illegal value for -P: \"" << part_how
                  << "\". Allowed values: dest source. Default: dest\n";
        return 1;
    }

    if(symmetric)
    {
        wholeGraph<symmetricVertex> G =
            readGraph<symmetricVertex>(iFile,symmetric,binary); //symmetric graph
        //intT n=G.n;
        //intT m=G.m;
        partitioned_graph<symmetricVertex> PG( G, numOfNode, start, part_src);

#if REUSE_DIST
        rd_bwd_src = new RDLogger<intT>(n, 1);
        rd_bwd_data = new RDLogger<intT>(n, ~size_t(0));
        rd_bwd_next = new RDLogger<intT>(n, ~size_t(0));
        rd_fwd_dest = new RDLogger<intT>(n, 1);
        rd_fwd_new = new RDLogger<intT>(n, ~size_t(0));
#endif

        for(int r=0; r<rounds; r++)
        {

            //PAPI_initial();         /*PAPI Event inital*/
            //PAPI_start_count();   /*start PAPI counters*/
            startTime();
            Compute(G,PG,start);
            nextTime("Running time");
            //PAPI_stop_count();   /*stop PAPI counters*/
            // PAPI_print();   /* PAPI results print*/
            timePrint();    /* BreakTime Details print*/
        }
        PG.del();
        G.del();

    }
    else
    {
        wholeGraph<asymmetricVertex> G =
            readGraph<asymmetricVertex>(iFile,symmetric,binary); //asymmetric graph
        partitioned_graph<asymmetricVertex> PG( G, numOfNode, start, part_src);
        if(G.transposed) G.transpose();

#if REUSE_DIST
        rd_bwd_src = new RDLogger<intT>(n, 1);
        rd_bwd_data = new RDLogger<intT>(n, ~size_t(0));
        rd_bwd_next = new RDLogger<intT>(n, ~size_t(0));
        rd_fwd_dest = new RDLogger<intT>(n, 1);
        rd_fwd_new = new RDLogger<intT>(n, ~size_t(0));
#endif

            PAPI_initial();         /*PAPI Event inital*/
        for(int r=0; r<rounds; r++)
        {
            PAPI_start_count();   /*start PAPI counters*/
            startTime();
            Compute(G,PG,start);
            nextTime("Running time");
            PAPI_stop_count();   /*stop PAPI counters*/
            PAPI_print();   /* PAPI results print*/
        }
        PAPI_end();   /* PAPI results print*/
        PG.del();
        G.del();

    }

#if REUSE_DIST
    rd_bwd_src->dump( "rd_bwd_src" );
    rd_bwd_data->dump( "rd_bwd_data" );
    rd_bwd_next->dump( "rd_bwd_next" );
    rd_fwd_dest->dump( "rd_fwd_dest" );
    rd_fwd_new->dump( "rd_fwd_new" );
#endif

    // PAPI_print();   /* PAPI results print*/
    //timeprint();    /* Time Details print*/
    //timePrint();    /* BreakTime Details print*/
    return 0;
}

