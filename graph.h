// -*- C++ -*-
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include "parallel.h"
#include "partitioner.h"
#include <assert.h>
//#include <type_traits>
#include <unistd.h>
#include <sched.h>
#include <errno.h>
#include <cstring>
#include <string>
#include <utility>
#include <algorithm>

#include <sys/mman.h>

#ifndef PULL_EDGES
#define PULL_EDGES 0
#endif
#ifndef EDGES_HILBERT
#define EDGES_HILBERT 0
#endif
#ifndef PARTITION_RANGE
#define PARTITION_RANGE 1
#endif

#if 0
#include "cilkpub/sort.h"
template<typename It, typename Cmp>
void mysort( It begin, It end, Cmp cmp )
{
    cilkpub::cilk_sort( begin, end, cmp );
}
#else
template<typename It, typename Cmp>
void mysort( It begin, It end, Cmp cmp )
{
    std::sort( begin, end, cmp );
}
#endif

using namespace std;

//#define PAGESIZE (32)
//#define PAGESIZE (2048)
#define PAGESIZE (4096)

// To be used in boolean expressions
// **************************************************************
//    ADJACENCY ARRAY REPRESENTATION
// **************************************************************

class symmetricVertex
{
private:
    intE* neighbors;
    intT degree;
#if !COMPRESSED_VERTICES
    intT id;
#endif
public:
    void del()
    {
        delete [] neighbors;
    }
    symmetricVertex() {}
    symmetricVertex(intE* n, intT d
#if !COMPRESSED_VERTICES
                    , intT vertex_id
#endif
                   ) : neighbors(n), degree(d)
#if !COMPRESSED_VERTICES
        , id(vertex_id)
#endif
    {
    }

#ifndef WEIGHTED
    intE getInNeighbor(intT j)
    {
        return neighbors[j];
    }
    intE getOutNeighbor(intT j)
    {
        return neighbors[j];
    }
    intE getInWeight(intT j)
    {
        return 1;
    }
    intE getOutWeight(intT j)
    {
        return 1;
    }
#else
    //weights are stored in the entry after the neighbor ID
    //so size of neighbor list is twice the degree
    intE getInNeighbor(intT j)
    {
        return neighbors[2*j];
    }
    intE getOutNeighbor(intT j)
    {
        return neighbors[2*j];
    }
    intE getInWeight(intT j)
    {
        return neighbors[2*j+1];
    }
    intE getOutWeight(intT j)
    {
        return neighbors[2*j+1];
    }
#endif
    intT getInDegree()
    {
        return degree;
    }
    intT getOutDegree()
    {
        return degree;
    }
    void setInNeighbors(intE* _i)
    {
        neighbors = _i;
    }
    void setOutNeighbors(intE* _i)
    {
        neighbors=_i;
    }
    void setInDegree(intT _d)
    {
        degree = _d;
    }
    void setOutDegree(intT _d)
    {
        degree = _d;
    }
#if !COMPRESSED_VERTICES
    void setInDegId(intT _vertex_id)
    {
        id = _vertex_id;
    }
    void setOutDegId(intT _vertex_id)
    {
        id= _vertex_id;
    }
    intT getInDegId()
    {
        return id;
    }
    intT getOutDegId()
    {
        return id;
    }
#endif
    intE* getInNeighborPtr()
    {
        return neighbors;
    }
    intE* getOutNeighborPtr()
    {
        return neighbors;
    }
    void flipEdges() {}
};

class asymmetricVertex
{
private:
    intE* inNeighbors;
    intE* outNeighbors;
    intT outDegree;
    intT inDegree;
#if !COMPRESSED_VERTICES
    intT inDegId;
    intT outDegId;
#endif

public:
    intE* getInNeighborPtr()
    {
        return inNeighbors;
    }
    intE* getOutNeighborPtr()
    {
        return outNeighbors;
    }
    void del()
    {
        delete [] inNeighbors;
        delete [] outNeighbors;
    }
    asymmetricVertex() {}
    asymmetricVertex(intE* iN, intE* oN, intT id, intT od
#if !COMPRESSED_VERTICES
                     , intT iId, intT oId
#endif
                    ) : inNeighbors(iN), outNeighbors(oN), inDegree(id), outDegree(od)
#if !COMPRESSED_VERTICES
        , inDegId(iId), outDegId(oId)
#endif
    {}
#ifndef WEIGHTED
    intE getInNeighbor(intT j)
    {
        return inNeighbors[j];
    }
    intE getOutNeighbor(intT j)
    {
        return outNeighbors[j];
    }
    intE getInWeight(intT j)
    {
        return 1;
    }
    intE getOutWeight(intT j)
    {
        return 1;
    }
#else
    intE getInNeighbor(intT j)
    {
        return inNeighbors[2*j];
    }
    intE getOutNeighbor(intT j)
    {
        return outNeighbors[2*j];
    }
    intE getInWeight(intT j)
    {
        return inNeighbors[2*j+1];
    }
    intE getOutWeight(intT j)
    {
        return outNeighbors[2*j+1];
    }
#endif
    intT getInDegree()
    {
        return inDegree;
    }
    intT getOutDegree()
    {
        return outDegree;
    }
    void setInNeighbors(intE* _i)
    {
        inNeighbors = _i;
    }
    void setOutNeighbors(intE* _i)
    {
        outNeighbors = _i;
    }
    void setInDegree(intT _d)
    {
        inDegree = _d;
    }
    void setOutDegree(intT _d)
    {
        outDegree = _d;

    }
#if !COMPRESSED_VERTICES
    void setInDegId(intT _vertex_id)
    {
        inDegId = _vertex_id;
    }
    void setOutDegId(intT _vertex_id)
    {
        outDegId= _vertex_id;
    }
    intT getInDegId()
    {
        return inDegId;
    }
    intT getOutDegId()
    {
        return outDegId;
    }
#endif
    void flipEdges()
    {
        swap(inNeighbors,outNeighbors);
        swap(inDegree,outDegree);
#if !COMPRESSED_VERTICES
        swap(inDegId,outDegId);
#endif
    }
};

#if PULL_EDGES
class Edge
{
private:
    intT src, dst;
#ifdef WEIGHTED
    intE weight;
#endif

public:
    Edge() { }
#ifdef WEIGHTED
    Edge( intT s, intT d, intE w ) : src( s ), dst( d ), weight( w ) { }
#else
    Edge( intT s, intT d, intE w ) : src( s ), dst( d ) { }
    Edge( intT s, intT d ) : src( s ), dst( d ) { }
#endif
    Edge( const Edge & e ) : src( e.src ), dst( e.dst )
#ifdef WEIGHTED
        , weight( e.weight )
#endif
    { }

    intE getSource() const
    {
        return src;
    }
    intE getDestination() const
    {
        return dst;
    }

#ifdef WEIGHTED
    intT getWeight() const
    {
        return weight;
    }
#else
    intT getWeight() const
    {
        return 0;
    }
#endif
};

template<typename T>
typename std::make_unsigned<T>::type roundUpPow2( T n_u )
{
    const unsigned int num_bits = sizeof(T) * 8;

    typename std::make_unsigned<T>::type n = n_u;
    --n;
    for( unsigned int shift=1; shift < num_bits; shift <<= 1 )
        n |= n >> shift;
    ++n;

    return n;
}

// Source code based on https://en.wikipedia.org/wiki/Hilbert_curve
class HilbertEdgeSort
{
    intT n;

public:
    HilbertEdgeSort( intT n_ ) : n( roundUpPow2(n_) ) { }

    bool operator () ( const Edge & l, const Edge & r ) const
    {
        return e2d( l ) < e2d( r );
    }

private:
    intT e2d( const Edge & e ) const
    {
        return xy2d( e.getSource(), e.getDestination() );
    }
    //convert (x,y) to d
    intT xy2d( intT x, intT y ) const
    {
        intT rx, ry, s, d=0;
        for (s=n/2; s>0; s/=2)
        {
            rx = (x & s) > 0;
            ry = (y & s) > 0;
            d += s * s * ((3 * rx) ^ ry);
            rot(s, &x, &y, rx, ry);
        }
        return d;
    }

    //convert d to (x,y)
    void d2xy(intT d, intT *x, intT *y) const
    {
        intT rx, ry, s, t=d;
        *x = *y = 0;
        for (s=1; s<n; s*=2)
        {
            rx = 1 & (t/2);
            ry = 1 & (t ^ rx);
            rot(s, x, y, rx, ry);
            *x += s * rx;
            *y += s * ry;
            t /= 4;
        }
    }
    static void rot( intT n, intT *x, intT *y, intT rx, intT ry )
    {
        if (ry == 0)
        {
            if (rx == 1)
            {
                *x = n-1 - *x;
                *y = n-1 - *y;
            }
            std::swap( *x, *y );
        }
    }
};

template<class Edge>
class EdgeList
{
private:
    Edge * edges;
    intE num_edges;
    intT num_vertices;
public:
    EdgeList( intE m=0, intT n=0) : num_edges( m ) , num_vertices( n )
    {
        edges = m > 0 ? new Edge[m] : (Edge *)0;
    }
    EdgeList( const EdgeList<Edge> & el ) : num_edges( el.num_edges ), num_vertices(el.num_vertices)
    {
        edges = new Edge[num_edges];
        std::copy( &el.edges[0], &el.edges[num_edges], edges );
    }
    EdgeList( EdgeList<Edge> && el ) : num_edges( el.num_edges ) ,num_vertices(el.num_vertices)
    {
        edges = el.edges;
        el.num_edges = 0;
        el.edges = 0;
        el.num_vertices=0;
    }
    ~EdgeList()
    {
        if( edges )
            delete[] edges;
    }

    const EdgeList<Edge> & operator = ( const EdgeList<Edge> & el )
    {
        num_edges = el.num_edges;
        num_vertices = el.num_vertices;
        edges = new Edge[num_edges];
        std::copy( &el.edges[0], &el.edges[num_edges], edges );
        return *this;
    }
    const EdgeList<Edge> & operator = ( EdgeList<Edge> && el )
    {
        num_edges = el.num_edges;
        num_vertices = el.num_vertices;
        edges = el.edges;
        el.num_edges = 0;
        el.edges = 0;
        el.num_vertices = 0;
        return *this;
    }

    typedef Edge * iterator;
    typedef const Edge * const_iterator;

    iterator begin()
    {
        return &edges[0];
    }
    iterator end()
    {
        return &edges[num_edges];
    }
    const_iterator begin() const
    {
        return &edges[0];
    }
    const_iterator end() const
    {
        return &edges[num_edges];
    }
    const_iterator cbegin() const
    {
        return &edges[0];
    }
    const_iterator cend() const
    {
        return &edges[num_edges];
    }

    size_t get_num_edges() const
    {
        return num_edges;
    }
    size_t get_num_vertices() const
    {
        return num_vertices;
    }

    Edge & operator[] ( long i )
    {
        return edges[i];
    }
    const Edge & operator[] ( long i ) const
    {
        return edges[i];
    }

    void hilbert_sort()
    {
        mysort( begin(), end(), HilbertEdgeSort(num_vertices) );
    }
};
#endif
//This is the function for hugepage mmap allocation
//Based on the NUMA Awareness

template <class vertex>
class wholeGraph
{
public:
    vertex *V;  //mmap
    intT n;
    intT m;
    intE* allocatedInplace; //mmap
    intE* inEdges;           //mmap
    intT* flags;             //mmap
    bool transposed;
    bool isSymmetric;

    wholeGraph() {}
    wholeGraph(intT nn, intT mm, bool issym)
        : n(nn), m(mm), isSymmetric(issym),
          flags(NULL), transposed(false)
    {
        V = new vertex [n];
#ifndef WEIGHTED
        allocatedInplace = new intE [m];
#else//WEIGHTED
        allocatedInplace = new intE [2*m];
#endif
        if(!isSymmetric)
        {
#ifndef WEIGHTED
            inEdges = new intE [m];
#else//WEIGHTED
            inEdges = new intE [2*m];
#endif
        }
        else
            inEdges = NULL;
    }

    void del()
    {
        if (flags != NULL)
            delete [] flags;

        if(allocatedInplace!=NULL) delete [] allocatedInplace;
        delete [] V;
        if(inEdges != NULL) delete [] inEdges;

    }


    void transpose()
    {
        if(!isSymmetric)
        {
            parallel_for(intT i=0; i<n; i++)
            V[i].flipEdges();
            transposed = !transposed;
        }
    }
};

template <class vertex>
class graph
{
public:
    vertex *V;  //mmap
    intT n;
    intT m;
    intE* allocatedInplace; //mmap
    intE* inEdges;           //mmap
    intT* flags;             //mmap
    bool transposed;
    intT n_src_start, n_src_end, n_dst_start, n_dst_end;
    bool isSymmetric;
    int numanode;

#if COMPRESSED_CSCCSR
    pair<intT,vertex>* CSCV;    //mmap
    pair<intT,vertex>* CSRV;    //mmap
    intT CSCVn;
    intT CSRVn;
#endif
#if COMPRESSED_VERTICES
    pair<intT,vertex> * CV;    //mmap
    intT CVn;
#endif
    graph() {}
    graph(intT nn, intT mm, bool issym,int pp)
        : n(nn), m(mm), isSymmetric(issym),
          flags(NULL), transposed(false),
          n_src_start(0), n_src_end(nn),
          n_dst_start(0), n_dst_end(nn)
          ,numanode(pp)
    {
        V = 0;
#ifndef WEIGHTED
        allocatedInplace = new intE [m];
#else
        allocatedInplace = new intE [2*m];
#endif
        if(!isSymmetric)
        {
#ifndef WEIGHTED
            inEdges = new intE [m];
#else
            inEdges = new intE [2*m];
#endif
        }
        else
            inEdges = NULL;
#if COMPRESSED_CSCCSR
        CSCV=0;
        CSRV=0;
        CSCVn=0;
        CSRVn=0;
#endif
#if COMPRESSED_VERTICES
        CV = 0;
        CVn =0;
#endif
    }

    void del()
    {
        if (flags != NULL)
            delete [] flags;

        if(allocatedInplace!=NULL) delete [] allocatedInplace;
        if(inEdges != NULL) delete [] inEdges;
#if COMPRESSED_CSCCSR
        if(CSCV) delete [] CSCV;
        if(CSRV) delete [] CSRV;
#endif
#if COMPRESSED_VERTICES
        if(CV)
            delete [] CV;
#endif
    }

    void reorder_vertices( intT * __restrict reorder );
//    asymmetricVertex * newV = new asymmetricVertex, n);

    void transpose()
    {
        if(!isSymmetric)
        {
#if COMPRESSED_CSCCSR
            if(CSCV)
            {
                parallel_for (intT i=0; i<CSCVn; i++)
                CSCV[i].second.flipEdges();
            }
            if(CSRV)
            {
                parallel_for (intT i=0; i<CSRVn; i++)
                CSRV[i].second.flipEdges();
            }
#endif
#if COMPRESSED_VERTICES
            if(CV)
            {
                parallel_for (intT i=0; i<CVn; i++)
                CV[i].second.flipEdges();
            }
#endif
            transposed = !transposed;
        }
    }

#if PULL_EDGES
#if DEST_EDGES
    EdgeList<Edge> build_edgelist() const
    {
        EdgeList<Edge> el (m,n);
        long k = 0;
        for( intT i=0; i<n; i++ )
        {
            vertex v = V[i];
            for( intT j=0; j < v.getOutDegree(); ++j )
            {
                intT d = v.getOutNeighbor( j );
#ifndef WEIGHTED
                el[k++] = Edge( i, d );
#else
                el[k++] = Edge( i, d, v.getOutWeight( j ) );
#endif
            }
        }
        assert( k == m );
        return el;
    }
#else//DEST_EDGES
    EdgeList<Edge> build_edgelist() const
    {
        EdgeList<Edge> el (m,n);
        long k = 0;
        for( intT i=0; i<n; i++ )
        {
            vertex v = V[i];
            for( intT j=0; j < v.getInDegree(); ++j )
            {
                intT d = v.getInNeighbor( j );
#ifndef WEIGHTED
                el[k++] = Edge( d, i );
#else
                el[k++] = Edge( d, i, v.getInWeight( j ) );
#endif
            }
        }
        assert( k == m );
        return el;
    }
#endif//DEST_EDGES
#endif//PULL_EDGES
};


template <class vertex>
class partitioned_graph
{
public:
    typedef vertex vertex_type;
    //typedef Allocator  allocator_type;
    vertex *Ver;
    partitioner partition;
    graph<vertex> * localG;
    intT m,n;
    bool source;
#if PULL_EDGES
    EdgeList<Edge> * localEdgeList;
#endif
    partitioned_graph( wholeGraph<vertex> & GA, int num_part, long start,
                       bool partition_source )
        : partition(num_part,GA.n),m(GA.m),n(GA.n),source(partition_source)
    {

        localG = new graph<vertex> [num_part];

#if PULL_EDGES
        localEdgeList = new EdgeList<Edge>[num_part];
#endif

        Ver = GA.V;
        // partition.set_num_elements(GA.n);

        //std::cerr << "partition graph: n=" << GA.n << " m=" << GA.m << std::endl;

#if PARTITION_RANGE
        // Normal path
        partitionByDegree( GA, num_part, partition.as_array(), sizeof(intT), partition_source );
        partition.compute_starts();

#if MIX
        {
            partitioner by2(2,GA.n);
            partitionByDegree( GA, 2, by2.as_array(), sizeof(intT), false );
            by2.compute_starts();
            graph<vertex> tempGA[2] =
            {
                PartitionByDest( GA, by2.start_of(0), by2.start_of(1), 0 ),
                PartitionByDest( GA, by2.start_of(1), by2.start_of(2), 2 )
            };
            partitionByDegree( tempGA[0], 2, partition.as_array(), sizeof(intT), true );
            partitionByDegree( tempGA[1], 2, &partition.as_array()[2], sizeof(intT), true );
            //partition.compute_starts();
            parallel_for_numa( int p = 0; p < 4; ++p )
            {
                intT src = p < 2 ? 0 : 1;
                intT dst = p & 1;
                // localG[p] = PartitionBySour( tempGA[src], partition.start_of(p), partition.start_of(p+1), p );
                intT start = ((p&1)==0) ? 0 : partition.as_array()[p-1];
                intT end = ((p&1)==0) ? partition.as_array()[p] : GA.n;
                localG[p] = PartitionBySour( tempGA[src], start, end, p );
                localG[p].n_dst_start = tempGA[src].n_dst_start;
                localG[p].n_dst_end = tempGA[src].n_dst_end;
            }
            tempGA[0].del();
            tempGA[1].del();
            partition.as_array()[0] = GA.n/4;
            partition.as_array()[1] = GA.n/4;
            partition.as_array()[2] = GA.n/4;
            partition.as_array()[3] = GA.n-3*(GA.n/4);
            partition.compute_starts();
        }

#else //MIX
        parallel_for( int p=0; p < partition.get_num_partitions(); ++p )
        {
            if( partition_source )
                localG[p] = PartitionBySour( GA, partition.start_of(p), partition.start_of(p+1),p);
            else
                localG[p] = PartitionByDest( GA, partition.start_of(p), partition.start_of(p+1),p);
#if PULL_EDGES
            localEdgeList[p] = localG[p].build_edgelist();
#if EDGES_HILBERT
            localEdgeList[p].hilbert_sort();
#endif
#endif
        }
#endif // MIX
#else // PARTITION_RANGE
        partitionApprox( GA, partition.as_array(), num_part );
        parallel_for( int p=0; p < partition.get_num_partitions(); ++p )
        {
            localG[p] = graphFilter( GA, partition.as_array(), p );
            intT s = sequence::reduce<intT>((intT)0,(intT)GA.n,addF<intT>(),IsPart(partition,p));
            partition.set_size( p, s );
        }
#endif
    }
    void del()
    {
        for( int p=0; p < partition.get_num_partitions(); ++p )
        {
            localG[p].del();

        }
        delete [] localG;
    }

    // Translate vertex id to partition
    int partition_of( intT vertex_id )
    {
        return partition.partition_of( vertex_id );
    }

    // Get a partition of the graph to operate on
    graph<vertex> get_partition( int p )
    {
        return localG[p];
    }
#if PULL_EDGES
    const EdgeList<Edge> & get_edge_list_partition( intT p )
    {
        return localEdgeList[p];
    }
#endif
    int get_num_partitions() const
    {
        return partition.get_num_partitions();
    }
    // Get the partitioner object constructed for this graph
    const partitioner & get_partitioner() const
    {
        return partition;
    }

    intT getFullOutDegree(intT v) const
    {
        return Ver[v].getOutDegree();
    }

    void transpose()
    {
        parallel_for(unsigned i=0; i<partition.get_num_partitions(); i++)
        {
            localG[i].transpose();
        }
    }

private:
// These functions are specific to the partitioned graph. No one needs
// to know about them, I think.
// graphFilter is polymer method, partitioend graph:
// Pull: each local vertex contains its own out-degree
// Push: each local vertex contains its own in-degree
// PartitionByDest :
// Pull: each local vertex contains its own in-degree
// Push each local vertex contains its own in-degree
// PartitionBySour:
// Pull: each local vertex contains its own out-degree
// Push: each local vertex contains its own out-degree
    void partitionByDegree( wholeGraph<vertex> GA, int numOfNode, intT *sizeArr,
                            int sizeOfOneEle, bool useOutDegree=false );
    void partitionApprox(wholeGraph<vertex> & GA, short * partitions, int np);
    graph<vertex> graphFilter( wholeGraph<vertex>& GA, int rangeLow, int rangeHi ,int numanode);
    graph<vertex> graphFilter( wholeGraph<vertex>& GA, short * partition, intT p );
    graph<vertex> PartitionByDest(wholeGraph<vertex>& GA, int rangeLow, int rangeHi ,int numanode);
    graph<vertex> PartitionBySour(wholeGraph<vertex>& GA, int rangeLow, int rangeHi ,int numanode);
};


// ======================================================================
// Graph Filtering (Graph Partitioned)
// ======================================================================
//This method is partitioned by outdegree, CSR, store the all out-degree for
//all local vertices of each partitions.
//And part indegree for whole vertices
// Place edge (u,v) in partition of u
template <class vertex>
graph<vertex> partitioned_graph<vertex>::PartitionBySour(wholeGraph<vertex> &GA, int rangeLow, int rangeHi,int numanode)
{
    //vertex *V = GA.getV();
    vertex *V = GA.V;
    const intT n = GA.n;
    bool isSymmetric = GA.isSymmetric;

    intT *counters = new intT [n];
    intT *offsets = new intT [n];
    intT *inCounters = new intT [n];
    intT *inOffsets = new intT [n];
    {
        parallel_for (intT i = 0; i < n; i++)
        {
            inCounters[i] = 0;
            counters[i] = 0;
            if(!isSymmetric)
            {
                intT d = V[i].getInDegree();
                for (intT j = 0; j < d; j++)
                {
                    intT ngh = V[i].getInNeighbor(j);
                    if (rangeLow <= ngh && ngh < rangeHi)
                        inCounters[i]++;
                }
            }
            if(rangeLow<=i && i<rangeHi)
            {
                intT ind = V[i].getOutDegree();
                counters[i] += ind;
            }
        }
    }
    intT totalSize = 0;
    intT totalInSize = 0;
    for (intT i = 0; i < n; i++)
    {
        offsets[i] = totalSize;
        totalSize += counters[i];
        if(!isSymmetric)
        {
            inOffsets[i] = totalInSize;
            totalInSize += inCounters[i];
        }
    }

    if(!isSymmetric)
        assert( totalSize == totalInSize );

    graph<vertex> FG(n, totalSize, isSymmetric, numanode);
    FG.V = new vertex [n];

    {
        parallel_for (intT i = 0; i < n; i++)
        {
            FG.V[i].setOutDegree(counters[i]);
            if(!isSymmetric)
                FG.V[i].setInDegree(inCounters[i]);
        }
    }

    delete [] counters;
    delete [] inCounters;

    intE *edges = FG.allocatedInplace;
    intE *inEdges = FG.inEdges;//0;

    {
        parallel_for (intT i = 0; i < n; i++)
        {
#ifndef WEIGHTED
            intE *localInEdges = &inEdges[inOffsets[i]];
#else
            intE *localInEdges = &inEdges[inOffsets[i]*2];
#endif
            intT incounter = 0;
            if(!isSymmetric)
            {
                intT d = V[i].getInDegree();
                for (intT j = 0; j < d; j++)
                {
#ifndef WEIGHTED
                    intT ngh = V[i].getInNeighbor(j);
#else
                    intT ngh = V[i].getInNeighbor(j);
                    intT wgh = V[i].getInWeight(j);
#endif
                    if (rangeLow <= ngh && ngh < rangeHi)
                    {
#ifndef WEIGHTED
                        localInEdges[incounter] = ngh;
#else
                        localInEdges[incounter*2]= ngh;
                        localInEdges[incounter*2+1]= wgh;
#endif

                        incounter++;
                    }
                }
                FG.V[i].setInNeighbors(localInEdges);
            }
#ifndef WEIGHTED
            intE *localEdges = &edges[offsets[i]];
#else
            intE *localEdges = &edges[offsets[i]*2];
#endif
            intT counter = 0;
            if (rangeLow <= i && i < rangeHi)
            {
                intT ind = V[i].getOutDegree();
                for (intT j = 0; j < ind; j++)
                {
#ifndef WEIGHTED
                    intT ngh = V[i].getOutNeighbor(j);
                    localEdges[counter] = ngh;
#else
                    intT ngh = V[i].getOutNeighbor(j);
                    intT wgh = V[i].getOutWeight(j);
                    localEdges[counter*2] = ngh;
                    localEdges[counter*2+1] = wgh;
#endif

                    counter++;
                }
                FG.V[i].setOutNeighbors(localEdges);
            }
        }
    }
#if COMPRESSED_CSCCSR
    intT nnz = 0, nnzi = 0, nnzo = 0, ne = 0;
    for( intT i=0; i < n; ++i )
    {
        if( FG.V[i].getInDegree() != 0 || FG.V[i].getOutDegree() != 0 )
        {
            ++nnz;
            if( FG.V[i].getInDegree() != 0 )
                ++nnzi;
            if( FG.V[i].getOutDegree() != 0 )
                ++nnzo;
        }
        ne += FG.V[i].getOutDegree();
    }
    FG.CSCVn=nnzi;
    FG.CSRVn=nnzo;

    FG.CSCV = new pair<intT,vertex> [FG.CSCVn];
    FG.CSRV = new pair<intT,vertex> [FG.CSRVn];

    intT k=0;
    for( intT i=0; i < n; ++i )
    {
        if( FG.V[i].getInDegree() != 0 )
        {
            FG.CSCV[k].first = i;
            FG.CSCV[k].second = FG.V[i];
            k++;
        }
    }
    intT j=0;
    for( intT i=0; i < n; ++i )
    {
        if( FG.V[i].getOutDegree() != 0 )
        {
            FG.CSRV[j].first = i;
            FG.CSRV[j].second = FG.V[i];
            j++;
        }
    }
    cerr << "CSCCSRCompressed graph n=" << n << " Compressed=" << nnz
         << " nnzi=" << nnzi << " nnzo=" << nnzo << " out-edges=" << ne << " compressed" << endl;

#endif // COMPRESSED_CSCCSR
#if COMPRESSED_VERTICES
    intT nnz = 0, nnzi = 0, nnzo = 0, ne = 0;
    for( intT i=0; i < n; ++i )
    {
        if( FG.V[i].getInDegree() != 0 || FG.V[i].getOutDegree() != 0 )
        {
            ++nnz;
            if( FG.V[i].getInDegree() != 0 )
                ++nnzi;
            if( FG.V[i].getOutDegree() != 0 )
                ++nnzo;
        }
        ne += FG.V[i].getOutDegree();
    }
    if( nnzo < n/2 || true )   // HEURISTIC TO USE THESE IN DENSE TRAVERSAL
    {
        FG.CVn = nnz;

        FG.CV = new pair<intT,vertex> [FG.CVn];

        intT k=0;
        for( intT i=0; i < n; ++i )
            if( FG.V[i].getInDegree() != 0 || FG.V[i].getOutDegree() != 0 )
            {
                FG.CV[k].first = i;
                FG.CV[k].second = FG.V[i];
                k++;
            }
        cerr << "Compressed graph n=" << n << " Compressed=" << FG.CVn
             << " nnzi=" << nnzi << " nnzo=" << nnzo << " out-edges=" << ne << " compressed" << endl;
    }
    else
    {
        cerr << "Compressed graph n=" << n << " nonzeroes=" << nnz
             << " nnzi=" << nnzi << " nnzo=" << nnzo << " out-edges=" << ne << " not built" << endl;
    }
#else
    //For push operator, empty vertex check of outdegree
    for( intT i=n-1, next_id=n; i >= 0; --i )
    {
        FG.V[i].setOutDegId(next_id);
        if( FG.V[i].getOutDegree() != 0 )
            next_id = i;
    }

    // Remember next vertex with non-zero indegree
    for( intT i=n-1, next_id=n; i >= 0; --i )
    {
        FG.V[i].setInDegId(next_id);
        if( FG.V[i].getInDegree() != 0 )
            next_id = i;
    }
#endif

    delete [] offsets;
    delete [] inOffsets;

    // We know that all vertices with outgoing edges are located within this
    // partition, so when iterating over sources of edges, we can limit the
    // range to:
    FG.n_src_start = rangeLow;
    FG.n_src_end = rangeHi;
    delete [] FG.V;
    return FG;//return graph<vertex>(newVertexSet, GA.n, GA.m);
}

//This method is designed for indegree, and out model use now
//all the vertices of graph will store part outdegree (CSR) for push
//local vertices of partitioend graph will store all in-degree (CSC) for pull
// Place edge (u,v) in partition of v
template <class vertex>
graph<vertex> partitioned_graph<vertex>::PartitionByDest(wholeGraph<vertex> &GA, int rangeLow, int rangeHi,int numanode)
{
    //vertex *V = GA.getV();
    vertex *V = GA.V;
    const intT n = GA.n;
    bool isSymmetric = GA.isSymmetric;
    intT *counters = new intT [n];
    intT *offsets = new intT [n];
    intT *inCounters = new intT [n];
    intT *inOffsets = new intT [n];
    {
        parallel_for (intT i = 0; i < n; i++)
        {
            counters[i] = 0;
            inCounters[i] = 0;
            intT d = V[i].getOutDegree();
            for (intT j = 0; j < d; j++)
            {
                intT ngh = V[i].getOutNeighbor(j);
                if (rangeLow <= ngh && ngh < rangeHi)
                    counters[i]++;
            }
            if(!isSymmetric)
            {
                if(rangeLow <= i && i < rangeHi)
                    inCounters[i] = V[i].getInDegree();
            }
        }
    }
    intT totalSize = 0;
    intT totalInSize = 0;
    for (intT i = 0; i < n; i++)
    {
        offsets[i] = totalSize;
        totalSize += counters[i];

        if(!isSymmetric)
        {
            inOffsets[i] = totalInSize;
            totalInSize += inCounters[i];
        }
    }

    //cerr<<"totalSize: "<<totalSize<<" and totalInSize: "<<totalInSize<<endl;
    if(!isSymmetric)
        assert( totalSize == totalInSize );

    graph<vertex> FG(n, totalSize, isSymmetric, numanode);
    FG.V = new vertex [n];
    {
      //  vertex * FGV = FG.getV();
        parallel_for (intT i = 0; i < n; i++)
        {
            FG.V[i].setOutDegree(counters[i]);
            if(!isSymmetric)
                FG.V[i].setInDegree(inCounters[i]);
        }
    }
    delete [] counters;
    delete [] inCounters;

    intE *edges = FG.allocatedInplace;
    intE *inEdges = FG.inEdges;

    {
        parallel_for (intT i = 0; i < n; i++)
        {
#ifndef WEIGHTED
            intE *localEdges = &edges[offsets[i]];
#else
            intE *localEdges = &edges[offsets[i]*2];
#endif
            intT counter = 0;
            intT d = V[i].getOutDegree();
            for (intT j = 0; j < d; j++)
            {
#ifndef WEIGHTED
                intT ngh = V[i].getOutNeighbor(j);
#else
                intT ngh = V[i].getOutNeighbor(j);
                intT wgh = V[i].getOutWeight(j);
#endif
                if (rangeLow <= ngh && ngh < rangeHi)
                {
#ifndef WEIGHTED
                    localEdges[counter] = ngh;
#else
                    localEdges[counter*2]= ngh;
                    localEdges[counter*2+1]= wgh;
#endif
                    counter++;
                }
            }

            FG.V[i].setOutNeighbors(localEdges);
            if(!isSymmetric)
            {
#ifndef WEIGHTED
                intE *localInEdges = &inEdges[inOffsets[i]];
#else
                intE *localInEdges = &inEdges[inOffsets[i]*2];
#endif
                intT incounter = 0;
                if (rangeLow <= i && i < rangeHi)
                {
                    d = V[i].getInDegree();
                    for (intT j = 0; j < d; j++)
                    {
#ifndef WEIGHTED
                        intT ngh = V[i].getInNeighbor(j);
                        localInEdges[incounter] = ngh;
#else
                        intT ngh = V[i].getInNeighbor(j);
                        intT wgh = V[i].getInWeight(j);

                        localInEdges[incounter*2] = ngh;
                        localInEdges[incounter*2+1] = wgh;
#endif

                        incounter++;
                    }
                    FG.V[i].setInNeighbors(localInEdges);
                }
            }
        }
    }
#if COMPRESSED_CSCCSR
    intT nnz = 0, nnzi = 0, nnzo = 0, ne = 0;
    for( intT i=0; i < n; ++i )
    {
        if( FG.V[i].getInDegree() != 0 || FG.V[i].getOutDegree() != 0 )
        {
            ++nnz;
            if( FG.V[i].getInDegree() != 0 )
                ++nnzi;
            if( FG.V[i].getOutDegree() != 0 )
                ++nnzo;
        }
        ne += FG.V[i].getOutDegree();
    }
    FG.CSCVn=nnzi;
    FG.CSRVn=nnzo;

    FG.CSCV = new pair<intT,vertex> [FG.CSCVn];
    FG.CSRV = new pair<intT,vertex> [FG.CSRVn];

    intT k=0;
    for( intT i=0; i < n; ++i )
    {
        if( FG.V[i].getInDegree() != 0 )
        {
            FG.CSCV[k].first = i;
            FG.CSCV[k].second = FG.V[i];
            k++;
        }
    }
    intT j=0;
    for( intT i=0; i < n; ++i )
    {
        if( FG.V[i].getOutDegree() != 0 )
        {
            FG.CSRV[j].first = i;
            FG.CSRV[j].second = FG.V[i];
            j++;
        }
    }
    cerr << "CSCCSRCompressed graph n=" << n << " Compressed=" << nnz
         << " nnzi=" << nnzi << " nnzo=" << nnzo << " out-edges=" << ne << " compressed" << endl;

#endif // COMPRESSED_CSCCSR
#if COMPRESSED_VERTICES
    intT nnz = 0, nnzi = 0, nnzo = 0, ne = 0;
    for( intT i=0; i < n; ++i )
    {
        if( FG.V[i].getInDegree() != 0 || FG.V[i].getOutDegree() != 0 )
        {
            ++nnz;
            if( FG.V[i].getInDegree() != 0 )
                ++nnzi;
            if( FG.V[i].getOutDegree() != 0 )
                ++nnzo;
        }
        ne += FG.V[i].getOutDegree();
    }
    if( nnzo < n/2 || true )   // HEURISTIC TO USE THESE IN DENSE TRAVERSAL
    {
        FG.CVn = nnz;

        FG.CV = new pair<intT,vertex> [FG.CVn];

        intT k=0;
        for( intT i=0; i < n; ++i )
            if( FG.V[i].getInDegree() != 0 || FG.V[i].getOutDegree() != 0 )
            {
                FG.CV[k].first = i;
                FG.CV[k].second = FG.V[i];
                k++;
            }
        cerr << "Compressed graph n=" << n << " Compressed=" << FG.CVn
             << " nnzi=" << nnzi << " nnzo=" << nnzo << " out-edges=" << ne << " compressed" << endl;
    }
    else
    {
        cerr << "Compressed graph n=" << n << " nonzeroes=" << nnz
             << " nnzi=" << nnzi << " nnzo=" << nnzo << " out-edges=" << ne << " not built" << endl;
    }
#else
    //For push operator, empty vertex check of outdegree
    for( intT i=n-1, next_id=n; i >= 0; --i )
    {
        FG.V[i].setOutDegId(next_id);
        if( FG.V[i].getOutDegree() != 0 )
            next_id = i;
    }

    // Remember next vertex with non-zero indegree
    for( intT i=n-1, next_id=n; i >= 0; --i )
    {
        FG.V[i].setInDegId(next_id);
        if( FG.V[i].getInDegree() != 0 )
            next_id = i;
    }
#endif

    delete [] offsets;
    delete [] inOffsets;
    // We know that all vertices with incoming edges are located within this
    // partition, so when iterating over destinations of edges, we can limit the
    // range to:
    FG.n_dst_start = rangeLow;
    FG.n_dst_end = rangeHi;
    delete [] FG.V;
    return FG;
}

//graph partition by each algorithm
template <class vertex>
graph<vertex> partitioned_graph<vertex>::graphFilter(wholeGraph<vertex> &GA, int rangeLow, int rangeHi,int numanode)
{
    vertex *V = GA.V;
    const intT n = GA.n;
    bool isSymmetric = GA.isSymmetric;
    intT *counters = new intT [n];
    intT *offsets = new intT [n];
    intT *inCounters = new intT [n];
    intT *inOffsets = new intT [n];
    {
        parallel_for (intT i = 0; i < n; i++)
        {
            intT d = V[i].getOutDegree();
            counters[i] = 0;
            inCounters[i] = 0;
            for (intT j = 0; j < d; j++)
            {
                intT ngh = V[i].getOutNeighbor(j);
                if (rangeLow <= ngh && ngh < rangeHi)
                    counters[i]++;
            }

            // This is redundant for symmetricVertex.
            if(!isSymmetric)
            {
                intT ind = V[i].getInDegree();
                for (intT j = 0; j < ind; j++)
                {
                    intT ngh = V[i].getInNeighbor(j);
                    if (rangeLow <= ngh && ngh < rangeHi)
                        inCounters[i]++;
                }
            }
        }
    }

    intT totalSize = 0;
    intT totalInSize = 0;
    for (intT i = 0; i < n; i++)
    {
        offsets[i] = totalSize;
        totalSize += counters[i];

        if(!isSymmetric)
        {
            inOffsets[i] = totalInSize;
            totalInSize += inCounters[i];
        }
    }

    graph<vertex> FG(n, totalSize, isSymmetric, numanode);
    FG.V = new vertex [n];

    {
        parallel_for (intT i = 0; i < n; i++)
        {
            FG.V[i].setOutDegree(counters[i]);
            if(!isSymmetric)
                FG.V[i].setInDegree(inCounters[i]);
        }
    }

    delete [] counters;
    delete [] inCounters;

    intE *edges = FG.allocatedInplace;
    intE *inEdges = FG.inEdges;

    // std::cerr << "(out)edges in partition: " << totalSize << std::endl;

    {
        parallel_for (intT i = 0; i < n; i++)
        {
#ifndef WEIGHTED
            intE *localEdges = &edges[offsets[i]];
#else
            intE *localEdges = &edges[offsets[i]*2];
#endif
            intT counter = 0;
            intT d = V[i].getOutDegree();
            for (intT j = 0; j < d; j++)
            {
#ifndef WEIGHTED
                intT ngh = V[i].getOutNeighbor(j);
#else
                intT ngh = V[i].getOutNeighbor(j);
                intT wgh = V[i].getOutWeight(j);
#endif
                if (rangeLow <= ngh && ngh < rangeHi)
                {
#ifndef WEIGHTED
                    localEdges[counter] = ngh;
#else
                    localEdges[counter*2]= ngh;
                    localEdges[counter*2+1]= wgh;
#endif
                    counter++;
                }
            }
            FG.V[i].setOutNeighbors(localEdges);
            if(!isSymmetric)
            {
#ifndef WEIGHTED
                intE *localInEdges = &inEdges[inOffsets[i]];
#else
                intE *localInEdges = &inEdges[inOffsets[i]*2];
#endif
                intT incounter = 0;
                intT ind = V[i].getInDegree();
                for (intT j = 0; j < ind; j++)
                {
#ifndef WEIGHTED
                    intT ngh = V[i].getInNeighbor(j);
#else
                    intT ngh = V[i].getInNeighbor(j);
                    intT wgh = V[i].getInWeight(j);
#endif
                    if (rangeLow <= ngh && ngh < rangeHi)
                    {
#ifndef WEIGHTED
                        localInEdges[incounter] = ngh;
#else
                        localInEdges[incounter*2] = ngh;
                        localInEdges[incounter*2+1] = wgh;
#endif
                        incounter++;

                    }
                }
                FG.V[i].setInNeighbors(localInEdges);
            }
        }
    }
#if COMPRESSED_CSCCSR
    intT nnz = 0, nnzi = 0, nnzo = 0, ne = 0;
    for( intT i=0; i < n; ++i )
    {
        if( FG.V[i].getInDegree() != 0 || FG.V[i].getOutDegree() != 0 )
        {
            ++nnz;
            if( FG.V[i].getInDegree() != 0 )
                ++nnzi;
            if( FG.V[i].getOutDegree() != 0 )
                ++nnzo;
        }
        ne += FG.V[i].getOutDegree();
    }
    FG.CSCVn=nnzi;
    FG.CSRVn=nnzo;

    FG.CSCV = new pair<intT,vertex> [FG.CSCVn];
    FG.CSRV = new pair<intT,vertex> [FG.CSRVn];


    intT k=0;
    for( intT i=0; i < n; ++i )
    {
        if( FG.V[i].getInDegree() != 0 )
        {
            FG.CSCV[k].first = i;
            FG.CSCV[k].second = FG.V[i];
            k++;
        }
    }
    intT j=0;
    for( intT i=0; i < n; ++i )
    {
        if( FG.V[i].getOutDegree() != 0 )
        {
            FG.CSRV[j].first = i;
            FG.CSRV[j].second = FG.V[i];
            j++;
        }
    }
    cerr << "CSCCSRCompressed graph n=" << n << " Compressed=" << nnz
         << " nnzi=" << nnzi << " nnzo=" << nnzo << " out-edges=" << ne << " compressed" << endl;

#endif // COMPRESSED_CSCCSR
#if COMPRESSED_VERTICES
    intT nnz = 0, nnzi = 0, nnzo = 0, ne = 0;
    for( intT i=0; i < n; ++i )
    {
        if( FG.V[i].getInDegree() != 0 || FG.V[i].getOutDegree() != 0 )
        {
            ++nnz;
            if( FG.V[i].getInDegree() != 0 )
                ++nnzi;
            if( FG.V[i].getOutDegree() != 0 )
                ++nnzo;
        }
        ne += FG.V[i].getOutDegree();
    }
    if( nnzo < n/2 || true )   // HEURISTIC TO USE THESE IN DENSE TRAVERSAL
    {
        FG.CVn = nnz;

        FG.CV = new pair<intT,vertex> [FG.CVn];

        intT k=0;
        for( intT i=0; i < n; ++i )
            if( FG.V[i].getInDegree() != 0 || FG.V[i].getOutDegree() != 0 )
            {
                FG.CV[k].first = i;
                FG.CV[k].second = FG.V[i];
                k++;
            }
        cerr << "Compressed graph n=" << n << " Compressed=" << FG.CVn
             << " nnzi=" << nnzi << " nnzo=" << nnzo << " out-edges=" << ne << " compressed" << endl;
    }
    else
    {
        cerr << "Compressed graph n=" << n << " nonzeroes=" << nnz
             << " nnzi=" << nnzi << " nnzo=" << nnzo << " out-edges=" << ne << " not built" << endl;
    }
#else
    //For push operator, empty vertex check of outdegree
    for( intT i=n-1, next_id=n; i >= 0; --i )
    {
        FG.V[i].setOutDegId(next_id);
        if( FG.V[i].getOutDegree() != 0 )
            next_id = i;
    }

    // Remember next vertex with non-zero indegree
    for( intT i=n-1, next_id=n; i >= 0; --i )
    {
        FG.V[i].setInDegId(next_id);
        if( FG.V[i].getInDegree() != 0 )
            next_id = i;
    }
#endif

    delete [] offsets;
    delete [] inOffsets;

    FG.n_src_start = rangeLow;
    FG.n_dst_start = rangeLow;
    delete [] FG.V;
    return FG;
}

// This routine places an edge in the home graph partition of the destination.
template <class vertex>
graph<vertex>
partitioned_graph<vertex>::graphFilter(wholeGraph<vertex> &GA,
                                       short * partition,
                                       intT p_num)
{
    vertex *V = GA.V;
    const intT n = GA.n;
    bool isSymmetric = GA.isSymmetric;
    intT *counters = new intT [n];
    intT *offsets = new intT [n];
    intT *inCounters = new intT [n];
    intT *inOffsets = new intT [n];
    parallel_for (intT i = 0; i < n; i++)
    {
        intT d = V[i].getOutDegree();
        counters[i] = 0;
        for (intT j = 0; j < d; j++)
        {
            intT ngh = V[i].getOutNeighbor(j);
            if( partition[ngh] == p_num )
                counters[i]++;
        }

        // This is redundant for symmetricVertex.
        if(!isSymmetric)
        {
            d = V[i].getInDegree();
            inCounters[i] = 0;
            for (intT j = 0; j < d; j++)
            {
                intT ngh = V[i].getInNeighbor(j);
                if( partition[ngh] == p_num )
                    inCounters[i]++;
            }
        }
    }

    intT totalSize = 0;
    intT totalInSize = 0;
    for (intT i = 0; i < n; i++)
    {
        offsets[i] = totalSize;
        totalSize += counters[i];

        inOffsets[i] = totalInSize;
        totalInSize += inCounters[i];
    }

    graph<vertex> FG(n, totalSize, isSymmetric, p_num);
    FG.V = new vertex [n];
    delete [] counters;
    delete [] inCounters;

    {
        parallel_for (intT i = 0; i < n; i++)
        {
            FG.V[i].setOutDegree(counters[i]);
            if(!isSymmetric)
                FG.V[i].setInDegree(inCounters[i]);
        }
    }
    intE *edges = FG.allocatedInplace;
    intE *inEdges = FG.inEdges;//0;
    if(!isSymmetric)
    {
        inEdges = new intE [totalInSize];
    }
    // std::cerr << "(out)edges in partition: " << totalSize << std::endl;

    {
        parallel_for (intT i = 0; i < n; i++)
        {
            intE *localEdges = &edges[offsets[i]];
            intT counter = 0;
            intT d = V[i].getOutDegree();
            for (intT j = 0; j < d; j++)
            {
                intT ngh = V[i].getOutNeighbor(j);
                if( partition[ngh] == p_num )
                {
                    localEdges[counter] = ngh;
                    counter++;
                }
            }
            FG.V[i].setOutNeighbors(localEdges);
            if(!isSymmetric)
            {
                intE *localInEdges = &inEdges[inOffsets[i]];
                counter = 0;
                d = V[i].getInDegree();
                for (intT j = 0; j < d; j++)
                {
                    intT ngh = V[i].getInNeighbor(j);
                    if( partition[ngh] == p_num )
                    {
                        localInEdges[counter] = ngh;
                        counter++;
                    }
                }

                FG.V[i].setInNeighbors(localInEdges);
            }
        }
    }

#if COMPRESSED_VERTICES
    intT nnz = 0, nnzi = 0, nnzo = 0, ne = 0;
    for( intT i=0; i < n; ++i )
    {
        if( FG.V[i].getInDegree() != 0 || FG.V[i].getOutDegree() != 0 )
        {
            ++nnz;
            if( FG.V[i].getInDegree() != 0 )
                ++nnzi;
            if( FG.V[i].getOutDegree() != 0 )
                ++nnzo;
        }
        ne += FG.V[i].getOutDegree();
    }

    if( nnzo < n/2 || true )   // HEURISTIC TO USE THESE IN DENSE TRAVERSAL
    {
        FG.CVn = nnz;

        FG.CV = new pair<intT,vertex> [FG.CVn];

        intT k=0;
        for( intT i=0; i < n; ++i )
            if( FG.V[i].getInDegree() != 0 || FG.V[i].getOutDegree() != 0 )
            {
                FG.CV[k].first = i;
                FG.CV[k].second = FG.V[i];
                k++;
            }
        cerr << "Compressed graph n=" << n << " Compressed=" << FG.CVn
             << " nnzi=" << nnzi << " nnzo=" << nnzo << " out-edges=" << ne << " compressed" << endl;
    }
    else
    {
        cerr << "Compressed graph n=" << n << " nonzeroes=" << nnz
             << " nnzi=" << nnzi << " nnzo=" << nnzo << " out-edges=" << ne << " not built" << endl;
    }
#else
    //For push operator, empty vertex check of outdegree
    for( intT i=n-1, next_id=n; i >= 0; --i )
    {
        FG.V[i].setOutDegId(next_id);
        if( FG.V[i].getOutDegree() != 0 )
            next_id = i;
    }

    // Remember next vertex with non-zero indegree
    for( intT i=n-1, next_id=n; i >= 0; --i )
    {
        FG.V[i].setInDegId(next_id);
        if( FG.V[i].getInDegree() != 0 )
            next_id = i;
    }
#endif
    delete [] FG.V;
    return FG;
}

// ======================================================================
// Graph Partitioning
// ======================================================================

//Graph partition for collection Shard size

template <class vertex>
void partitioned_graph<vertex>::partitionByDegree(wholeGraph<vertex> GA, int numOfNode, intT *sizeArr, int sizeOfOneEle, bool useOutDegree/*=false*/)
{
    const intT n = GA.n;
    intT *degrees = new intT [n];

    //int arraySize = n / numOfNode;

    if (useOutDegree)
    {
        {
            parallel_for(intT i = 0; i < n; i++) degrees[i] = GA.V[i].getOutDegree();
        }
    }
    else
    {
        {
            parallel_for(intT i = 0; i < n; i++) degrees[i] = GA.V[i].getInDegree();
        }
    }

    intT edges[numOfNode];
    for (int i = 0; i < numOfNode; i++)
    {
        edges[i] = 0;
        sizeArr[i] = 0;
    }

    intT totalDegree = 0;
    for (intT i = 0; i < n; i++)
    {
        totalDegree += degrees[i];
    }

    intT averageDegree = totalDegree / numOfNode;
    //printf("average is %d\n", averageDegree);
    intT counter = 0;
    intT tmpSizeCounter = 0;
    for (intT i = 0; i < n; i++)
    {
        edges[counter]+=degrees[i];
        sizeArr[counter]++;
        tmpSizeCounter++;
        if (edges[counter]>=averageDegree && counter <numOfNode-1)
        {
            counter++;
            tmpSizeCounter=0;
        }
    }
    //assert( counter+1 == numOfNode );

    delete [] degrees;

}

template <class vertex>
void partitioned_graph<vertex>::partitionApprox(wholeGraph<vertex> & GA, short * partitions, int np)
{
    const intT n = GA.n;
    intT *load = new intT [np];
    intT *cost = new intT [np];

    for( int i=0; i < np; ++i ) load[i] = 0;
    parallel_for( int i=0; i < n; ++i ) partitions[i] = -1;

    for( int i=0; i < n; ++i )
    {
        vertex & V = GA.V[i];
        if( partitions[i] < 0 )
        {
            if( 0 && V.getOutDegree() > 2 * np )   // highly connected - spread
            {
                cerr << " high " << V.getOutDegree() << endl;
                intT lp = np - 1;
                for( int j=0; j < np; ++j )
                    if( load[j] < load[lp] )
                        lp = j;
                if( partitions[i] < 0 )
                {
                    partitions[i] = lp;
                    load[lp]++;
                }
                for( int j=0; j < V.getOutDegree(); ++j )
                {
                    int jj = V.getOutNeighbor(j);
                    if( partitions[jj] < 0 )
                    {
                        intT lp = 0;
                        for( int k=1; k < np; ++k )
                            if( load[k] < load[lp] )
                                lp = k;
                        partitions[jj] = lp;
                        load[lp]++;
                    }
                }
            }
            else     // few connections, dense
            {
                int p = -1;
                for( int j=0; j < np; ++j ) cost[j] = 0;
                for( int j=0; j < V.getOutDegree(); ++j )
                {
                    int jj = V.getOutNeighbor(j);
                    if( partitions[jj] >= 0 )
                        cost[partitions[jj]]++;
                }
                // Find partition with highest cost
                p = 0;
                for( int j=1; j < np; ++j )
                    // if( cost[j] > cost[p] )
                    if( cost[j] > 0 && cost[j] < cost[p] ) // least non-zero cost
                        p = j;
                // Find partition with least load
                intT lp = 0;
                for( int j=1; j < np; ++j )
                    if( load[j] < load[lp] )
                        lp = j;
                if( cost[p] == 0 )
                    p = lp;

                // Place on p
                assert( p >= 0 && p < np );
                if( partitions[i] < 0 )
                {
                    partitions[i] = lp;
                    load[lp]++;
                }
                for( int j=0; j < V.getOutDegree(); ++j )
                {
                    int jj = V.getOutNeighbor(j);
                    if( partitions[jj] < 0 )
                    {
                        partitions[jj] = p;
                        load[p]++;
                    }
                }
            }
        }
    }

    int nz = 0, no = 0, na = 0;
    for( int i=0; i < n; ++i )
    {
        vertex & V = GA.V[i];
        assert( partitions[i] >= 0 );
        for( int j=0; j < np; ++j ) cost[j] = 0;
        for( int j=0; j < V.getOutDegree(); ++j )
        {
            int jj = V.getOutNeighbor(j);
            assert( partitions[jj] >= 0 );
            cost[partitions[jj]]++;
        }
        for( int j=0; j < np; ++j )
        {
            if( cost[j] == 0 ) ++nz;
            else if( cost[j] == 1 ) ++no;
            else ++na;
        }
    }
    //cerr << "partitioning: m=" << GA.m << " no=" << no << " nz=" << nz << " others=" << na << endl;

    // Rewriting graph
    {
        intT s = 0;
        for( int j=0; j < np; ++j )
            s += load[j];
        assert( s == n );
    }
    for( int j=0; j < np; ++j )
        cost[j] = load[j];

    sequence::plusScan(cost, cost, np);

    intT * shuffle = new intT [n];
    parallel_for( int i=0; i < n; ++i )
    {
        intT w = partitions[i];
        shuffle[i] = __sync_fetch_and_add(&cost[w],1);
        assert( 0 <= shuffle[i] && shuffle[i] < n );
    }

    GA.reorder_vertices( shuffle );

    delete [] shuffle;
    delete [] load;
    delete [] cost;
}

#if 0
// Not applicable - partitioning edges instead of vertices
template <class vertex>
void partitioned_graph<vertex>::partitionGraphLab(wholeGraph<vertex> & GA, short * partitions, int np)
{
    intT n = GA.n;
    intT * assigned = new intT, n);
    size_t * present = new size_t, n);
    intT * load = new intT,np);

    parallel_for( int i=0; i < n; ++i )
    {
        assigned[i] = 0;
        present[i] = 0;
    }
    for( int i=0; i < np; ++i )
    {
        load[i] = 0;
    }

    for( int i=0; i < n; ++i )
    {
        int d = GA.V[i].getOutDegree();
        for( int j=0; j < d; ++j )
        {
            int u = i;
            int v = GA.V[i].getOutNeighbor(j);

            size_t Au = present[u];
            size_t Av = present[v];

            size_t Auv = Au & Av;
            if( Auv )
            {
                // select from Auv;
            }
            else if( Au && Av )
            {
                if( GA.V[u].getOutDegree() - assigned[u]
                        > GA.V[v].getOutDegree() - assigned[v] )
                    Auv = Au;
                else
                    Auv = Av;
            }
            else if( !Au || !Av )
            {
                Auv = Au | Av;
            }
            else
                Auv = ( size_t(1) << np ) - 1; // any

            int wc = -1;
            size_t Ap = 1;
            for( int w=0; w < np; ++w, Ap <<= 1 )
                if( Auv & Ap )
                {
                    if( wc >= 0 && load[wc] > load[w] )
                        wc = w;
                }

            assert( wc >= 0 );
            Ap = size_t(1) << wc;
            present[u] |= Ap;
            present[v] |= Ap;
            assigned[u]++;
            assigned[v]++;
            load[wc]++;
        }
    }

    free(assigned);
    free(present);
    free(load);
}
#endif

template<>
void graph<symmetricVertex>::reorder_vertices( intT * __restrict reorder )
{
    symmetricVertex * newV = new symmetricVertex [n];

#ifndef WEIGHTED
    intE * newE = new intE [m];
#else
    intE * newE = new intE [2*m];
#endif

    intT * start = new intT [n];
    parallel_for( intT invi=0; invi < n; ++invi )
    {
        intT i = reorder[invi];
        intT d = V[invi].getOutDegree();
        start[i] = d;

        newV[invi].setOutDegree( 0 );
        newV[invi].setOutNeighbors( 0 );
    }
    sequence::plusScan(start, start, n);

    /*parallel_*/for( intT invi=0; invi < n; ++invi )
    {
        intT i = reorder[invi];
        intT d = V[invi].getOutDegree();
        intE * ne = &newE[start[i]];
        if( d > 1000 )
        {
            parallel_for( intT j=0; j < d; ++j )
            ne[j] = reorder[V[invi].getOutNeighbor(j)];
            mysort( &ne[0], &ne[d], std::less<intT>() );
        }
        else
        {
            for( intT j=0; j < d; ++j )
                ne[j] = reorder[V[invi].getOutNeighbor(j)];
            mysort( &ne[0], &ne[d], std::less<intT>() );
        }
        newV[i].setOutDegree( d );
        newV[i].setOutNeighbors( ne );
    }

    delete [] start;
    delete [] V;
    delete [] allocatedInplace;
    V = newV;
    allocatedInplace = newE;
}

template<>
void graph<asymmetricVertex>::reorder_vertices( intT * __restrict reorder )
{
    asymmetricVertex * newV = new asymmetricVertex [n];

#ifndef WEIGHTED
    intE * newInE = new intE [m];
    intE * newOutE = new intE [m];
#else
    intE * newInE = new intE [2*m];
    intE * newOutE = new intE [2*m];
#endif

    intT * startOut = new intT [n];
    intT * startIn = new intT [n];
    parallel_for( intT invi=0; invi < n; ++invi )
    {
        intT i = reorder[invi];
        intT d = V[invi].getOutDegree();
        startOut[i] = d;
        d = V[invi].getInDegree();
        startIn[i] = d;
    }
    sequence::plusScan(startOut, startOut, n);
    sequence::plusScan(startIn, startIn, n);


    // Do not shuffle edges around - will be written to disk
    parallel_for( intT invi=0; invi < n; ++invi )
    {
        intT i = reorder[invi];
        intE * oe = &newOutE[startOut[i]];
        intT od = V[invi].getOutDegree();
        if( od > 1000 )
        {
            parallel_for( intT j=0; j < od; ++j )
            mysort( &oe[0], &oe[od], std::less<intT>() );
        }
        else
        {
            for( intT j=0; j < od; ++j )
                oe[j] = reorder[V[invi].getOutNeighbor(j)];
            mysort( &oe[0], &oe[od], std::less<intT>() );
        }
        newV[i].setOutDegree( od );
        newV[i].setOutNeighbors( oe );

        intE * ie = &newInE[startIn[i]];
        intT id = V[invi].getInDegree();
        if( id > 1000 )
        {
            parallel_for( intT j=0; j < id; ++j )
            ie[j] = reorder[V[invi].getInNeighbor(j)];
            mysort( &ie[0], &ie[id], std::less<intT>() );
        }
        else
        {
            for( intT j=0; j < id; ++j )
                ie[j] = reorder[V[invi].getInNeighbor(j)];
            mysort( &ie[0], &ie[id], std::less<intT>() );
        }
        newV[i].setInDegree( id );
        newV[i].setInNeighbors( ie );
    }

    delete [] startIn;
    delete [] startOut;
    delete [] V;
    delete [] allocatedInplace;
    delete [] inEdges;
    V = newV;
    allocatedInplace = newOutE;
    inEdges = newInE;
}

