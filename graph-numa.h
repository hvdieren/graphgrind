// -*- C++ -*-
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include "parallel.h"
#include <assert.h>
#include "numa_page_check.h"
#include "mm.h"
//#include <type_traits>
#include <unistd.h>
#include <sched.h>
#include <errno.h>
#include <cstring>
#include <string>
#include <utility>
#include <algorithm>

#include <sys/mman.h>
#include <numaif.h>
#include <numa.h>

#ifndef VERTEX
#define VERTEX 0
#endif

#ifndef VERTEX_BASED
#define VERTEX_BASED 0
#endif

#ifndef COMPRESSED_VERTICES
#define COMPRESSED_VERTICES 0
#endif

#ifndef PULL_EDGES
#define PULL_EDGES 0
#endif
#ifndef EDGES_HILBERT
#define EDGES_HILBERT 0
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

// wholeGraph for whole graph loading
// and sparse iteration graph traversal
// uses NUMA interleave to allocate
template <class vertex>
class wholeGraph
{
public:
    mmap_ptr<vertex> V;  //mmap
    intT n;
    intT m;
    mmap_ptr<intT> flags;
    mmap_ptr<intE> allocatedInplace;
    mmap_ptr<intE> inEdges;
    bool transposed;
    bool isSymmetric;

    wholeGraph() {}
    wholeGraph(intT nn, intT mm, bool issym)
        : n(nn), m(mm), isSymmetric(issym),
          transposed(false)
    {

//NUMA_AWARE and Ligra_normal without partition
           V.Interleave_allocate(n);
#ifndef WEIGHTED
	   allocatedInplace.Interleave_allocate(m);
#else//WEIGHTED
           allocatedInplace.Interleave_allocate(2*m);
#endif
            if(!isSymmetric)
            {
#ifndef WEIGHTED
                inEdges.Interleave_allocate(m);
#else//WEIGHTED
                inEdges.Interleave_allocate(2*m);
#endif
            }
    }


    void del()
    {
        flags.del();
        allocatedInplace.del();
        V.del();
        inEdges.del();
    }

    void reorder_vertices( intT * __restrict reorder );
//    asymmetricVertex * newV = new asymmetricVertex [n];

    void transpose()
    {
        if(!isSymmetric)
        {
            parallel_for(intT i=0; i<n; i++)
		V[i].flipEdges();
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

// graph for partitioned_graph to do graph partitioning
// each huge array use NUMA-Aware node allocation
// to corresponding node localG[i] to Node i
// Dense iteration use this
template <class vertex>
class graph
{
public:
    mmap_ptr<vertex> V;  //mmap
    intT n;
    intT m;
    mmap_ptr<intT> flags;
    mmap_ptr<intE> allocatedInplace;
    mmap_ptr<intE> inEdges;
    bool transposed;
    intT n_src_start, n_src_end, n_dst_start, n_dst_end;
    bool isSymmetric;
    int numanode;

//Compressed partitioned_graph after partitioning,
//remove the empty vertices,
//restruct the vertext set,
//CSC for backwark and CSR for forward;
#if COMPRESSED_CSCCSR
    mmap_ptr< pair<intT,vertex> > CSCV;
    mmap_ptr< pair<intT,vertex> > CSRV;
    intT CSCVn;
    intT CSRVn;
#endif
#if COMPRESSED_VERTICES
    mmap_ptr< <intT,vertex> > CV;    //mmap
    intT CVn;
#endif
    graph() {}
    graph(intT nn, intT mm, bool issym,int pp)
        : n(nn), m(mm), isSymmetric(issym),
          transposed(false),
          n_src_start(0), n_src_end(nn),
          n_dst_start(0), n_dst_end(nn)
          ,numanode(pp)
    {

//NUMA_AWARE and Ligra_normal without partition
#ifndef WEIGHTED
	   allocatedInplace.local_allocate(m,numanode);
#else//WEIGHTED
           allocatedInplace.local_allocate(2*m,numanode);
#endif
            if(!isSymmetric)
            {
#ifndef WEIGHTED
                inEdges.local_allocate(m,numanode);
#else//WEIGHTED
                inEdges.local_allocate(2*m,numanode);
#endif
            }
#if COMPRESSED_CSCCSR
        CSCVn=0;
        CSRVn=0;
#endif
#if COMPRESSED_VERTICES
        CVn = 0;
#endif
    }

    void del()
    {
        flags.del();
        allocatedInplace.del();
        inEdges.del();
        //V.del();
#if COMPRESSED_CSCCSR
        CSCV.del();
        CSRV.del();
#endif
#if COMPRESSED_VERTICES
        CV.del();
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

//Graph partitioning, contain partitioned graph,
//partitioner value
//Select there partitioning method,
//Mix partition (only four partition), first partitioning 
//as two partition by souce and then repartitioning by destination
//from two to four 
//PartitionBySour : Partition By source, all the 
//home vertices contain out-degrees, improve locality 
//for backward, BFS, BC, Component and PageRank
//PartitionByDest : Partitino By destination, all the
//home vertices contain in-degrees, improve locality
//for forward, PageRankDelta, SPMV, BP, BellmanFord
template <class vertex>
class partitioned_graph
{
public:
    typedef vertex vertex_type;

    partitioner partition;
    partitioner vertex_partition;
    graph<vertex> * localG;
    intT m,n;
    bool source;
private:
    // All variables should be private
#if PULL_EDGES
    EdgeList<Edge> * localEdgeList;
#endif
    wholeGraph<vertex> & GA;

public:
    partitioned_graph( wholeGraph<vertex> & GA_, int num_part,
                       bool partition_source )
        : partition(num_part,GA_.n),m(GA_.m),vertex_partition(num_part,GA_.n),n(GA_.n),source(partition_source),
	  GA(GA_)
    {

        localG = new graph<vertex> [num_part];
        const int perNode = partition.get_num_per_node_partitions();
#if PULL_EDGES
        localEdgeList = new EdgeList<Edge>[num_part];
#endif
        //std::cerr << "partition graph: n=" << GA.n << " m=" << GA.m << std::endl;

#if PARTITION_RANGE
        // Normal path
#if MIX
        {
            int ini=2;
            int mix = num_part/ini;
	    assert( mix * ini == num_part );
            partitioner by2(ini,GA.n);
            partitionByDegree( GA, ini, by2.as_array(), true );
            by2.compute_starts();

            wholeGraph<vertex>* tempGA = new wholeGraph<vertex> [ini];
            parallel_for ( int i =0 ; i < ini; i++)
                tempGA[i] = PartitionBySourW( GA, by2.start_of(i), by2.start_of(i+1), 0 );

            partitioner* per_src = new partitioner [ini];
	    for( int src=0; src < ini; ++src ) {
                per_src[src] = partitioner(mix,GA.n);
		partitionByDegree( tempGA[src], mix, per_src[src].as_array(), false );
		per_src[src].compute_starts();
            }
            
                parallel_for_numa( int i=0; i < num_numa_node; ++i )
                {
                  parallel_for( int p = perNode*i; p < perNode*(i+1); ++p )
                  {
                    int src=p % ini; // p / mix; // iterates from 0..ini; likely within NUMA node
                    int k=p / ini; // p % mix; // iterates from 0..mix; across NUMA node
                    localG[p] = PartitionByDest( tempGA[src], per_src[src].start_of(k), per_src[src].start_of(k+1), i );
                  }
                }
            
		for(int i=0; i<mix; i++) {
		    int b = 0;
		    for( int j=0; j < ini; ++j )
			b += per_src[j].as_array()[i];
		    b = b/ini;
		    int rem = b;
		    for( int j=0; j < ini-1; ++j ) {
			partition.as_array()[i*ini+j] = b/ini;
			rem -= b/ini;
		    }
		    partition.as_array()[i*ini+ini-1] = rem;
		}
		partition.as_array()[num_part] = GA.n;
                
                partition.compute_starts();
            
                for (int i = 0; i < ini; i++)
                   tempGA[i].del();
                delete [] per_src;
        }
#else //MIX
#if VERTEX
        partitionByVertex( GA, num_part, partition.as_array(), partition_source );
        partition.compute_starts();
#else
        partitionByDegree( GA, num_part, partition.as_array(), partition_source );
        partition.compute_starts();
#endif
        partitionByVertex( GA, num_part, vertex_partition.as_array(), partition_source );
        vertex_partition.compute_starts();
#if CPU_PARTITION
        partition.compute_vertexstarts();
#endif
//#pragma cilk numa(strict)
            parallel_for_numa( int i=0; i < num_numa_node; ++i )   //same loop with allocation
            {
                parallel_for( int p = perNode*i; p < perNode*(i+1); ++p )
                {
                    if( partition_source )
                        localG[p] = PartitionBySour( GA, partition.start_of(p), partition.start_of(p+1),i);
                    else
                        localG[p] = PartitionByDest( GA, partition.start_of(p), partition.start_of(p+1),i);
#if PULL_EDGES
                    localEdgeList[p] = localG[p].build_edgelist();
#if EDGES_HILBERT
                    localEdgeList[p].hilbert_sort();
#endif
#endif
                    //numa_page_check<vertex,intT>(localG[p].V,n,p);
                }
            }
#endif
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
#if VERTEX_BASED
    const partitioner & get_partitioner() const
    {
        return vertex_partition;
    }
#else
    const partitioner & get_partitioner() const
    {
        return partition;
    }
#endif
    wholeGraph<vertex> & getWholeGraph() { return GA; }

    void transpose()
    {
	// Do transpositions on NUMA node storing graph partition
       // parallel_for_numa(unsigned i=0; i<partition.get_num_partitions(); i++)
        const int perNode = partition.get_num_per_node_partitions();
        parallel_for_numa( int i=0; i < num_numa_node; ++i )
        {
          parallel_for( int p = perNode*i; p < perNode*(i+1); ++p )
            localG[p].transpose();
        }
	GA.transpose();
    }
    bool transposed() const {
	return GA.transposed; // Assume status of partitions matches whole graph
    }

private:
// These functions are specific to the partitioned graph. No one needs
// to know about them, I think.
// graphFilter is polymer method, partitioend graph:
// Backward: each local vertex contains its own out-degree
// Forward: each local vertex contains its own in-degree
// PartitionByDest :
// Backward: each local vertex contains its own in-degree
// Forward:  each local vertex contains its own in-degree
// PartitionBySour:
// Backward: each local vertex contains its own out-degree
// Forward: each local vertex contains its own out-degree
    void partitionByDegree( wholeGraph<vertex> GA, int numOfNode, intT *sizeArr,
                            bool useOutDegree );
    void partitionByVertex( wholeGraph<vertex> GA, int numOfNode, intT *sizeArr,
                            bool useOutDegree );
    void partitionApprox(wholeGraph<vertex> & GA, short * partitions, int np);
    graph<vertex> graphFilter( wholeGraph<vertex>& GA, int rangeLow, int rangeHi ,int numanode);
    graph<vertex> graphFilter( wholeGraph<vertex>& GA, short * partition, intT p );
    graph<vertex> PartitionByDest(wholeGraph<vertex>& GA, int rangeLow, int rangeHi ,int numanode);
    graph<vertex> PartitionBySour(wholeGraph<vertex>& GA, int rangeLow, int rangeHi ,int numanode);
    wholeGraph<vertex> PartitionBySourW(wholeGraph<vertex>& GA, int rangeLow, int rangeHi ,int numanode);
    wholeGraph<vertex> PartitionByDestW(wholeGraph<vertex>& GA, int rangeLow, int rangeHi ,int numanode);
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
                counters[i] = V[i].getOutDegree();
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
    FG.V.local_allocate(n,numanode);

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
    intT nnz = 0, nnzi = 0, nnzo = 0, ne = 0, ni = 0;
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
    }
    FG.CSCVn=nnzi;
    FG.CSRVn=nnzo;

    FG.CSCV.local_allocate(FG.CSCVn,numanode);
    FG.CSRV.local_allocate(FG.CSRVn,numanode);

    intT k=0;
    for( intT i=0; i < n; ++i )
    {
        if( FG.V[i].getInDegree() != 0 )
        {
            FG.CSCV[k].first = i;
            FG.CSCV[k].second = FG.V[i];
            ni += FG.V[i].getInDegree();
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
            ne += FG.V[i].getOutDegree();
            j++;
        }
    }
#if 0
    cerr << "CSCCSRCompressed graph n=" << n << " Compressed=" << nnz
         << " CSC=" << nnzi << " CSR=" << nnzo << " out-edges=" << ne 
         << " in-edges" << ni << " compressed" << endl;
#endif // COMPRESSED_CSCCSR
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

        FG.CV.local_allocate(FG.CVn,numanode);


        intT k=0;
        for( intT i=0; i < n; ++i )
            if( FG.V[i].getInDegree() != 0 || FG.V[i].getOutDegree() != 0 )
            {
                FG.CV[k].first = i;
                FG.CV[k].second = FG.V[i];
                k++;
            }
        //cerr << "Compressed graph n=" << n << " Compressed=" << FG.CVn
         //    << " nnzi=" << nnzi << " nnzo=" << nnzo << " out-edges=" << ne << " compressed" << endl;
    }
    else
    {
        //cerr << "Compressed graph n=" << n << " nonzeroes=" << nnz
          //   << " nnzi=" << nnzi << " nnzo=" << nnzo << " out-edges=" << ne << " not built" << endl;
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
    FG.V.del();
    return FG;//return graph<vertex>(newVertexSet, GA.n, GA.m);
}
//This is used for MIX method 
//get the wholeGraph for next source partition
template <class vertex>
wholeGraph<vertex> partitioned_graph<vertex>::PartitionByDestW(wholeGraph<vertex> &GA, int rangeLow, int rangeHi,int numanode)
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

    if(!isSymmetric)
        assert( totalSize == totalInSize );

    wholeGraph<vertex> FG(n, totalSize, isSymmetric);
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

    delete [] offsets;
    delete [] inOffsets;
    return FG;
}

template <class vertex>
wholeGraph<vertex> partitioned_graph<vertex>::PartitionBySourW(wholeGraph<vertex> &GA, int rangeLow, int rangeHi,int numanode)
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
            counters[i] = 0;
            inCounters[i] = 0;
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
              if(rangeLow <= i && i < rangeHi)
                  counters[i] = V[i].getOutDegree();
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

    wholeGraph<vertex> FG(n, totalSize, isSymmetric);
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

    delete [] offsets;
    delete [] inOffsets;
    // We know that all vertices with incoming edges are located within this
    // partition, so when iterating over destinations of edges, we can limit the
    // range to:
    //FG.n_dst_start = rangeLow;
    //FG.n_dst_end = rangeHi;
    return FG;
}
//This method is designed for indegree, and out model use now
//all the vertices of graph will store part outdegree (CSR) for push
//local vertices of partitioend graph will store all in-degree (CSC) for pull
// Place edge (u,v) in partition of v
template <class vertex>
graph<vertex> partitioned_graph<vertex>::PartitionByDest(wholeGraph<vertex> &GA, int rangeLow, int rangeHi,int numanode)
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

    if(!isSymmetric)
        assert( totalSize == totalInSize );

    graph<vertex> FG(n, totalSize, isSymmetric, numanode);
    FG.V.local_allocate(n,numanode);
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
    intT nnz = 0, nnzi = 0, nnzo = 0, ne = 0, ni=0;
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
    }
    FG.CSCVn=nnzi;
    FG.CSRVn=nnzo;

    FG.CSCV.local_allocate(FG.CSCVn,numanode);
    FG.CSRV.local_allocate(FG.CSRVn,numanode);

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
            ne += FG.V[i].getOutDegree();
            j++;
        }
    }
#if 0
    cerr << "CSCCSRCompressed graph n=" << n << " Compressed=" << nnz
         << " CSC=" << nnzi << " CSR=" << nnzo << " out-edges=" << ne 
         << " compressed" << endl;
#endif // COMPRESSED_CSCCSR
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
        FG.CV.local_allocate(FG.CVn,numanode);

        intT k=0;
        for( intT i=0; i < n; ++i )
            if( FG.V[i].getInDegree() != 0 || FG.V[i].getOutDegree() != 0 )
            {
                FG.CV[k].first = i;
                FG.CV[k].second = FG.V[i];
                k++;
            }
        //cerr << "Compressed graph n=" << n << " Compressed=" << FG.CVn
          //   << " nnzi=" << nnzi << " nnzo=" << nnzo << " out-edges=" << ne << " compressed" << endl;
    }
    else
    {
        //cerr << "Compressed graph n=" << n << " nonzeroes=" << nnz
          //   << " nnzi=" << nnzi << " nnzo=" << nnzo << " out-edges=" << ne << " not built" << endl;
    }
#endif
#if !COMPRESSED_VERTICES && !COMPRESSED_CSCCSR
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
#if COMPRESSED_VERTICES || COMPRESSED_CSCCSR
    FG.V.del();
#endif
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
    FG.V.local_allocate(n,numanode);


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

    FG.CSCV.local_allocate(FG.CSCVn,numanode);
    FG.CSRV.local_allocate(FG.CSRVn,numanode);

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

        FG.CV.local_allocate(FG.CVn,numanode);

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
    FG.V.del();
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
    FG.V.local_allocate(n,p_num);
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

        FG.CV.local_allocate(FG.CVn,numanode);

        intT k=0;
        for( intT i=0; i < n; ++i )
            if( FG.V[i].getInDegree() != 0 || FG.V[i].getOutDegree() != 0 )
            {
                FG.CV[k].first = i;
                FG.CV[k].second = FG.V[i];
                k++;
            }
        //cerr << "Compressed graph n=" << n << " Compressed=" << FG.CVn
          //   << " nnzi=" << nnzi << " nnzo=" << nnzo << " out-edges=" << ne << " compressed" << endl;
    }
    else
    {
       // cerr << "Compressed graph n=" << n << " nonzeroes=" << nnz
         //    << " nnzi=" << nnzi << " nnzo=" << nnzo << " out-edges=" << ne << " not built" << endl;
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
    FG.V.del();
    return FG;
}

// ======================================================================
// Graph Partitioning
// ======================================================================
//vertex balancing partition

template <class vertex>
void partitioned_graph<vertex>::partitionByVertex(wholeGraph<vertex> GA, int numOfNode, intT *sizeArr,  bool useOutDegree)
{
    intT chunck=0;
    for ( intT i=0; i < numOfNode-1; ++i)
    {
        sizeArr[i]=GA.n/numOfNode;
        chunck+=sizeArr[i];
    }
    sizeArr[numOfNode-1]=GA.n-chunck;
}

template <class vertex>
void partitioned_graph<vertex>::partitionByDegree(wholeGraph<vertex> GA, int numOfNode, intT *sizeArr,  bool useOutDegree)
{
    const intT n = GA.n;
    intT *degrees = new intT [n];
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
    
    intT averageDegree = GA.m / numOfNode;
    int counter = 0;
    for (intT i = 0; i < n; i++)
    {
        edges[counter]+=degrees[i];
        sizeArr[counter]++;
        if (edges[counter]>=averageDegree && counter <numOfNode-1)
        {
            counter++;
        }
    }
    assert( counter+1 == numOfNode );
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
    //symmetricVertex * newV = new symmetricVertex [n];
      mmap_ptr<symmetricVertex> newV(n,0);
#ifndef WEIGHTED
    //intE * newE = new intE [m];
      mmap_ptr<intE> newE(m,0);
#else
    //intE * newE = new intE [2*m];
     mmap_ptr<intE> newE(2*m,0);
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
    allocatedInplace.del();
    V = newV;
    allocatedInplace = newE;
}

template<>
void graph<asymmetricVertex>::reorder_vertices( intT * __restrict reorder )
{
    //aymmetricVertex * newV = new asymmetricVertex [n];
     mmap_ptr<asymmetricVertex> newV (n,0);
#ifndef WEIGHTED
      mmap_ptr<intE> newInE(m,0);
      mmap_ptr<intE> newOutE(m,0);
    //intE * newInE = new intE [m];
    //intE * newOutE = new intE [m];
#else
    mmap_ptr<intE> newInE(2*m,0);
    mmap_ptr<intE> newOutE(2*m,0);
    //intE * newInE = new intE [2*m];
    //intE * newOutE = new intE [2*m];
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
    allocatedInplace.del();//delete [] allocatedInplace;
    inEdges.del(); //delete [] inEdges;
    V = newV;
    allocatedInplace = newOutE;
    inEdges = newInE;
}

