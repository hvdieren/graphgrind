// -*- C++ -*-
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include "parallel.h"
#include <assert.h>
#include <unistd.h>
#include <sched.h>
#include <errno.h>
#include <cstring>
#include <utility>
#include <algorithm>

#include <numa.h>
#include <numaif.h>
static int num_numa_node=numa_num_configured_nodes();
#ifndef CPU_PARTITION
#define CPU_PARTITION 0
#endif

#if CPU_PARTITION
#define loop(vname,part,PerNode,code) {                                 \
      parallel_for_numa( int i=0; i < num_numa_node; ++i ) {            \
	_Pragma( STRINGIFY(cilk grainsize = _SCAN_BSIZE) ) parallel_for(intT vname=part.vertexstart_of(PerNode*i); vname < part.vertexstart_of(PerNode*(i+1)); ++vname) {  \
		code;                                                   \
	}                                                               \
      }                   						\
}
#else
#define loop(vname,part,PerNode,code) {                                 \
      parallel_for_numa( int i=0; i < num_numa_node; ++i ) {            \
	_Pragma( STRINGIFY(cilk grainsize = _SCAN_BSIZE) ) parallel_for(intT vname=part.start_of(PerNode*i); vname < part.start_of(PerNode*(i+1)); ++vname) {  \
		code;                                                   \
	}                                                               \
      }                   						\
}
#endif

#ifndef PARTITION_RANGE
#define PARTITION_RANGE 1
#endif
//This is the function for hugepage mmap allocation
//Based on the NUMA Awareness
#if PARTITION_RANGE
class partitioner
{
    intT num_partitions;
    intT * partition;
#if CPU_PARTITION
    intT * vstarts;
#endif
    intT * starts;
    int num_per_node;
public:
    // Deep copy semantics: every copy gets a new array
    partitioner() : num_partitions( 0 ), partition( 0 ), starts( 0 ), 
#if CPU_PARTITION
vstarts ( 0 ), 
#endif
num_per_node(0) { }
    partitioner( intT n, intT e ) : num_partitions( n )
    {
        partition = new intT [num_partitions+1];
        starts = new intT [num_partitions+1];
#if CPU_PARTITION
        vstarts = new intT [num_partitions+1];
#endif
        partition[num_partitions] = e;
        num_per_node = num_partitions/num_numa_node;
    }
    partitioner( const partitioner & p ) : num_partitions( p.num_partitions )
    {
        partition = new intT [num_partitions+1];
        starts = new intT [num_partitions+1];
        std::copy( &p.partition[0], &p.partition[num_partitions+1], partition );
        std::copy( &p.starts[0], &p.starts[num_partitions+1], starts );
#if CPU_PARTITION
        vstarts = new intT [num_partitions+1];
        std::copy( &p.vstarts[0], &p.vstarts[num_partitions+1], vstarts );
#endif
        num_per_node = num_partitions/num_numa_node;
    }
    const partitioner & operator = ( const partitioner & p )
    {
        if( partition )
            delete [] partition;
        if( starts )
            delete [] starts;
#if CPU_PARTITION
        if( vstarts )
            delete [] vstarts;
#endif
        num_partitions = p.num_partitions;
        num_per_node = p.num_partitions/num_numa_node;
        partition = new intT [num_partitions+1];
        starts = new intT [num_partitions+1];
        std::copy( &p.partition[0], &p.partition[num_partitions+1], partition );
        std::copy( &p.starts[0], &p.starts[num_partitions+1], starts );
#if CPU_PARTITION
        vstarts = new intT [num_partitions+1];
        std::copy( &p.vstarts[0], &p.vstarts[num_partitions+1], vstarts );
#endif
        return *this;
    }
    ~partitioner()
    {

        if( partition )
            delete [] partition;
        if( starts )
            delete [] starts;
#if CPU_PARTITION
        if( vstarts )
            delete [] vstarts;
#endif
    }
    // For easy interfacing with partitionByDegree()
    intT * as_array()
    {
        return partition;
    }
    int get_num_per_node_partitions() const
    {
        return num_per_node;
    }

    int get_num_partitions() const
    {
        return num_partitions;
    }
    intT get_num_elements() const
    {
        return partition[num_partitions];
    }
    intT set_num_elements(intT i)
    {
        return partition[num_partitions]=i;
    }

    // Translate vertex id to partition
    intT partition_of( intT vertex_id ) const
    {
        intT n = 0;
        for( intT p=0; p < num_partitions; ++p )
        {
            n += partition[p];
            if( vertex_id < n )
                return p;
        }
        abort(); // should not occur unless vertex_id is out of range
    }

    //Get the size of each partition
    intT get_size(intT i) const
    {
        return partition[i];
    }

    //get the start number of each partition
    void compute_starts()
    {
        intT startID=0;
        for( intT i=0; i <= num_partitions; i++ )
        {
            starts[i] = startID;
            startID += partition[i];
        }
    }
    intT start_of(intT i) const
    {
        return starts[i];
    }
#if CPU_PARTITION
    intT vertexstart_of(intT i) const
    {
        return vstarts[i];
    }
    void compute_vertexstarts()
    {
        intT startID=0;
        for( intT i=0; i < num_partitions; i++ )
        {
            vstarts[i] = startID;
            startID += partition[num_partitions]/num_partitions;
        }
        vstarts[num_partitions] = partition[num_partitions];
    }
#endif
    // Get offset of vertex id within its partition
    intT offset_of( intT vertex_id ) const
    {
        intT n = 0;
        for( intT p=0; p < num_partitions; ++p )
        {
            n += partition[p];
            if( vertex_id < n )
                return vertex_id - (n-partition[p]);
        }
        abort(); // should not occur unless vertex_id is out of range
    }

    /* Fancy C++ style iterator
        typedef intT * iterator;
        iterator begin() const { return &partition[0]; }
        iterator end() const { return &partition;}
    */
};
#else
class partitioner
{
    intT num_partitions, elements;
    short * partition;
    intT * size;

public:
    // Deep copy semantics: every copy gets a new array
    partitioner() : num_partitions( 0 ), partition( 0 ), size( 0 ),
        elements( 0 ) { }
    partitioner( intT n, intT e ) : num_partitions( n ), elements( e )
    {
        partition = new short [e];
        size = new intT [num_partitions+1];
    }
    partitioner( const partitioner & p ) : num_partitions( p.num_partitions ),
        elements( p.elements )
    {
        partition = new short [elements];
        size = new intT [num_partitions+1];
        std::copy( &p.partition[0], &p.partition[elements], partition );
        std::copy( &p.size[0], &p.size[num_partitions+1], size );
    }
    const partitioner & operator = ( const partitioner & p )
    {
        if( partition )
            delete [] partition;
        if( starts )
            delete [] starts;
        num_partitions = p.num_partitions;
        elements = p.elements;
        partition = new short [elements];
        size = new intT [num_partitions+1];
        std::copy( &p.partition[0], &p.partition[elements], partition );
        std::copy( &p.size[0], &p.size[num_partitions+1], size );
        return *this;
    }
    ~partitioner()
    {

        if( partition )
            delete [] partition;
        if( starts )
            delete [] starts;
    }
    short *as_array()
    {
        return partition;
    }
    int get_num_partitions() const
    {
        return num_partitions;
    }
    intT get_num_elements() const
    {
        return elements;
    }
    intT set_num_elements(intT i)
    {
        return elements = i;
    }
    void set_size(intT p, intT s)
    {
        size[p] = s;
    }

    // Translate vertex id to partition
    intT partition_of( intT vertex_id ) const
    {
        return partition[vertex_id];
    }

    //Get the size of each partition
    intT get_size(intT i) const
    {
        return size[i];
    }
};
#endif

struct IsPart
{
    partitioner & part;
    short p;

    IsPart( partitioner & part_, short p_ ) : part( part_ ), p( p_ ) { }
    intT operator() ( intT i )
    {
        return part.partition_of(i)==p ? 1 : 0;
    }
};
