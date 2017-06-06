// -*- C++ -*-
#include <stdlib.h>
#include "parallel.h"
#include "partitioner.h"
#include <assert.h>
#include <unistd.h>
#include <errno.h>
#include <string>
#include <algorithm>
#include <sys/mman.h>
#include <numa.h>
#include <numaif.h>
#include <vector>
#define MAP_HUGETLB 0x40000
//#define FLAGS (MAP_PRIVATE|MAP_ANON|MAP_HUGETLB)
#define FLAGS (MAP_PRIVATE|MAP_ANON)
#define PROTECTED (PROT_WRITE|PROT_READ)
#define bflags (MPOL_MF_MOVE)
#define mflag (MPOL_PREFERRED)
// Linux only supports two megabyte pages
//Align totalSize 1 << 21 2M 1G 1024*2048=2097152;
intptr_t page_size = intptr_t(1)<<21;
//Align page for part_allocation
intptr_t small_size = intptr_t(4)<<10;
using namespace std;
static double mmap_alloc=0;
static double del_time=0;
template <typename T>
class mmap_ptr
{
    size_t totalSize;
    void * mem;
public:
    mmap_ptr():mem(0),totalSize(0) {}

    mmap_ptr(const partitioner & part)
    {
        // Constructor intended for frontiers
        // and algorithm-specific vertex arrays
        part_allocate(part);
    }
    mmap_ptr(size_t elements)
    {
        // Constructor intended for whole graph's edge array. 
        // It does a page-by-page
        // interleaved allocation
        Interleave_allocate(elements);
    }

    mmap_ptr(size_t elements, size_t numa_node)   // NUMA-local allocation
    {
        // Constructor intended for partitioned graphs.
        local_allocate(elements,numa_node);
    }

    void part_allocate(const partitioner &part)
    {
        //timer part_alloc;
       // part_alloc.start();
        if( totalSize !=0 || mem != 0 )
        {
            cerr<<"partitioner already allocated"<<'\n';
            abort();
        }
        totalSize = part.get_num_elements()*sizeof(T);
        if((totalSize % page_size) !=0)
        {
           totalSize = (((totalSize+page_size-1)/ page_size)) * page_size;
	}
        mem = mmap( 0, totalSize, PROTECTED, FLAGS ,0, 0);
        if( mem == (void *)-1 || mem ==(void *)0 )
        {
            std::cerr << "part mmap failed: " << strerror(errno) << ", size " << totalSize << '\n';
            exit(1);
        }
        // cout << "mmap: mem=" << mem << " size=" << totalSize << "\n";
	//For frontier and algorithm data array, if not use 
	//NUMA aware allocation, not use mbind, only mmap       
        const int partNum = part.get_num_partitions();
        const int perNode = part.get_num_per_node_partitions();
        intptr_t pmem = reinterpret_cast<intptr_t>(mem);
        //Try to do the vector special allocation
        for ( int p =0 ; p < num_numa_node; ++p)
        {
            for( int i = perNode*p; i < perNode*(p+1); ++i )
            { 
                //This function use too many time during huge array 
                //to do special allocation use the mbind()
                size_t size = part.get_size(i)*sizeof(T);
                intptr_t pmem_rounded = round_page(pmem,small_size,size);
                bind_pages(reinterpret_cast<void*>(pmem_rounded),size,mflag,p);
#if 0
                cout << "mem=" << mem << " rnd=" << (void*)pmem_rounded
				 << " size=" <<size
				 << " totalSize=" <<totalSize
				 << " Page_Size=" <<page_size
				 << " end of range= " << (void*)((char*)mem+totalSize)
				 << "\n";
#endif
                pmem +=size;
            }
        }
       //mmap_alloc+=part_alloc.next();
    }

    void Interleave_allocate(size_t elements)
    {
        if( totalSize!=0 ||mem != 0 )
        {
            cerr<<"Interleave already allocated\n";
            abort();
        }

        totalSize = elements*sizeof(T);
        if((totalSize % page_size) !=0)
        {
           totalSize = (((totalSize+page_size-1)/ page_size)) * page_size;
	}
        mem = mmap( 0, totalSize, PROTECTED, FLAGS ,0, 0);
        if( mem == (void *)-1|| mem==(void*)0)
        {
            std::cerr << "numa interleave mmap failed: " << strerror(errno) << ", size " << totalSize << '\n';
            exit(1);
        }
        interleave_pages(mem,totalSize);
    }

    void local_allocate(size_t elements, int numa_node)
    {
        if( totalSize!=0 ||mem != 0 )
        {
            cerr<<"local NUMA already allocated\n";
            abort();
        }
        totalSize = elements*sizeof(T);
        if((totalSize % page_size) !=0)
        {
           totalSize = (((totalSize+page_size-1)/ page_size)) * page_size;

        }
        mem = mmap( NULL, totalSize, PROTECTED,FLAGS, 0, 0);
        if( mem == (void *)-1|| mem==(void*)0)
        {
            std::cerr << "numa-node mmap failed: " << strerror(errno) << ", size " << totalSize << '\n';
            exit(1);
        }
        bind_pages(mem,totalSize,mflag,numa_node);
    }

private:
    intptr_t round_page(intptr_t ptr, intptr_t pagesize, size_t & sz)
    {
        intptr_t ret = (ptr + (pagesize-1)) & ~intptr_t(pagesize - 1);
        assert( (ret & (pagesize - 1)) == 0 );

	intptr_t size = sz;
	size += (pagesize-1) - ((ptr + pagesize-1) & intptr_t(pagesize -1));
        size = (size + pagesize-1) & ~intptr_t(pagesize - 1);
        assert( (size & (pagesize - 1)) == 0 );
	sz = size;
        return ret;
    }

    void bind_pages(void * mem, size_t size,int policy, int numa_node)
    {
        struct bitmask *bmp;
        bmp = numa_allocate_nodemask();
        numa_bitmask_setbit(bmp, numa_node);
        if (mem == (void *)-1)
            mem = NULL;
        else
            dombind(mem, size, policy, bmp);
        numa_bitmask_free(bmp);
    }

    void interleave_pages(void * mem, size_t size)
    {
        struct bitmask *bmp;
        bmp = numa_allocate_nodemask();
        numa_bitmask_setall(bmp);
        if (mem == (void *)-1)
            mem = NULL;
        else
            dombind(mem, size, MPOL_INTERLEAVE, bmp);
        numa_bitmask_free(bmp);
    }

    static void dombind(void *mem, size_t size, int policy, struct bitmask *bmp)
    {
       if(mbind(mem, size, policy, bmp ? bmp->maskp : NULL, bmp ? bmp->size : 0, 0 )<0)
            std::cerr << "mbind failed: " << strerror(errno) 
		      << " address " << (size_t*)mem
		      << ", size " << size << '\n';
    }

public:
    void del()
    {
        timer del;
        del.start();
        if( mem )
        {
            int munmapres = munmap (mem,totalSize);
            if(munmapres == -1)
            {
                cerr<<"munmap failed "<<errno<<" "<<strerror(errno)
                    <<" address "<<mem
                    <<" and size "<<totalSize<<endl;
                abort();
            }
        }
        mem = 0;
        totalSize=0;
        del_time+=del.next();
    }

    operator T * ()
    {
        return reinterpret_cast<T *>(mem);
    }
    operator const T * () const
    {
        return reinterpret_cast<T *>(mem) ;
    }

    operator bool () const
    {
        return mem != 0;
    }

    size_t get_bytes() const
    {
        return totalSize;
    }

    T *get() const
    {
        return reinterpret_cast<T *>( mem );
    }
};

