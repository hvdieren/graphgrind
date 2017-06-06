# include <stdio.h>
# include <unistd.h>
# include <math.h>
# include <float.h>
# include <limits.h>
# include <sys/time.h>

# include <errno.h>
# include <string.h>

#if NUMA

#include <numaif.h>
#include <numa.h>
#define ALIGN(p,a) (void*)(((long)p)&~(long)((a)-1))
template<class OT,class intT>
int numa_page_check(OT* a, intT sample_size,intT node ){
	
        long pagesize = sysconf( _SC_PAGESIZE );
	//printf ("System page size is %ld KiB\n", pagesize/1024);

	long array_size = sample_size * sizeof(a[0]);
	long npages = (array_size + pagesize-1) / pagesize;

	void **page_list = (void **)malloc(sizeof(void *)*npages);
	int *status = (int *)malloc(sizeof(int *)*npages);

	long stride = pagesize / sizeof(a[0]);
	long p, e;
	int i;	

        for( i=0, p=0; i < npages; ++i, p += stride ) {
	    page_list[i] = ALIGN(&a[p], pagesize);
	    status[i] = 0;
	}

	e = numa_move_pages(0, npages, page_list, NULL, status, MPOL_MF_MOVE);

	if( e != 0 )
	   printf( "move_pages validation error on a[]: %s\n", strerror(errno) );
	
        //printf( "a[0]=%p a[last]=%p\n", (void*)&a[0], (void*)&a[sample_size] );
        
        for( i=0; i < npages; ++i ) {
            //if (status[i]!=node) {printf("Now in node %d,not %d\n",status[i],node);}//exit(1);
            if (status[i]!=node) exit(1);
	    //printf( "page %d of a[] on node %d, req=%p\n", i, status[i], page_list[i] );
	    if( status[i] < 0 )
		exit(1);
	}
     return 0;
}
#else // ! NUMA
template<class OT,class intT>
int numa_page_check(OT* a, intT sample_size,intT node ){
    return 0;
}
#endif // NUMA
