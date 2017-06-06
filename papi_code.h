#include <cilk/cilk.h>
#include <cilk/cilk_api.h>
#include <cilk/reducer.h>
long long **values;
#include <papi.h>
#define NUM_EVENTS 4
int native = 0x0;
int n_threads;
bool init=false;
bool *init_start;
int retval;
int *EventSet;
char local_DRAM[]= "OFFCORE_RESPONSE_0:ANY_REQUEST:LLC_MISS_LOCAL:SNP_NONE:SNP_NOT_NEEDED:SNP_MISS:SNP_NO_FWD:u=0:k=0";
char remote_DRAM[]="OFFCORE_RESPONSE_1:ANY_REQUEST:LLC_MISS_REMOTE:SNP_NONE:SNP_NOT_NEEDED:SNP_MISS:SNP_NO_FWD:u=0:k=0";
char ins_count[]="INSTRUCTION_RETIRED:u=0:k=0";
char TLB[]="perf::PERF_COUNT_HW_CACHE_DTLB:MISS";
int event_codes[NUM_EVENTS];
#define START 1
#define STOP 2

/*----------------------------Functions related to PAPI event counting------------------------------------*/
void on_all_workers_help( size_t i, size_t n, volatile bool * flags, int stage )
{
    if(n>1)
    {
        if( i < n-1 )
        {
            cilk_spawn on_all_workers_help( i+1, n, flags, stage);
        }
    }
    int id= __cilkrts_get_worker_number();
    
    if(n>1)
    {
        if( i == n-2 )
        {
            for (int j=0; j<n; j++)
                flags[j]=true;
        }
        else
        {
           // printf("started busy waiting\n");
            while( !flags[i]); // busy wait
            //printf("ended busy wait\n");
        }
    }
    if(stage == START)
    {
        if(!init_start[id]){
            init_start[id]=true;
           if(PAPI_create_eventset(&EventSet[id])!= PAPI_OK)
              printf("error creating eventset: id = %d \n", id);
           if (PAPI_add_events(EventSet[id], event_codes, NUM_EVENTS) != PAPI_OK)
            printf("event add error: id %d \n", id);
       }
        if(PAPI_start(EventSet[id])!=PAPI_OK) 	//call papi_start(eventset), call papi_stop(eventset, values)
            printf("start error: id =%d \n", id);

    }

    if (stage==STOP)
    {
        if(PAPI_stop(EventSet[id], values[id])!=PAPI_OK)
            printf("stop error: iter: %lu Thread: %d \n", i, id);
    }

    cilk_sync;
}
/*---------------------------------------------------------------------------*/
void start_on_all_workers()
{
    int n = __cilkrts_get_nworkers();
    n_threads=n;
   volatile bool flags[n];
    for( int i=0; i < n; ++i )
    {
        flags[i] = false;
    }
    if(!init){
       init=true;
       init_start= new bool[n];
       EventSet= new int[n];
       values = new long long*[n];
       for( int i=0; i < n; ++i )
       {
           EventSet[i]= PAPI_NULL;
           values[i] = new long long[NUM_EVENTS];
           init_start[i]=false;
        }
    }
    on_all_workers_help( 0, n, flags, START);
}
/*---------------------------------------------------------------------------*/
void stop_on_all_workers()
{
    int n = __cilkrts_get_nworkers();
    volatile bool flags_stop[n];
    for( int i=0; i < n; ++i )
        flags_stop[i] = false;
    on_all_workers_help( 0, n, flags_stop, STOP );
}
/*---------------------------------------------------------------------------*/
inline void PAPI_start_count()
{
    start_on_all_workers( );
}
/*---------------------------------------------------------------------------*/
inline void PAPI_stop_count()
{
    stop_on_all_workers( );
}

/*---------------------------------------------------------------------------*/

/*-------------------------------------------------MAIN()----------------------------------------------*/
/*-----------------------------------------------------------------------------------------------------*/
inline void PAPI_initial()
{

    /* PAPI library initialisations*/
    retval = PAPI_library_init(PAPI_VER_CURRENT);
    if (retval != PAPI_VER_CURRENT)
    {
        printf("PAPI library init error! \n");
        exit(1);
    }
    if (PAPI_thread_init(pthread_self) != PAPI_OK)
    {
        printf("thread init error\n");
        exit(1);
    }
    /*------------Map PAPI event names to hardware event codes--------*/
#if 0
    if (PAPI_event_name_to_code(ins_count,&event_codes[0]) != PAPI_OK)
    {
        printf("event 3 name error");
        exit(1);
    }
#endif
    if (PAPI_event_name_to_code(local_DRAM,&event_codes[0]) != PAPI_OK)
    {
        printf("event 1 name error\n");
        exit(1);
    }
    if (PAPI_event_name_to_code(remote_DRAM,&event_codes[1]) != PAPI_OK)
    {
        printf("event 2 name error");
        exit(1);
    }
    if (PAPI_event_name_to_code(ins_count,&event_codes[2]) != PAPI_OK)
    {
        printf("event 3 name error");
        exit(1);
    }
    if (PAPI_event_name_to_code(TLB,&event_codes[3]) != PAPI_OK)
    {
        printf("event 4 name error");
        exit(1);
    }
    /*-----------------------------------------------------------------*/
}
inline void PAPI_print()
{

    /*-------- print values of PAPI counters for all threads------*/
    printf("Threads	L3_MISS_Local	L3_MISS_REMOTE	Instruction_count TLB\n");
#if THREADS_CACHE
    for (int k=0; k<n_threads; k++)
    {
    printf("%d\t ", k);   /*L3_MISS*/
    printf("%llu\t", values[k][0]);   /*L3_MISS*/
    printf("%llu\t", values[k][1]);   /*TLB_MISS*/
    printf("%llu\t", values[k][2]);   /*Ins_MISS*/
    printf("%llu\n", values[k][3]);   /*TLB_MISS*/
    printf("\n");
    }
#else
    long long TLB_total=0;
    long long Ins_total=0;
    long long L3_local=0;
    long long L3_remote=0;
    for (int k=0; k<n_threads; k++)
    {
#if 0
          Ins_total+=values[k][0];
#endif
          L3_local+=values[k][0];
          Ins_total+=values[k][2];
          L3_remote+=values[k][1];
          TLB_total+=values[k][3];
    }
    printf("%llu\t", L3_local);   /*L3_MISS*/
    printf("%llu\t", L3_remote);   /*TLB_MISS*/
    printf("%llu\t", Ins_total);   /*Ins_MISS*/
    printf("%llu\n", TLB_total);   /*TLB_MISS*/
    printf("\n");
#endif
    for (int k=0; k<n_threads; k++)
    {
          values[k][0]=0;
          values[k][2]=0;
          values[k][1]=0;
          values[k][3]=0;
	}

    //delete [] values;
}
inline void PAPI_end(){
    delete [] init_start;
    delete [] EventSet;
    delete [] values;
}

