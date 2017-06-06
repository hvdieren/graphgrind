#define the long 64--bit for large scale graph
INTT = -DLONG 
#INTE = -DEDGELONG

#MACHINE:M
M=jacob

#cilkview tool using the envirment
#LID= /var/shared/cilktools-4225/include/cilktools
#LID= -I$(HOME)/../../usr/include -I$(HOME)/../../usr/lib64 -L$(HOME)/../../usr/lib64

#compilers
ifdef CILK
PCC = g++
PCFLAGS = -fcilkplus -lcilkrts -O2 -DCILK $(INTT) $(INTE)
PLFLAGS = -fcilkplus -lcilkrts

else ifdef MKLROOT
PCC = icpc
PCFLAGS = -O3 -DCILKP $(INTT) $(INTE)

else ifdef OPENMP
PCC = g++
PCFLAGS = -fopenmp -O3 -DOPENMP $(INTT) $(INTE)

#here is the clang for cilk
else

#ifeq ($USER,"jsun")
SWANCCDIR=/home/jsun/swanModel/$(M)/llvm-clang/bin
SWANRTDIR=${HOME}/swanModel/$(M)/cilk-swan
#else
#SWANCCDIR=$(HOME)/work/research/llvm/bin
#SWANRTDIR=$(HOME)/work/research/llvm
#endif

PCC = $(SWANCCDIR)/clang++
PCFLAGS = -fcilkplus -lcilkrts -O3 -DCILK $(INTT) $(INTE) -I $(SWANRTDIR)/include -L $(SWANRTDIR)/lib -ldl
endif

#PCFLAGS += -I./cilkpub_v105/include
COMMON= ligra.h papi_code.h graph.h utils.h IO.h parallel.h gettime.h quickSort.h parseCommandLine.h mm.h partitioner.h graph-numa.h ligra-numa.h

ALL= BFS BC Components PageRank PageRankDelta BellmanFord SPMV BP
#LIBS_I_NEED= -lnuma -DMIX=0 -DCOMPRESSED_VERTICES=1 -DCOMPRESSED_CSCCSR=0
#LIBS_I_NEED= -lnuma -DMIX=0 -DCOMPRESSED_CSCCSR=0
#LIBS_I_NEED= -lnuma -DMIX=1 -DCOMPRESSED_CSCCSR=1
#LIBS_I_NEED= -lnuma -DMIX=0 -DCOMPRESSED_CSCCSR=1 
#LIBS_I_NEED= -lnuma -DMIX=0 -DCPU_PARTITION=0 -DVERTEX_BASED=1 -DCOMPRESSED_CSCCSR=1
LIBS_I_NEED= -lnuma -DMIX=0 -DCPU_PARTITION=1 -DVERTEX_BASED=0 -DCOMPRESSED_CSCCSR=1
#LIBS_I_NEED= -lnuma -DMIX=0 -DCPU_PARTITION=1 -DPART96=1 -DCOMPRESSED_CSCCSR=1
#LIBS_I_NEED= -lnuma -DMIX=0 -DCPU_PARTITION=1 -DVERTEX_BASED=0 -DCOMPRESSED_CSCCSR=0
all: $(ALL)

#other option
CLIDOPT += -std=c++11 
CACHEOPT += -lpapi
NUMAOPT += -DNUMA=1
#When using clang:
# OPT += -Wall -Wno-cilk-loop-control-var-modification

% : %.C $(COMMON)
	$(PCC) $(PCFLAGS) $(CFLAGS) $(CACHEOPT) $(NUMAOPT) $(OPT) -o $@ $< $(LIBS_I_NEED)
#       $(PCC) $(OPT) $(PCFLAGS) $(CFLAGS) -o $@ $< $(LIBS_I_NEED) $(NUMAOPT)
#	$(PCC) $(OPT) $(PCFLAGS) $(CFLAGS) -o $@ $< $(LIBS_I_NEED) $(NUMAOPT) $(CLIDOPT)

.PHONY : clean

clean :
	rm -f *.o $(ALL)
