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
#define VERTEX 1
#if NUMA
#include "ligra-numa.h"
#else
#include "ligra.h"
#endif

typedef float fType;

struct BC_F
{
    fType* NumPaths;
    bool* Visited;
    static const bool use_cache = true;
    struct cache_t
    {
        fType value;
    };
    BC_F(fType* _NumPaths, bool* _Visited) :
        NumPaths(_NumPaths), Visited(_Visited) {}
    inline bool update(intT s, intT d)  //Update function for forward phase
    {
        fType oldV = NumPaths[d];
        NumPaths[d] += NumPaths[s];
        return oldV == (fType)0;
    }
    inline bool updateAtomic (intT s, intT d)   //atomic Update, basically an add
    {
        volatile fType oldV, newV;
        do
        {
            oldV = NumPaths[d];
            newV = oldV + NumPaths[s];
        }
        while(!CAS(&NumPaths[d],oldV,newV));
        return oldV == (fType)0;
    }
    inline void create_cache(cache_t &cache, intT d)
    {
        cache.value = NumPaths[d];
    }
    inline bool update(cache_t &cache, intT s)
    {
        fType oldV = cache.value;
        cache.value += NumPaths[s];
        return oldV == (fType)0;
    }

    inline void commit_cache(cache_t &cache, intT d)
    {
	// Cache used only when vertex accessed sequentially
	NumPaths[d] = cache.value;
    }
    inline bool cond (intT d)
    {
        return Visited[d] == false;    //check if visited
    }
};

struct BC_Back_F
{
    fType* Dependencies;
    bool* Visited;
    BC_Back_F(fType* _Dependencies, bool* _Visited) :
        Dependencies(_Dependencies), Visited(_Visited) {}
    static const bool use_cache = true;//true;
    struct cache_t
    {
        fType value;
    };
    inline bool update(intT s, intT d)  //Update function for backwards phase
    {
        fType oldV = Dependencies[d];
        Dependencies[d] += Dependencies[s];
        return oldV == (fType)0;
    }
    inline bool updateAtomic (intT s, intT d)   //atomic Update
    {
        volatile fType oldV, newV;
        do
        {
            oldV = Dependencies[d];
            newV = oldV + Dependencies[s];
        }
        while(!CAS(&Dependencies[d],oldV,newV));
        return oldV == (fType)0;
    }
    inline void create_cache(cache_t &cache, intT d)
    {
        cache.value = Dependencies[d];
    }
    inline bool update(cache_t &cache, intT s)
    {
        fType oldV = cache.value;
        cache.value += Dependencies[s];
        return oldV == (fType)0;
    }

    inline void commit_cache(cache_t &cache, intT d)
    {
        Dependencies[d] = cache.value;
    }
    inline bool cond (intT d)
    {
        return Visited[d] == false;    //check if visited
    }
};

//vertex map function to mark visited vertexSubset
struct BC_Vertex_F
{
    bool* Visited;
    BC_Vertex_F(bool* _Visited) : Visited(_Visited) {}
    inline bool operator() (intT i)
    {
        Visited[i] = true;
        return true;
    }
};

//vertex map function (used on backwards phase) to mark visited vertexSubset
//and add to Dependencies score
struct BC_Back_Vertex_F
{
    bool* Visited;
    fType* Dependencies, *inverseNumPaths;
    BC_Back_Vertex_F(bool* _Visited, fType* _Dependencies, fType* _inverseNumPaths) :
        Visited(_Visited), Dependencies(_Dependencies), inverseNumPaths(_inverseNumPaths) {}
    inline bool operator() (intT i)
    {
        Visited[i] = true;
        Dependencies[i] += inverseNumPaths[i];
        return true;
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

    mmap_ptr<fType> NumPaths;
    NumPaths.part_allocate (part);
    mmap_ptr<bool> Visited;
    Visited.part_allocate (part);

    loop(j,part,perNode,NumPaths[j]=0.0);
    loop(j,part,perNode,Visited[j]=false);

    NumPaths[start] = 1.0;
    Visited[start] = true;

    //vertexSubset Frontier(n,start);
    partitioned_vertices Frontier=partitioned_vertices::create(n,start,GA.getWholeGraph().V[start].getOutDegree());

    vector<partitioned_vertices> Levels;
    Levels.push_back(Frontier);

    intT round = 0;
    timer t1,t2;
    while(!Frontier.isEmpty())  //first phase
    {
        round++;
        partitioned_vertices output=edgeMap(GA,Frontier,BC_F(NumPaths,Visited),m/20);
        vertexMap(part,output, BC_Vertex_F(Visited)); //mark visited
        Levels.push_back(output); //save frontier onto Levels
        Frontier = output;
    }
    mmap_ptr<fType> Dependencies;
    Dependencies.part_allocate (part);
    
    loop(j,part,perNode,Dependencies[j]=0.0);

    //invert numpaths
    mmap_ptr<fType> inverseNumPaths;
    inverseNumPaths = NumPaths;
    loop(j,part,perNode,inverseNumPaths[j]=1/inverseNumPaths[j]);

    Levels[round].del();
    //reuse Visited
    loop(j,part,perNode, Visited[j]=false);
    Frontier = Levels[round-1];
    vertexMap(part,Frontier,BC_Back_Vertex_F(Visited,Dependencies,inverseNumPaths));

    //tranpose graph
    GA.transpose();

    for(intT r=round-2; r>=0; r--) //backwards phase
    {

        partitioned_vertices output=edgeMap(GA,Frontier,BC_Back_F(Dependencies,Visited), m/20);
        output.del();
        Frontier.del();
        Frontier = Levels[r]; //gets frontier from Levels array
        //vertex map to mark visited and update Dependencies scores
        vertexMap(part,Frontier,BC_Back_Vertex_F(Visited,Dependencies,inverseNumPaths));
    }

    Frontier.del();
    //Update dependencies scores
    loop(j,part,perNode, Dependencies[j]=(Dependencies[j]-inverseNumPaths[j])/inverseNumPaths[j]);
    inverseNumPaths.del(); //free(inverseNumPaths);
    Visited.del(); //free(Visited);
    Dependencies.del(); //free(Dependencies);
}
