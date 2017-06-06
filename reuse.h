// -*- C++ -*-

#ifndef REUSE_H
#define REUSE_H

#include <fstream>
#include <cstdio>
#include <cstdlib>


template<typename IntTy>
class RDLogger {
    struct node {
	node *l;
	node *r;
	node *p;
	size_t timestamp;
	size_t volume;
    };
    
    IntTy n;

    // Old variables
    // size_t lru_fill;
    // IntTy *lru;
    // bool *present;

    // New variables
    node *root;
    node *prealloc;
    size_t timestamp;

    static const size_t absent = 0;

    size_t *freq;
    size_t nsteps;

public:
    RDLogger(size_t n_, size_t nsteps_)
	: n(n_), // lru_fill(0),
	  root(nullptr), timestamp(0), nsteps(nsteps_) {
	freq = new size_t[n]();
	// lru = new IntTy[n];
	// present = new bool[n]();

	prealloc = new node[n];
	for( size_t i=0; i < n; ++i )
	    prealloc[i].timestamp = absent;
    }
    ~RDLogger() {
	delete[] freq;
	// delete[] present;
	// delete[] lru;

	delete[] prealloc;
    }

    void access(IntTy idx) {
	if( nsteps == 0 )
	    return;

	assert( 0 <= idx && idx < n );

#if 0
	size_t oldrd = ~size_t(0);

	if( present[idx] ) {
	    // std::cout << "fill=" << lru_fill << ", access " << idx << "\n";
	    size_t pos = 0;
	    IntTy c = idx;
	    while( pos < lru_fill ) {
		// std::cout << "LRU[" << pos << "] = " << lru[pos] << "\n";
		if( lru[pos] == idx ) {
		    freq[pos]++;  // reuse distance pos observed
		    oldrd = pos;
		    lru[pos] = c; // update final LRU entry: pos moved to top
		    // std::cout << idx << " RD=" << pos << " freq " << freq[pos] << "\n";
		    break; // currently same as return
		} else {
		    std::swap( lru[pos], c ); // push LRU contents down
		}
		++pos;
	    }

	    // Not found on stack. No reuse distance recorded
	    if( pos == lru_fill ) {
		lru[lru_fill++] = c;
		// std::cout << idx << " RD=inf\n";
	    }
	} else {
	    for( size_t i=lru_fill; i > 0; --i )
		lru[i] = lru[i-1];
	    lru[0] = idx;
	    ++lru_fill;
	    present[idx] = true;
	    // std::cout << idx << " RD=inf\n";
	}
#endif

	size_t last = prealloc[idx].timestamp;
	node * p = &prealloc[idx];
	++timestamp;
	
	if( p->timestamp != absent ) {
	    size_t rd = FindRD( p->timestamp );
	    // assert( rd == oldrd );
	    freq[rd]++;
	    Erase( p ); // Includes Splay(p);
	} else
	    ; // assert( oldrd == ~size_t(0) );
	
	p->l = nullptr;
	p->r = nullptr;
	p->p = nullptr;
	p->timestamp = timestamp;
	p->volume = 1;
	Insert(p);
    }

    void step() {
	// if( lru_fill > 0 ) {
	if( root ) {
	    std::cout << "******************** STEP/" << /*lru_fill*/ calculateVolume(root) << " ********************\n";
	    // std::fill( &present[0], &present[n], false );
	    parallel_for( size_t i=0; i < n; ++i )
		prealloc[i].timestamp = absent;
	    // lru_fill = 0;
	    root = 0;

	    if( nsteps > 0 )
		--nsteps;
	}
    }

    void dump(const char *fname) {
	std::ofstream f(fname);
	dump(f);
    }
    
    ostream & dump(ostream & os) {
	for( size_t i=0; i < n; ++i )
	    os  << freq[i] << '\n';
	return os;
    }

private:

    size_t calculateVolume(node *P) {
	return P ? (1 + (P->l ? P->l->volume : 0) + (P->r ? P->r->volume : 0))
	    : 0;
    }
    
    // Based on http://codeforces.com/blog/entry/18462
    // extended with calculation of volume of nodes in sub-tree
    void rightRotate(node *P) {
	node *T=P->l;
	node *B=T->r;
	node *D=P->p;
	if(D) {
	    if(D->r==P) D->r=T;
	    else D->l=T;
	}
	if(B)
	    B->p=P;
	T->p=D;
	T->r=P;
	
	P->p=T;
	P->l=B;

	// D->volume = unchanged
	// B->volume = unchanged
	P->volume = calculateVolume(P);
	T->volume = calculateVolume(T);
    }

    void leftRotate(node *P) {
	node *T=P->r;
	node *B=T->l;
	node *D=P->p;
	if(D) {
	    if(D->r==P) D->r=T;
	    else D->l=T;
	}
	if(B)
	    B->p=P;
	T->p=D;
	T->l=P;
	
	P->p=T;
	P->r=B;

	// D->volume = unchanged
	// B->volume = unchanged
	P->volume = calculateVolume(P);
	T->volume = calculateVolume(T);
    }

    void Splay(node *T) {
	while(true) {
	    node *p=T->p;
	    if(!p) break;
	    node *pp=p->p;
	    if(!pp)//Zig
	    {
		if(p->l==T)
		    rightRotate(p);
		else
		    leftRotate(p);
		break;
	    }
	    if(pp->l==p)
	    {
		if(p->l==T)
		{//ZigZig
		    rightRotate(pp);
		    rightRotate(p);
		}
		else
		{//ZigZag
		    leftRotate(p);
		    rightRotate(pp);
		}
	    }
	    else
	    {
		if(p->l==T)
		{//ZigZag
		    rightRotate(p);
		    leftRotate(pp);
		}
		else
		{//ZigZig
		    leftRotate(pp);
		    leftRotate(p);
		}
	    }
	}
	root=T;
    }

    void Insert(node *new_node) {
	if(!root) {
	    root = new_node;
	    return;
	}
	node *P=root;
	while(true) {
	    if( P->timestamp == timestamp ) break; // not multiset
	    if( timestamp < (P->timestamp) ) {
		if(P->l)
		    P=P->l;
		else
		{
		    P->l = new_node;
		    P->l->p=P;
		    P->l->r=nullptr;
		    P->l->l=nullptr;
		    // P->l->v=v;
		    P=P->l;
		    break;
		}
	    }
	    else
	    {
		if(P->r) P=P->r;
		else
		{
		    P->r = new_node;
		    P->r->p=P;
		    P->r->r=NULL;
		    P->r->l=NULL;
		    // P->r->v=v;
		    P=P->r;
		    break;
		}
	    }
	}
	Splay(P);
    }

    size_t FindRD(size_t timestamp) {
	// Note: only called if item is present in tree
	if(!root) return 0; // should not occur
	node *P=root;
	size_t volume = 0;
	while(P) {
	    if(P->timestamp == timestamp) {
		if( P->r )
		    volume += P->r->volume;
		break;
	    }
	    if(timestamp<(P->timestamp)) {
		if(P->l) {
		    volume += P->volume - P->l->volume;
		    P=P->l;
		} else {
		    if( P->r )
			volume += P->r->volume;
		    break;
		}
	    } else {
		if(P->r)
		    P=P->r;
		else
		    break;
	    }
	}

	Splay(P); // P == &prealloc[v];

	return volume;
    }

/*
    node* Find(IntTy v) {
	// Could also use direct access provided that we initialise
	// all elements of pre-alloc, and reset them on step(), to indicate
	// absence from tree.
	if(!root) return NULL;
	node *P=root;
	while(P) {
	    if(P->v==v)
		break;
	    if(v<(P->v)) {
		if(P->l)
		    P=P->l;
		else
		    break;
	    } else {
		if(P->r)
		    P=P->r;
		else
		    break;
	    }
	}
	Splay(P);
	if(P->v==v) return P;
	else return NULL;
    }
*/

    bool Erase(node *N) {
	if(!N) return false;
	Splay(N); //check once more; now root == N
	node *P=N->l;
	if(!P) {
	    root = N->r;
	    if( root )
		root->p = NULL;
	    N->timestamp = absent; // dealloc
	    return true;
	}
	while(P->r) {
	    --P->volume; // will move one node out of sub-tree
	    P = P->r;
	}
	if(N->r) {
	    P->r=N->r;
	    N->r->p=P;
	}
	root=N->l;
	root->p=NULL;
	root->volume = calculateVolume(root);
	N->timestamp = absent; // dealloc
	return true;
    }

};

#endif // REUSE_H
