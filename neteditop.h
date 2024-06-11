#ifndef _NETEDITOPT_H_
#define _NETEDITOPT_H_

// #include "rtree.h"
#include "network.h"
// #include "tools.h"
// #include "costs.h"
// #include "randnets.h"


// Network neighbourhood iterator using edfit operation; in situ
// next - for accessing next neigbour of a network; returns false -> no more neighbours
// reset - sets the current network to be the source
class EditOp
{	

protected:
	Network *source;

public:  	

	EditOp() {}

	// Initialize 
	virtual void init(Network *net);

	// Get next neighbour; returns false if there is none
	virtual bool next()=0;	

	// Search in the neighbourhood of the current network
	virtual void reset()=0;

};

// NNI - only for edges with both tree nodes
class NNI: public EditOp
{

	int cycle;     // cycle=0 -> search for next; cycle=1,2,3 - in the rotation of curnode
	NODEID curnode;  // current (cand.) rotation node
	NODEID *ac,*bc,*cc;      // child addresses  (in leftchild/rightchild)
	NODEID *ap,*bp,*cp;   // parent addresses (in parent or retparent)  

public: 

	bool next();
	void reset();

};

// TailMove - only for edges with both tree nodes
class TailMove: public EditOp
{

	NODEID u; // treenode
	NODEID v; // a child of unode
	NODEID s;
	NODEID t;
	NODEID p; // reversing
	NODEID q;
	
	bool moved;  // moved is true if the move was executed
	int networkclass;  // treechild, relaxed or general
	bool timeconsistent; // if networks have to be time consistent

	bool *vreachable;
	Clusters *guideclusters;	
	Clusters *guidetree;	

	void move(NODEID u, NODEID v, NODEID s, NODEID t);

public: 

	TailMove(int _networkclass, Clusters *_guideclusters=NULL, Clusters *_guidetree=NULL) : EditOp(), vreachable(NULL), networkclass(_networkclass), guideclusters(_guideclusters), guidetree(_guidetree) {}

	bool next();
	void reset();

};

#endif