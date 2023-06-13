#ifndef _HILL_CLIMB_H
#define _HILL_CLIMB_H

#include "rtree.h"
#include "network.h"
#include "tools.h"
#include "costs.h"


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
	bool limittotreechild;  // true - moves are limited to stay in TC class

	bool *vreachable;

	// NODEID *ac,*bc,*cc;      // child addresses  (in leftchild/rightchild)
	// NODEID *ap,*bp,*cp;   // parent addresses (in parent or retparent)  

	void move(NODEID u, NODEID v, NODEID s, NODEID t);

public: 

	TailMove(bool _limittotreechild) : EditOp(), vreachable(NULL), limittotreechild(_limittotreechild) {}

	bool next();
	void reset();


};



class HillClimb
{
protected:

	vector<RootedTree*> &genetrees;

public:

	// Initializes Gene Trees
	HillClimb(vector<RootedTree*> &gtvec): 
		genetrees(gtvec) {}

	// Executes hill climbing using edit operation
	// Starts from net, net is modified
	// Returns optimal cost 
	// TODO: additional info (stats, more optimal solutions, etc.) 
	double climb(EditOp &op, Network *net, CostFun &costfun, NetworkHCStats &nhcstats, bool usenaive_oe, int runnaiveleqrt_t);

};

#endif