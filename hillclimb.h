#ifndef _HILL_CLIMB_H
#define _HILL_CLIMB_H

#include "rtree.h"
#include "network.h"
#include "tools.h"
#include "costs.h"
#include "netgen.h"

#define ODTBASENAME "odt"

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
	double climb(EditOp &op, Network *net, CostFun &costfun, 
		ClimbStats &nhcstats, 
		bool usenaive, 
		int runnaiveleqrt_t, 		
		int hcmaximprovements = 0,
		int hcstopclimb = 0,
		float displaytreesampling=0,
		bool cutwhendtimproved=false,
		bool stopatcostdefined=false,
		float stopatcost=0
		);

};

void supnetheuristic(		
		vector<RootedTree*> &gtvec,				
		NetGen *netgenerator,
		EditOp *op,
		CostFun *costfun,

		int printstats,				
		int hcstopinit,
		int hcstopclimb,
		bool usenaive, 
		int runnaiveleqrt_t, 		
		int hcmaximprovements,		
		vector<ClimbStatsGlobal*> globalstatsrarr,
		bool cutwhendtimproved,
		int multipleoptima=0,
		bool stopatcostdefined=false,
		float stopatcost=0,
		DagScoredQueue* scoreddags=NULL // do not report
		);

void iterativeretinsertionoptimizer(		
		vector<RootedTree*> &gtvec,				
		Network *startnet,		
		CostFun *costfun,		
		int printstats,						
		bool usenaive, 
		int runnaiveleqrt_t, 				
		ClimbStatsGlobal *globalstats, // could be sampler
		bool cutwhendtimproved,
		int networkclass, 
		int timeconsistency, 
		Clusters *guideclusters=NULL, 
		Clusters *guidetree=NULL, 
		string _retid=""
		);

#endif