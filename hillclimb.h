#ifndef _HILL_CLIMB_H
#define _HILL_CLIMB_H

#include "rtree.h"
#include "network.h"
#include "tools.h"
#include "dagset.h"


// Network neighbourhood iterator using edit operation; in situ
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
	SPID curnode;  // current (cand.) rotation node
	SPID *ac,*bc,*cc;      // child addresses  (in leftchild/rightchild)
	SPID *ap,*bp,*cp;   // parent addresses (in parent or retparent)  

public: 

	bool next();
	void reset();

};

// TailMove - only for edges with both tree nodes
class TailMove: public EditOp
{

	SPID u; // treenode
	SPID v; // a child of unode
	SPID s;
	SPID t;
	SPID p; // reversing
	SPID q;
	
	bool moved;  // moved is true if the move was executed
	bool limittotreechild;  // true - moves are limited to stay in TC class

	bool *vreachable;

	// SPID *ac,*bc,*cc;      // child addresses  (in leftchild/rightchild)
	// SPID *ap,*bp,*cp;   // parent addresses (in parent or retparent)  

	void move(SPID u, SPID v, SPID s, SPID t);

public: 

	TailMove(bool _limittotreechild) : EditOp(), vreachable(NULL), limittotreechild(_limittotreechild) {}

	bool next();
	void reset();


};

class NetworkHCStats
{

	DagSet *dagset;

	double optcost;
	long improvements;
	double hctime;
	double mergetime;
	long steps; 
	int startingnets;
	int topnetworks; 

public:

	NetworkHCStats();
	~NetworkHCStats();

	void step() { steps++; }

	int add(Dag &n) { 		
		dagset->add(n);		
		topnetworks = dagset->size();
		return 1; 		
	}

	void setcost(double cost) { 
		improvements++;
		dagset->clear(); 
		optcost = cost;
	}

	void save(string file)
	{
		dagset->save(file);
	}
	

	void finalize()
	{
		hctime = gettime()-hctime;
	}

	void start()
	{
		hctime = gettime();
	}

	// Save dat file
	void savedat(string file);

	// Print stats
	void print(bool global=false);

	// Merge HC results
	int merge(NetworkHCStats &nhc, int printstats);
	

}; 


class HillClimb
{
protected:

	vector<RootedTree*> &genetrees;
	int verbose;	

public:

	// Initializes Gene Trees
	HillClimb(vector<RootedTree*> &gtvec, int _verbose = 0): 
		genetrees(gtvec), verbose(_verbose) {}

	// Executes hill climbing using edit operation
	// Starts from net, net is modified
	// Returns optimal cost 
	// TODO: additional info (stats, more optimal solutions, etc.) 
	double climb(EditOp &op, Network *net, int costfunc, NetworkHCStats &nhcstats);

};

#endif