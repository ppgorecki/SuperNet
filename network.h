

#ifndef _NETWORK_H_
#define _NETWORK_H_

#include "dag.h"
#include "rtree.h"

#define NT_GENERAL 2   // no limits
#define NT_CLASS1 1    // int node has >=1 tree node/leaf child
#define NT_TREECHLD 0  // int node has <=1 ret. child

// Bit represents reticulation switching
// Range 0..2^{rt-1}, where rt is the number of reticulation nodes
// Current limit is rt<=64

typedef unsigned long DISPLAYTREEID;


void initbitmask();

class Network: public Dag
{

	DISPLAYTREEID displaytreemaxid;

	SPID _gendisplaytree(DISPLAYTREEID id, RootedTree *t, SPID i, SPID iparent, SPID &freeint);

	virtual ostream& printdebstats(ostream&s);

	void _getreachablefrom(SPID v, bool *reachable, bool *visited);
	void _getreachableto(SPID v, bool *reachable, bool *visited);

	void initdid() 
	{
		if (rt > 8*sizeof(DISPLAYTREEID)) 
	    {
	      cout << "Network has too many reticulation nodes (" << rt << "). The limit is " << 8*sizeof(DISPLAYTREEID) << "." << endl;    
	      exit(-1);
	    }
		displaytreemaxid = 1 << rt;		
	}

public:
	
	Network(char *s, double weight=1.0): Dag(s,weight) { initdid(); }

	Network(string s, double weight=1.0): Dag(s,weight) { initdid(); }

	Network(Dag *d, SPID v, SPID p, SPID w, SPID q, string retid, double dagweight=1.0):
		Dag(d, v, p, w, q, retid, dagweight)
	{
		initdid();
	}

	// Copy network and add random reticulation node
	// retid - reticulation label; if empty then the label will be automatically assigned to #NUM
	//     if r==NULL then the label will be "#NUM" where NUM is the number of reticulations
	// type - 0 - TC insertion
	// type - 1 - non-TC type 1
	// type - 2 - general
	// 
	// Returns NULL if such a network cannot be created
	

	Network *addrandreticulation(string retid, int networktype);
	
	// Generate id'th display tree; 
	// id encodes switching reticulation nodes and should be in range 0..2^{rt-1}	
	// id encodes reticulation schema -> for long size is 8 -> 64 bits
	// If t is given the display tree will be generated using t,
	// otherwise new RootedTree is allocated
	// Returns NULL if there is no display tree 
	// Typical usage:
	//   RootedTree *t = NULL;
	//   DISPLAYTREEID tid = 0;
	//   while ((t=n->gendisplaytree(tid,t))!=NULL)   { ... use display tree t ... }    
	RootedTree* gendisplaytree(DISPLAYTREEID id, RootedTree *t = NULL);

	// Compute min cost vs. gene trees
	double odtnaivecost(vector<RootedTree*> &genetrees, int costfunc, DISPLAYTREEID &optid);

	// Mark nodes reachable from v (including v)	
	void getreachablefrom(SPID v, bool *reachable);

	// Mark nodes w such that v is reachble from w (including v)	
	void getreachableto(SPID v, bool *reachable);

	friend class NNI;
	friend class TailMove;

};

#endif