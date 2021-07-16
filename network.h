

#ifndef _NETWORK_H_
#define _NETWORK_H_

#include "dag.h"
#include "rtree.h"

// Bit represents reticulation switching
// Range 0..2^{rt-1}, where rt is the number of reticulation nodes
// Current limit is rt<=64

typedef unsigned long DISPLAYTREEID;


void initbitmask();

class Network: public Dag
{

	DISPLAYTREEID displaytreemaxid;

	SPID _gendisplaytree(DISPLAYTREEID id, RootedTree *t, SPID i, SPID iparent, SPID &freeint);

public:
	
	Network(char *s, double weight=1.0): Dag(s,weight) 
	{
		if (rt > 8*sizeof(DISPLAYTREEID)) 
	    {
	      cout << "Network has too many reticulation nodes (" << rt << "). The limit is " << 8*sizeof(DISPLAYTREEID) << "." << endl;    
	      exit(-1);
	    }
		displaytreemaxid = 1 << rt;		
	}

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
};

#endif