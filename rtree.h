#ifndef _RTREE_H
#define _RTREE_H

#include <string.h>

// If LCATAB set lca queries are stored in O(n^2) array
// Otherwise LCA queries are computed directly using depth
// TODO: test which one is more efficient
// TODO: lca queries in O(1)

// LCATAB seems not helpful when display trees are often used
// #define LCATAB 

#include "tools.h"
#include "clusters.h"
#include "dag.h"
#include "network.h"
#include "bitcluster.h"


class Network;
class CostFun;

class RootedTree: public Dag
{

protected:

  // Set depth of nodes 
  void _setdepth(NODEID i, int dpt);



  // Insert clusters of children node n into TreeClusters
  GTCluster* _getgtccluster(NODEID n, TreeClusters *gtc);

  // Used in sorted representation to compute min leaf labels
  SPID _setmm_repr(NODEID i, SPID* mm);

	// Used in sorted representation of tree
	void _repr(SPID *t, NODEID i, int &pos, SPID* mm);

	void inittreestr();

	// Variant of _cost without precomputing needed structures
	double _cost(RootedTree &speciestree, NODEID* lcamap, CostFun &costfun);

	bool depthinitialized;

#ifdef LCATAB
	bool lcatabinitialized;
#endif



public:

	NODEID *depth;
	GTCluster **gtcclusters; // node clusters

	int treeid;

#ifdef LCATAB
  NODEID **lcatab;  
  NODEID _lca(NODEID a, NODEID b);
#endif

  // Gen tree from str
	RootedTree(char *s, double weight=1.0);

	RootedTree(string s, double weight=1.0): RootedTree(strdup(s.c_str()), weight) { }

	// Gen caterpillar tree from labs; cherry from the first two labels
	RootedTree(int _lf, NODEID* labels, double weight=1.0); 

	virtual ~RootedTree()
  	{
 	   if (depth) delete[] depth; 	   
	}

  // Init depth array, i.e., distance from the root where depth[root]=0
	void initdepth();

	// Get all clusters of nodes
	// Now used only to compute the root cluster in quasiconsensus 
	// TODO replace it by GTClusters
	NODEID** getspclusterrepr();

	// Compute cost using species tree
	double cost(RootedTree &speciestree, CostFun &cost);

	// Computes ODT via naive approach
	// Returns the cost and id of an optimal display tree via optid
	// double odtnaivecost(Network &network, CostFun &costfun, DISPLAYTREEID &optid);
	

 //  // Compute duplication cost 
 // 	long costduplication(RootedTree &speciestree, NODEID *lcamap);

	// // Compute deep-coalescence cost 
 // 	long costdeepcoalx(RootedTree &speciestree, NODEID *lcamap);

 // 	// Compute duplication-loss cost 
 // 	long costduplicationloss(RootedTree &speciestree, NODEID *lcamap);

 // 	// Compute Robinson-Foulds cost 
 // 	long costrobinsonfoulds(RootedTree &speciestree, NODEID *lcamap);

 // 	// Compute loss cost 
 // 	long costloss(RootedTree &speciestree, NODEID *lcamap);
  
 	// Insert clusters into Tree Clusters
 	void setspclusters(TreeClusters *gtc);

 	// all clusters 
 	bitcluster* getbitclusterrepr();

 	// Compute lca mapping to the species tree
 	// Allocates array if lcamap = NULL
  NODEID* getlcamapping(RootedTree &speciestree, NODEID* lcamap = NULL);

  // Compute lca between nodes a and b
 	NODEID lca(NODEID a, NODEID b);

 	// Initialize lca query computations
 	void initlca();

 	// Generate unique (sorted) representation of a tree 
  SPID* repr(SPID *t = NULL);

  void setid(int id) { treeid = id; }
  int getid() { return treeid; }
 
  ostream& printsubtrees(ostream &s)
  {
  	for (NODEID i=0;i<nn;i++)
  		printsubtree(s,i) << endl;
  	return s;
  }

  friend class Network;

  virtual ostream& printdeb(ostream&s, int gse, string tn="");

  friend class CFDuplication;
  friend class CFDeepCoalescence;
  friend class CFDeepCoalescenceEdge;
  friend class CFLoss;
  friend class CFDuplicationLoss;
};


// Gen random species tree from defined set of species 
RootedTree *randspeciestree();

string randspeciestreestr();

#endif