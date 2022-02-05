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

#define COSTDUPLICATIONLOSS 0   // dl
#define COSTDUPLICATION 1    // duplication
#define COSTDEEPCOAL 2   // deep coalescence
#define COSTROBINSONFOULDS 3   // robinson-foulds
#define COSTLOSS 4    // loss

class Network;

class RootedTree: public Dag
{

protected:

  // Set depth of nodes 
  void _setdepth(SPID i, int dpt);



  // Insert clusters of children node n into TreeClusters
  GTCluster* _getgtccluster(SPID n, TreeClusters *gtc);

  // Used in sorted representation to compute min leaf labels
  SPID _setmm_repr(SPID i, SPID* mm);

	// Used in sorted representation of tree
	void _repr(SPID *t, SPID i, int &pos, SPID* mm);

	void inittreestr();

	// Variant of _cost without precomputing needed structures
	double _cost(RootedTree &speciestree, SPID* lcamap, int costfunc);

	bool depthinitialized;

#ifdef LCATAB
	bool lcatabinitialized;
#endif



public:

	SPID *depth;
	GTCluster **gtcclusters; // node clusters

#ifdef LCATAB
  SPID **lcatab;  
  SPID _lca(SPID a, SPID b);
#endif

  // Gen tree from str
	RootedTree(char *s, double weight=1.0);

	RootedTree(string s, double weight=1.0): RootedTree(strdup(s.c_str()), weight) { }

	// Gen caterpillar tree from labs; cherry from the first two labels
	RootedTree(int _lf, SPID* labels, double weight=1.0); 

	virtual ~RootedTree()
  	{
 	   if (depth) delete[] depth; 	   
	}

  // Init depth array, i.e., distance from the root where depth[root]=0
	void initdepth();

	// Get all clusters of nodes
	// Now used only to compute the root cluster in quasiconsensus 
	// TODO replace it by GTClusters
	SPID** getspclusterrepr();

	// Compute cost using species tree
	double cost(RootedTree &speciestree, int costfunc);

	// Computes ODT via naive approach
	// Returns the cost and id of an optimal display tree via optid
	double odtnaivecost(Network &network, int costfunc, DISPLAYTREEID &optid);
	

  // Compute duplication cost 
 	long costduplication(RootedTree &speciestree, SPID *lcamap);

	// Compute deep-coalescence cost 
 	long costdeepcoal(RootedTree &speciestree, SPID *lcamap);

 	// Compute duplication-loss cost 
 	long costduplicationloss(RootedTree &speciestree, SPID *lcamap);

 	// Compute Robinson-Foulds cost 
 	long costrobinsonfoulds(RootedTree &speciestree, SPID *lcamap);

 	// Compute loss cost 
 	long costloss(RootedTree &speciestree, SPID *lcamap);
  
 	// Insert clusters into Tree Clusters
 	void setspclusters(TreeClusters *gtc);

 	// Compute lca mapping to the species tree
 	// Allocates array if lcamap = NULL
  SPID* getlcamapping(RootedTree &speciestree, SPID* lcamap = NULL);

  // Compute lca between nodes a and b
 	SPID lca(SPID a, SPID b);

 	// Initialize lca query computations
 	void initlca();

 	// Generate unique (sorted) representation of a tree as a vector
  SPID* repr(SPID *t = NULL);

  ostream& printrepr(ostream&s);

  friend class Network;

  virtual ostream& printdeb(ostream&s, int gse, string tn="");

  
};


// Gen random species tree from defined set of species 
RootedTree *randspeciestree();

string randspeciestreestr();

#endif