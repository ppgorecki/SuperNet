

#ifndef _NETWORK_H_
#define _NETWORK_H_

#include "dag.h"
#include "bb.h"
#include "rtree.h"

class ContractedNetwork;
class CostFun;

#define NT_GENERAL 2   // no limits
#define NT_CLASS1 1    // int node has >=1 tree node/leaf child
#define NT_TREECHLD 0  // int node has <=1 ret. child

#ifndef MAXRTNODES
	#define MAXRTNODES 64
#endif 

class AdaptiveBB;

#if MAXRTNODES < 33

	typedef uint64_t RETUSAGE;
	#define _allleft 6148914691236517205   // bits    ...101010101
	#define _allright 12297829382473034410   // bits  ...010101010

	#define emptyretusage(r) r=0
	#define unionretusage(r,a,b) r=((a)|(b))

	// rtids: 0..rt-1

	#define addleftretusage(r, rtid)  r|=(1<<((rtid)*2))    // left -> via parent
	#define addrightretusage(r, rtid)  r|=(1<<((rtid)*2+1)) // right -> via retparent

	#define leftret(r, rtid)  ((r)&(1<<((rtid)*2))) 
	#define rightret(r, rtid)  ((r)&(1<<((rtid)*2+1))) 

	#define conflicted(r) (((((r)&_allright)>>1) & ((r)&_allleft))!=0)

#elif MAXRTNODES < 65

	typedef struct { uint64_t lft, rgh; } RETUSAGE;

	#define emptyretusage(r) { r.lft=0; r.rgh=0; }
	#define unionretusage(r,a,b) { r.rgh=a.rgh|b.rgh; r.lft=a.lft|b.lft; }

	// rtids: 0..rt-1
	#define addleftretusage(r, rtid)  r.lft|=1<<(rtid)   
	#define addrightretusage(r, rtid)  r.rgh|=1<<(rtid)  

	#define leftret(r, rtid)  ((r.lft)&(1<<(rtid)))
	#define rightret(r, rtid)  ((r.rgh)&(1<<(rtid)))

	#define conflicted(r) ((r.lft & r.rgh)!=0)

 	ostream& operator<<(ostream& os, const RETUSAGE& r); 

#else

 	// TODO more generic retusages

#endif

// Bit represents reticulation switching
// Range 0..2^{rt-1}, where rt is the number of reticulation nodes
// Current limit is rt<=64
typedef unsigned long DISPLAYTREEID;

void initbitmask();

class Network: public Dag
{

protected:

	SPID _gendisplaytree(DISPLAYTREEID id, RootedTree *t, SPID i, SPID iparent, SPID &freeint);

	virtual ostream& printdebstats(ostream&s);

	void _getreachablefrom(SPID v, bool *reachable, bool *visited);
	void _getreachableto(SPID v, bool *reachable, bool *visited);

	virtual bool _skiprtedge(SPID i, SPID iparent, DISPLAYTREEID id);
	virtual void initdid();

public:
	
	Network(char *s, double weight=1.0): Dag(s,weight) { initdid(); }
	Network(string s, double weight=1.0): Dag(s,weight) { initdid(); }
	Network(Dag *d, SPID v, SPID p, SPID w, SPID q, string retid, double dagweight=1.0):
		Dag(d, v, p, w, q, retid, dagweight)
	{
		initdid();
	}

	Network(Network &net, int shallowcopy) : Dag(net,shallowcopy) { initdid(); }

	SPID tchild_contract(uint8_t mark[], SPID leftc[], SPID rightc[]);

	// Max id + 1 of a display tree 
	DISPLAYTREEID displaytreemaxid() { return 1 << rtcount(); }
	
	// Copy network and add random reticulation node
	// retid - reticulation label; if empty then the label will be automatically assigned to #NUM
	//     if r==NULL then the label will be "#NUM" where NUM is the number of reticulations
	// type - 0 - TC insertion
	// type - 1 - non-TC type 1
	// type - 2 - general
	// 
	// Returns NULL if such a network cannot be created
	Network *addrandreticulation(string retid, int networktype, bool unform=false);
		
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


	virtual RootedTree* gendisplaytree(DISPLAYTREEID id, RootedTree *t);

	virtual ~Network() {}

	// Compute min cost vs. gene trees via enumeration of all display trees
	double odtcostnaive(vector<RootedTree*> &genetrees, CostFun &costfun);

	// Compute min cost vs. gene tree via enumeration of all display trees
	double odtcostnaive(RootedTree *genetree, CostFun &costfun);

	// Compute min cost vs. gene trees via DP&BB alg.	
	// Only for DC
	double odtcostdpbb(vector<RootedTree*> &genetrees, CostFun &costfun, int runnaiveleqrt);

	// Compute min cost vs. gene trees 	
	double odtcost(vector<RootedTree*> &genetrees, CostFun &costfun, bool usenaive, int runnaiveleqrt);
	
	// Mark nodes reachable from v (including v)	
	void getreachablefrom(SPID v, bool *reachable);

	// Mark nodes w such that v is reachble from w (including v)	
	void getreachableto(SPID v, bool *reachable);	

	// approx DCE(G,N) via DP
	COSTT approxmindce(RootedTree &genetree, CostFun &costfun);
	COSTT approxmindceusage(RootedTree &genetree, RETUSAGE &retusage, CostFun &costfun);

	// exact DCE(G,M) via BB
	COSTT mindce(RootedTree &genetree, 
		int runnaiveleqrt, 
		CostFun &costfun,
		BBTreeStats *bbtreestats=NULL, 
		COSTT bbstartscore=0,
		bool bbstartscoredefined=false);
	
	// nodetype==1 -> count visible leaves 
	// nodetype==2 -> count all visible nodes
	ostream& visibilenodestats(int nodetypes, ostream&s);

	friend class NNI;
	friend class TailMove;	
	friend class DP_DC;

};


#endif

