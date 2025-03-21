
#ifndef _NETWORK_H_
#define _NETWORK_H_

#include "dag.h"
#include "bb.h"
#include "rtree.h"
#include "treespace.h"
#include "stats.h"

class ContractedNetwork;
class CostFun;
class BBTreeStats;


#ifndef NODTCACHE
	#define DTCACHE
#endif

#ifndef MAXRTNODES
	#define MAXRTNODES 64
#endif 


class RootedTree;
class Network;

typedef vector<RootedTree*> VecRootedTree;
typedef vector<Network*> VecNetwork;

class AdaptiveBB;


// Todo: Optimize bit operators

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
// Current limit is rt==64
typedef unsigned long DISPLAYTREEID;

void initbitmask();
struct SNode;
class TreeSpace;

class Network: public Dag
{

protected:

	NODEID _gendisplaytree(DISPLAYTREEID id, RootedTree *t, NODEID i, NODEID iparent, NODEID &freeint);
	struct SNode *_gendisplaytree2(DISPLAYTREEID id, SNode *t, NODEID i, NODEID iparent, TreeSpace *tsp);

	virtual ostream& printdebstats(ostream&s);

	void _getreachablefrom(NODEID v, bool *reachable, bool *visited);
	void _getreachableto(NODEID v, bool *reachable, bool *visited);

	virtual bool _skiprtedge(NODEID i, NODEID iparent, DISPLAYTREEID id);
	virtual void initdid();

	void _getdestconnectorsforedge(int networkclass, NODEID v, int &dsrclen, NODEID* dsrc);
	int _getsrcconnectors(int networkclass, NODEID *esrc);
	int _getdestconnectors(int networkclass, int len, NODEID *esrc, NODEID* dsrc, bool uniform);

public:
	
	Network(char *s, double weight=1.0): Dag(s,weight) { initdid(); }
	Network(string s, double weight=1.0): Dag(s,weight) { initdid(); }
	Network(Dag *d, NODEID v, NODEID p, NODEID w, NODEID q, string retid, double dagweight=1.0):
		Dag(d, v, p, w, q, retid, dagweight)
	{
		initdid();
	}

	Network(Dag &dag, bool shallowcopy): Dag(dag, shallowcopy) { initdid(); }

	NODEID tchild_contract(uint8_t mark[], NODEID leftc[], NODEID rightc[]);

	// Max id + 1 of a display tree 
	DISPLAYTREEID displaytreemaxid() { return 1 << rtcount(); }
	
	// Copy network and add random reticulation node
	// retid - reticulation label; if empty then the label will be automatically assigned to #NUM
	//     if r==NULL then the label will be "#NUM" where NUM is the number of reticulations
	// type - 0 - TC insertion
	// type - 1 - non-TC type 1
	// type - 2 - general

	// uniform=true - each pair of candidate edges is equally probable
	// uniform=false - take random src edge -> take random reachable dest edge -> connect
	// 
	// Returns NULL if such a network cannot be created

	Network *addrandreticulation(string retid, int networkclass, bool uniform=false);
		
	
	// Generate id'th display tree; 
	// id encodes switching reticulation nodes and should be in range 0..2^{rt-1}	
	// id encodes reticulation schema -> for long size is 8 -> 64 bits
	// If t is given the display tree will be generated using t,
	// otherwise new RootedTree is allocated
	// Returns NULL if there is no display tree 
	// Typical usage:
	//   RootedTree *t = NULL;
	//   DISPLAYTREEID tid = 0;


	virtual RootedTree* gendisplaytree(DISPLAYTREEID id, RootedTree *t);
	virtual SNode* gendisplaytree2(DISPLAYTREEID id, SNode *t, TreeSpace *tsp);

	virtual ~Network() {}

	// Compute min cost vs. gene trees via enumeration of all display trees
	double odtcostnaive(vector<RootedTree*> &genetrees, CostFun &costfun, struct ODTStats &odtstats, float sampling = 0,
		bool cutwhendtimproved = false,
		double externalbestcost = 0);

	// Compute min cost vs. gene tree via enumeration of all display trees
	double odtcostnaive(RootedTree *genetree, CostFun &costfun, ODTStats &odtstats, float sampling = 0);
	
	double odtcostnaiverev(
		vector<RootedTree*> &genetrees, 
		CostFun &costfun, 
		ODTStats &odtstats, 	
		bool cutwhendtimproved, 
		double externalbestcost);
	
	// Compute min cost vs. gene trees via DP&BB alg.	
	// Only for DC
	double odtcostdpbb(vector<RootedTree*> &genetrees, CostFun &costfun, int runnaiveleqrt_t, ODTStats &odtstats);

	// Compute min cost vs. gene trees 	
	double odtcost(vector<RootedTree*> &genetrees, CostFun &costfun, bool usenaive_oe, int runnaiveleqrt_t, ODTStats &odtstats, float displaytreesampling = 0,
		bool cutwhendtimproved = false,
		double externalbestcost = 0);
	
	// Mark nodes reachable from v (including v)	
	void getreachablefrom(NODEID v, bool *reachable);

	// Mark nodes w such that v is reachble from w (including v)	
	void getreachableto(NODEID v, bool *reachable);	

	// Return cluster = visible leaf labels 
	NODEID *getcluster(NODEID v);

	// approx DCE(G,N) via DP
	COSTT approxmindce(RootedTree &genetree, CostFun &costfun);
	COSTT approxmindceusage(RootedTree &genetree, RETUSAGE &retusage, CostFun &costfun);

	// exact DCE(G,M) via BB
	COSTT mindce(RootedTree &genetree, 
		int runnaiveleqrt_t, 
		CostFun &costfun,
		ODTStats &odtstats,
		BBTreeStats *bbtreestats=NULL, 
		COSTT bbstartscore=0,
		bool bbstartscoredefined=false);

	
	// nodetype==1 -> count visible leaves 
	// nodetype==2 -> count all visible nodes
	ostream& visibilenodestats(int nodetypes, ostream&s);

	friend class NNI;
	friend class TailMove;	

	friend class DP_DC;
	friend class NetworkRetIterator;

};


class NetworkRetIterator
{

	protected:
		Network &net;

		int pair; 
		NODEID *esrc;   
		NODEID *dsrc;   
		int esrclen;
		int currentpair;
		int currentesrc;
		int dsrclen;
		string retid;
		int networkclass;
		int timeconsistency;
		Clusters *guidetree;
		Clusters *guideclusters;


	public: 
		NetworkRetIterator(Network &n, int _networkclass, int _timeconsistency, Clusters *_guideclusters=NULL, Clusters *_guidetree=NULL, string _retid="") : 
			net(n), pair(0), retid(_retid), networkclass(_networkclass), timeconsistency(_timeconsistency), guideclusters(_guideclusters), guidetree(_guidetree)
		{
	 		int len=net.nn+net.rt; 
	 		esrc = new NODEID[len];        
	 		esrclen = net._getsrcconnectors(networkclass, esrc);	
	 		dsrc = new NODEID[len*2]; // cand pairs	
			currentpair = 0;
			currentesrc = -1;		
			dsrclen = 0;						
		}
		
		Network *_next()
		{					
			while (dsrclen==currentpair)
			{
				currentesrc++;
				if (currentesrc>=esrclen) return NULL;
				dsrclen = 0;
				net._getdestconnectorsforedge(networkclass, esrc[currentesrc], dsrclen, dsrc);	
				currentpair = 0;		  		
			}	

			NODEID p, v = dsrc[currentpair];
			NODEID q, w = dsrc[currentpair+1];		
			currentpair+=2;

#define getvencnet(nod,par) if (nod==net.root) { par=MAXNODEID; } else  { if (nod<0) { nod=-nod; par = net.retparent[nod]; } else { par = net.parent[nod]; } }

			getvencnet(v,p);		
			getvencnet(w,q);

			return new Network(&net, v, p, w, q, retid);						
		}

		Network *next()
		{
			Network *n;
			
			while ((n = _next())!=NULL)
			{
				if (!n->checktimeconsistency(timeconsistency)) continue;
				if (guideclusters && !n->hasclusters(guideclusters)) continue;
				if (guidetree && !n->hastreeclusters(guidetree)) continue;
				return n;
			}

			return NULL;

		}
		
		~NetworkRetIterator() { 
			delete[] esrc; 
			delete[] dsrc; 
		}
	
};

#endif

