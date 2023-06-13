#ifndef _CONTRNET_H_
#define _CONTRNET_H_


#include "network.h"

class ContractedNetwork: public Network
{

protected: 
	NODEID *mapdn;
	NODEID *mapup;
	RETUSAGE retdeleted; 
	DISPLAYTREEID *localbitmask;

	string _newickrepr(NODEID v, NODEID p);
	void _checkmaps(NODEID v, NODEID p);
	void _init();	
	void propagate_maps(NODEID v);
	void initlocalbitmask(DISPLAYTREEID id);

	virtual bool _skiprtedge(NODEID i, NODEID iparent, DISPLAYTREEID id);

	public:

		ContractedNetwork(char *s, double weight=1.0): Network(s,weight) { _init(); }
		ContractedNetwork(string s, double weight=1.0): Network(s,weight) { _init(); }		
		
		// Create copy ContractedNetwork from a Network	
		// node arrays are not duplicated if shallow is true
		ContractedNetwork(Network &net, int shallowcopy) : Network(net,shallowcopy) { _init(); }	

		// Copy of contracted network
		// Note that mapdn/mapup are always allocated		
		ContractedNetwork(ContractedNetwork &net, int shallowcopy);

		virtual ~ContractedNetwork();

		void contract(RETUSAGE &retcontract);
		void contract2(RETUSAGE &retcontract);
		void gendot(ostream &s);
		void contractedge(NODEID c, NODEID p);
		void contractabove(NODEID v);
		void gendotcontracted(ostream &s);
		void checkmaps();

		virtual RootedTree* gendisplaytree(DISPLAYTREEID id, RootedTree *t);
		virtual SNode* gendisplaytree2(DISPLAYTREEID id, SNode *t, TreeSpace *tsp);

		// Reticulation count
		NODEID rtcount();

		string newickrepr();

		NODEID getconflictedreticulation(RETUSAGE &retusage);
		
		NODEID getparent(NODEID v) { return mapup[parent[v]]; }
		NODEID getleftchild(NODEID v) { return mapdn[leftchild[v]]; }
		NODEID getrightchild(NODEID v) { return mapdn[rightchild[v]]; }
		NODEID getretchild(NODEID v) { return mapdn[retchild[v]]; }
		  
		bool getnodeiter(NODEID &i);

		ostream& printdebstats(ostream&s);
  
};

#endif
