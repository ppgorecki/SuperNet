#ifndef _CONTRNET_H_
#define _CONTRNET_H_


#include "network.h"



class ContractedNetwork: public Network
{

protected: 
	SPID *mapdn;
	SPID *mapup;
	RETUSAGE retdeleted; 
	DISPLAYTREEID *localbitmask;

	string _newickrepr(SPID v, SPID p);
	void _checkmaps(SPID v, SPID p);
	void _init();	
	void propagate_maps(SPID v);

	virtual bool _skiprtedge(SPID i, SPID iparent, DISPLAYTREEID id);

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
		void contractedge(SPID c, SPID p);
		void contractabove(SPID v);
		void gendotcontracted(ostream &s);
		void checkmaps();

		virtual RootedTree* gendisplaytree(DISPLAYTREEID id, RootedTree *t);

		// Reticulation count
		SPID rtcount();

		string newickrepr();

		SPID getconflictedreticulation(RETUSAGE &retusage);
		
		SPID getparent(SPID v) { return mapup[parent[v]]; }
		SPID getleftchild(SPID v) { return mapdn[leftchild[v]]; }
		SPID getrightchild(SPID v) { return mapdn[rightchild[v]]; }
		SPID getretchild(SPID v) { return mapdn[retchild[v]]; }
		  
		bool getnodeiter(SPID &i);
  
};

#endif
