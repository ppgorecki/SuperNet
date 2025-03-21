	
#include <random>

#include "rtree.h"
#include "network.h"
#include "bb.h"
#include "dp.h"
#include "costs.h"
#include "treespace.h"

// #define _PRINTNAIVESTATS_

extern int verbosealg;

DISPLAYTREEID bitmask[8*sizeof(long)]; 

void initbitmask()
{
	bitmask[0]=1;
	for (int i=1; i<8*sizeof(long); i++)
		bitmask[i]=2*bitmask[i-1];
}

ostream& Network::printdebstats(ostream&s)
{
	return Dag::printdebstats(s) << " displaytreemaxid=" << displaytreemaxid();
}


// Generate id'th display tree
// id encodes reticulation schema -> for long size is 8 -> 64 bits
RootedTree* Network::gendisplaytree(DISPLAYTREEID id, RootedTree *t)
{

	if (id>=displaytreemaxid()) return NULL;
	if (!t)	
		t = new RootedTree(lf, lab);				
	
	NODEID freeint = lf;
	// cerr << "==============" << id << " " << rtcount() << endl; 
	t->root = _gendisplaytree(id, t, getroot(), MAXNODEID, freeint); 
	t->parent[t->root] = MAXNODEID;
	t->depthinitialized = false;

	if (t->root!=t->size()-1)
	{		
		// t->printdeb(cerr, 0, "") ;
		cerr << "Root Error? Expected size " << t->size()-1 << ". Found " << t->root << " getroot=" << getroot() << endl;
		cerr << *this << endl;
		cerr << "tid=" << id << endl;
		// cerr << *t << endl;
		exit(-1);
		return NULL;
	}	
	return t;
}


// Generate id'th display tree
// id encodes reticulation schema -> for long size is 8 -> 64 bits
SNode* Network::gendisplaytree2(DISPLAYTREEID id, SNode *t, TreeSpace *treespace)
{
	if (id>=displaytreemaxid()) return NULL;

	SNode *dtroot = _gendisplaytree2(id, t, getroot(), MAXNODEID, treespace); 

	return dtroot;
}

SNode *Network::_gendisplaytree2(DISPLAYTREEID id, SNode *t, NODEID i, NODEID iparent, TreeSpace *treespace) 
{ 
	
	// leaf
	if (i < lf) { 		
		return treespace->leaf(lab[i]); 
	}

	// skip current edge
	if (i>=rtstartid)
	{				
		if (_skiprtedge(i, iparent, id)) return NULL;		
	}

	NODEID cleft = MAXNODEID;
	getchild(i, cleft);

	
	SNode* cleftgen = _gendisplaytree2(id, t, cleft, i, treespace);	
	
	NODEID cright = cleft;
	if (!getchild(i, cright))
	{
		// reticulation only
		return cleftgen;  // maybe MAXNODEID (ignore)
	}

	SNode *crightgen = _gendisplaytree2(id, t, cright, i, treespace);

	if (!cleftgen) return crightgen;  // maybe MAXNODEID
	if (!crightgen) return cleftgen;  // single node

	return treespace->find(cleftgen,crightgen);
}

bool Network::_skiprtedge(NODEID i, NODEID iparent, DISPLAYTREEID id)
{
	return ((iparent == retparent[i]) == bool(id & bitmask[i-rtstartid]));
}

NODEID Network::_gendisplaytree(DISPLAYTREEID id, RootedTree *t, NODEID i, NODEID iparent, NODEID &freeint) 
{ 
	
	// leaf
	if (i < lf) return i; 

	// skip current edge
	if (i>=rtstartid)
	{		
		//cout << "RT check!" << i << " " << rtstartid << endl;
		if (_skiprtedge(i, iparent, id)) return MAXNODEID;		
	}
		

	NODEID cleft = MAXNODEID;
	getchild(i,cleft);
	NODEID cleftgen = _gendisplaytree(id, t, cleft, i, freeint);

	// cout << (int)i << " cleft" << (int)cleft << " GD-LEFT " << (int)cleftgen << endl;

	NODEID cright = cleft;
	if (!getchild(i,cright))
	{
		// reticulation only
		return cleftgen;  // maybe MAXNODEID (ignore)
	}

	NODEID crightgen = _gendisplaytree(id, t, cright, i, freeint);

	if (cleftgen==MAXNODEID) return crightgen;  // maybe MAXNODEID
	if (crightgen==MAXNODEID) return cleftgen;  // single node

	// two nodes to be connected
	
	t->parent[cleftgen] = t->parent[crightgen] = freeint;
	t->leftchild[freeint] = cleftgen;
	t->rightchild[freeint] = crightgen;

	return freeint++;
}

// Marks nodes reachable from v (including v)	
void Network::getreachablefrom(NODEID v, bool *reachable)
{
	bool visited[nn];
	for (NODEID i=0; i<nn; i++) reachable[i] = visited[i] = false;
	_getreachablefrom(v, reachable, visited);
}

// Print how many nodes are reachable from all nodes
// supnet -n '((((c)#B,b))#A,(#A,(#B,a)))' -eL
// 1 1 1 4 4 8 9 2 5 
// supnet -n '((a,b),c)' -eL 
// 1 1 1 3 5
// supnet -n '((a,b),c);(a,b)' -eL
// 1 1 1 3 5 
// 1 1 3 
// nodetype==1 -> count visible leaves 
// nodetype==2 -> count all visible nodes

ostream& Network::visibilenodestats(int nodetypes, ostream&s)
{
	bool *vreachable = new bool[nn];	
	int last = nn;
	if (nodetypes==1) last=lf;
	for (NODEID i=0; i<nn; i++)		
	{	
		getreachablefrom(i, vreachable);
		int cnt = 0;
		for (NODEID j=0; j<last; j++)		
			if (vreachable[j]) cnt++;
		s << cnt << " ";
	}
	return s;

}

void Network::_getreachablefrom(NODEID v, bool *reachable, bool *visited)
{
	if (visited[v]) return;
	reachable[v] = true;
	if (v<lf) return;
	if (v>=rtstartid) _getreachablefrom(retchild[v], reachable, visited);
	else
	{	
		_getreachablefrom(leftchild[v],reachable, visited);
		_getreachablefrom(rightchild[v],reachable, visited);
	}
}

// Marks nodes w such that v is reachable from w (including v)	
void Network::getreachableto(NODEID v, bool *reachable)
{
	bool visited[nn];
	for (NODEID i=0; i<nn; i++) reachable[i] = visited[i] = false;
	_getreachableto(v,reachable, visited);
}

void Network::_getreachableto(NODEID v, bool *reachable, bool *visited)
{
	if (visited[v]) return;
	reachable[v] = true;
	if (v==root) return;
	_getreachableto(parent[v], reachable, visited);
	if (v>=rtstartid)  
		_getreachableto(retparent[v],reachable, visited);		
}

double Network::odtcost(vector<RootedTree*> &genetrees, CostFun &costfun, bool usenaive_oe, int runnaiveleqrt_t, ODTStats &odtstats, 
	float displaytreesampling,
	bool cutwhendtimproved,
	double externalbestcost
	)
{	
	double res;

	if ((displaytreesampling>0) || usenaive_oe || (runnaiveleqrt_t>rtcount()))
	{		
		//if (displaytreesampling>0)
		res = odtcostnaive(genetrees, costfun, odtstats, displaytreesampling, cutwhendtimproved, externalbestcost);
		// else
		// 	res = odtcostnaiverev(genetrees, costfun, odtstats, cutwhendtimproved, externalbestcost);


		if (verbosealg>=5) 
		{ 
			cout << "ODT-cost-naive " << genetrees.size() << " tree(s);"
				 << "rt=" << rtcount();

			if (!usenaive_oe && runnaiveleqrt_t>rtcount())
				cout << "(-t)";
			else
				cout << "(usenaive)";

			cout << " cost=" << res
				<< endl;
		}
	}
	else
	{
		res = odtcostdpbb(genetrees, costfun, runnaiveleqrt_t, odtstats);
		if (verbosealg>=5) 
		{ 
			cout << "ODT-cost-BB " << genetrees.size() << " tree(s)"
				 << "; rt=" << rtcount()	
				 << " cost=" << res
				 << endl;		
		}
	}
	return res;
}

// only for DC and DCE
double Network::odtcostdpbb(vector<RootedTree*> &genetrees, CostFun &costfun, int runnaiveleqrt_t, ODTStats &odtstats)
{
	int ct = costfun.costtype();
	if (ct!=COSTDEEPCOAL && ct!=COSTDEEPCOALEDGE)
	{
		cerr << "DP&BB cost computation only for DC or DCE (use -CDC or -CDCE)" << endl;
		exit(-1);
	}

	BBTreeStats bbstats;

	double cost = 0;
	if (ct==COSTDEEPCOAL)
	{
		cost = 2*genetrees.size(); // adjust DCE to DC
		for (int gt=0; gt<genetrees.size(); gt++)
			cost+=mindce(*genetrees[gt], runnaiveleqrt_t, costfun, odtstats, &bbstats)-2*genetrees[gt]->lf; // adjust DCE to DC
	}
	else
	{
		for (int gt=0; gt<genetrees.size(); gt++)
			cost+=mindce(*genetrees[gt], runnaiveleqrt_t, costfun, odtstats, &bbstats); 
	}
    return cost;        
}

/*
 *  Compute cost via naive enumeration with optional displaytreesampling for a single gene tree
 */ 
double Network::odtcostnaive(RootedTree *genetree, CostFun &costfun, ODTStats &odtstats, 
	float displaytreesampling)
{
	vector<RootedTree*> genetrees;
	genetrees.push_back(genetree);
	return odtcostnaive(genetrees, costfun, odtstats, displaytreesampling);
}

extern bool print_repr_inodtnaive;

/*
	-profiler
	-cache D.Trees vs. costs, wczesniej test czestoci uzycia d.trees	
	-lca queries w czasie stałym po liniowym preprocessing (C/C++?)
*/

#ifdef _PRINTNAIVESTATS_
clock_t dtgen =0, dtgen2 =0, lcagen = 0, depthgen = 0,  lcamapgen=0, costgen=0;
clock_t treespacegen = 0, treespaceinit=0, btreprgen=0;
long odtcnt = 0;
#endif

long gcnt=0;

extern TreeSpace *globaltreespace;

/*
 *  Compute cost via naive enumeration with optional displaytreesampling 
 */ 

std::random_device rd;
std::mt19937 gen(rd());


double Network::odtcostnaiverev(
	vector<RootedTree*> &genetrees, 
	CostFun &costfun, 
	ODTStats &odtstats, 	
	bool cutwhendtimproved, 
	double externalbestcost)
{
    RootedTree tr(lf, lab);	
    RootedTree *t = &tr;
    DISPLAYTREEID tid = 0; // id of display tree
    double mincost = 0;
    NODEID *lcamaps[genetrees.size()];
    double gtcost;
    double lb;
    
    odtstats.startnaive();

    // https://en.cppreference.com/w/cpp/numeric/random/exponential_distribution
        
    
    
    // cout << "LBN:" << externalbestcost << " :" << endl;
    // for (int gt=0; gt<genetrees.size(); gt++)
    // 	cout << lb[gt] << " ";    
    // cout << endl;
   
	int lbcnt = 0;
	
#ifdef _PRINTNAIVESTATS_
	clock_t start, end;
	double stoplbratio = 1.0;

	#define STARTCLOCK start = clock();
	#define TERMCLOCK(var) var += clock() - start;
#else
	#define STARTCLOCK
	#define TERMCLOCK(var)
#endif  

	SNode *s = NULL;

	long dtgtcnt = 0;

	for (int gt=0; gt<genetrees.size(); gt++)
    {    		
    	RootedTree *genetree = genetrees[gt];
    	lb = costfun.lowerboundnet(*genetree, *this);    			    	    	
    	gtcost = lb - 1;
    	int treeid = genetree->getid();
	#ifdef DTCACHE
 		while (1) 
    	{         

		
			STARTCLOCK;
    		if ((s = gendisplaytree2(tid, s, globaltreespace))==NULL) break;    
    		TERMCLOCK(dtgen2); 
    		odtstats.displaytreecnt++;
    		COSTT curcost = s->cost[treeid] - 2*genetree->lf+2;

    		if (gtcost>curcost)
    		{
    			gtcost = curcost;
    			if (gtcost == lb)
    				break; // lower bound reached
    		}
    	}
    	mincost += gtcost;
    			

	
	#else
		exit(-1);
	#endif
       
    } // for 

#ifdef DTCACHE
    	globaltreespace->treecompleted(); // inform treespace 
#endif    	

	odtstats.stopnaive();

	if (verbosealg>=6)
	{
		cout << " --odtcostnaive cost=" << mincost 
			 << " dt/gt_cost_pairs=" << dtgtcnt << endl;
	}
    return mincost;        
}




double Network::odtcostnaive(vector<RootedTree*> &genetrees, 
	CostFun &costfun, 
	ODTStats &odtstats, 
	float displaytreesampling, 
	bool cutwhendtimproved, 
	double externalbestcost)
{
    RootedTree tr(lf, lab);	
    RootedTree *t = &tr;
    DISPLAYTREEID tid = 0; // id of display tree
    double mincost = 0;
    NODEID *lcamaps[genetrees.size()];
    double gtcost[genetrees.size()];
    double lb[genetrees.size()];

	std::exponential_distribution<> expdistr(displaytreesampling);
    	
    
    odtstats.startnaive();

    // https://en.cppreference.com/w/cpp/numeric/random/exponential_distribution
        
    for (int gt=0; gt<genetrees.size(); gt++)
    {    		
    	RootedTree *genetree = genetrees[gt];
    	lb[gt]=costfun.lowerboundnet(*genetree, *this);    			    	    	
    	gtcost[gt] = 0;
    }   
    
    // cout << "LBN:" << externalbestcost << " :" << endl;
    // for (int gt=0; gt<genetrees.size(); gt++)
    // 	cout << lb[gt] << " ";    
    // cout << endl;
   
	int lbcnt = 0;
	
#ifdef _PRINTNAIVESTATS_
	clock_t start, end;
	double stoplbratio = 1.0;

	#define STARTCLOCK start = clock();
	#define TERMCLOCK(var) var += clock() - start;
#else
	#define STARTCLOCK
	#define TERMCLOCK(var)
#endif  

	SNode *s = NULL;

	long dtgtcnt = 0;

    while (1) 
    {         

#ifdef DTCACHE
		STARTCLOCK;
    	if ((s = gendisplaytree2(tid, s, globaltreespace))==NULL) break;    
    	TERMCLOCK(dtgen2); 
#else
		STARTCLOCK;
    	if ((t = gendisplaytree(tid,t))==NULL) break;            

    	TERMCLOCK(dtgen);         	
    	STARTCLOCK;
    	t->initlca();
    	TERMCLOCK(lcagen);    	
    	STARTCLOCK;
    	t->initdepth(); 				    	
    	TERMCLOCK(depthgen);
		
#endif  
		odtstats.displaytreecnt++;

    	// Compute cost of the current tree vs all gene trees
    	lbcnt = 0;
    	for (int gt=0; gt<genetrees.size(); gt++)
    	{   	
    	 	// Check if the lower bound is already reached		
			if (tid && (gtcost[gt]==lb[gt]))
			{
				lbcnt++;
				continue; // skip, lower bound is reached
			}
			RootedTree *genetree = genetrees[gt];   

			dtgtcnt++;

#ifdef DTCACHE			
    		COSTT curcost = s->cost[genetree->getid()] - 2*genetree->lf+2;
#else    		
    		 
			NODEID *lcamap = lcamaps[gt];
    		STARTCLOCK;	
    		if (tid) 
    		{
    			    				    	    			
    			// TODO: optimize resuing leaf-maps (force exactreespaceecies==1)
    			genetree->getlcamapping(*t,lcamaps[gt]);    // overwrite previous    		
    		}
    		else  
    		{	    		
    			lcamaps[gt] = genetrees[gt]->getlcamapping(*t); // init lcamaps	
    		}
    		TERMCLOCK(lcamapgen);    			
   
    		STARTCLOCK;
    		COSTT curcost = costfun.compute(*genetree, *t, lcamaps[gt]);
    		TERMCLOCK(costgen);
#endif     		

    		if (!tid || (gtcost[gt] > curcost)) 
    		{ 
    			mincost += (curcost - gtcost[gt]);
    			
    			gtcost[gt] = curcost;  

    		// 	if (tid>0)
    		// 	{
    		// 	cout << "LBN: tid=" << tid << " ext=" << externalbestcost << " cur=" << mincost << "  -->";
    		// 	for (int gt=0; gt<genetrees.size(); gt++)
    		// 		cout << lb[gt] << ":" << gtcost[gt] << " ";    
    		// 	cout << endl;
    		// }

    			if (curcost==lb[gt])  // check again lb
    			{
#ifdef _PRINTNAIVESTATS_    				
    				stoplbratio -= (1.0-1.0*tid/displaytreemaxid())/genetrees.size();
#endif    				
    				lbcnt++;
    				continue;
    			}
    		}

    		if (print_repr_inodtnaive)
    		{    			
    			t->printrepr() << " " << curcost << endl;
    		}
    	}

#ifdef DTCACHE
    	globaltreespace->treecompleted(); // inform treespace 
#endif    	

    	if (lbcnt == genetrees.size())
    	{
    		// all is lower bound (ODT optimum)    		
    		break;
    	}
    	// DT Sampler
    	if (displaytreesampling)
    	{
    		int nexttid=(int)floor(tid+expdistr(gen));
    		if (nexttid==tid) tid++;
    		else tid = nexttid;

    		if (cutwhendtimproved && (mincost < externalbestcost))
    		{    		    			
    			break;
    		}
    		
    	}
    	else    
    		tid++;
    }

#ifdef _PRINTNAIVESTATS_
    if (odtcnt%100==0)
    {
    double sum = dtgen+dtgen2+lcagen+depthgen+lcamapgen+costgen+treespacegen+treespaceinit+btreprgen;    

	#ifdef DTCACHE
    	cout << "DTCACHE"
	#else
    	cout << "LCACOMP"
	#endif    
    		<< std::fixed << std::setw(11) << std::setprecision(4) 
    		<< odtcnt++ 
    		<< " DT=" << dtgen/sum
    		<< " DT2=" << dtgen2/sum
    		<< " LCA=" << lcagen/sum
    		<< " DEP=" << depthgen/sum
    		<< " LCAMAP=" << lcamapgen/sum
    		<< " COST="<< costgen/sum 
    		// << " TSP="<< treespacegen/sum 
    		// << " TSPINIT="<< treespaceinit/sum 
    		<< " BTREPR="<< btreprgen/sum 
    		<< " total=" << sum/1000000 << "sec" 
    		<< " lbcnt=" << lbcnt 
    		<< " lbratio=" << stoplbratio
    		<< " CST=";
    // for (int gt=0; gt<genetrees.size(); gt++)
    // 	cout << gtcost[gt] << "/" << lb[gt] << " ";
   	cout << endl;
   }
#endif   

    mincost = 0;
	for (int gt=0; gt<genetrees.size(); gt++)
		mincost += gtcost[gt];

#ifndef DTCACHE
	// free lcamaps
	for (int gt=0; gt<genetrees.size(); gt++)
		delete[] lcamaps[gt];
#endif

#ifdef _VISNETS_

	//cout << "=================================" << endl;
	cout << "Naive" <<  ++gcnt << " [";
	for (int gt=0; gt<genetrees.size(); gt++) cout << genetrees[gt]->getid() << " ";
	cout << "] ";
	cout << "] "  << endl;
	cout << *this << " ::: ";
	printrepr();

	// cout << "=================================" << endl;

	// ofstream myfile;
	// std::stringstream fname;
	// fname << "tg" << gcnt << ".dot";
 //    myfile.open (fname.str());
 //    myfile << "digraph {" << endl;
 //    myfile << "l [label=\"";
 //    printrepr(myfile);
 //    myfile << "\"]" << endl;
 //    printdot(myfile,1);
 //    myfile << "}" << endl;
 //    myfile.close();

 //    cout << fname.str() << endl;

 //	cout << endl;
#endif

	odtstats.stopnaive();

	if (verbosealg>=6)
	{
		cout << " --odtcostnaive cost=" << mincost 
			 << " dt/gt_cost_pairs=" << dtgtcnt << endl;
	}
    return mincost;        
}


int Network::_getsrcconnectors(int networkclass, NODEID *esrc)
{
	// int len=nn+rt-1; 
	// NODEID esrc[len];        // source 	
	int len=0;

	// Src edges (from node i to parent)
	for (NODEID i = 0; i<rtstartid; i++) 		
		esrc[len++]=i;

	// Add also reticulation edges 
	if (networkclass==NET_GENERAL)
		for (NODEID i = rtstartid; i < nn; i++) 
		{	
			esrc[len++]=i;  // i to parent
			esrc[len++]=-i; // i to retparent
		}

	return len; 

}

#define getvenc(nod,par) if (nod==root) { par=MAXNODEID; } else  { if (nod<0) { nod=-nod; par = retparent[nod]; } else { par = parent[nod]; } }


void Network::_getdestconnectorsforedge(int networkclass, NODEID v, int &dsrclen, NODEID* dsrc)
{

	bool reachable[size()];	
	NODEID vsrc = v;				
	NODEID p;
	getvenc(v,p);

	// src edge (v,p)
	// gen reachable 	
	getreachableto(v, reachable);


	for (NODEID i = 0; i<nn; i++)			
		if (!reachable[i]) 
		{
			NODEID w = i; // candidate
			if (networkclass==NET_TREECHILD)
			{
				// w, parent of w and sib w must be tree nodes/leaves
				if (w<rtstartid && parent[w]<rtstartid)
					if (sibling(w)<rtstartid)
					{
						dsrc[dsrclen++]=vsrc;			
						dsrc[dsrclen++]=w;						
					}
			}
			else if (networkclass==NET_CLASS1RELAXED)
			{										
				// w is a tree node or ret
				// par tree node ->  sibling not ret. 
				if (parent[w]>=rtstartid || sibling(w)<rtstartid)
				{
					
					dsrc[dsrclen++]=vsrc;					
					dsrc[dsrclen++]=w;						
				}
				

				if (w>=rtstartid) 
				{
					// check second parent
					if (retparent[w]>=rtstartid || retsibling(w)<rtstartid)
					{
						dsrc[dsrclen++]=vsrc;						
						dsrc[dsrclen++]=-w;						
					}
				}										
			}
			else 
				{					
					dsrc[dsrclen++]=vsrc;		
					dsrc[dsrclen++]=w;						
					if (w>=rtstartid) 
					{
						dsrc[dsrclen++]=vsrc;		
						dsrc[dsrclen++]=-w;	
					}
				}				
		}		
}

int Network::_getdestconnectors(int networkclass, int len, NODEID *esrc, NODEID* dsrc, bool takefirstnonempty)
{
		
	int dsrclen = 0;
	for (auto i = 0; i<len; i++)
	{		
		_getdestconnectorsforedge(networkclass, esrc[i], dsrclen, dsrc);	
		if (takefirstnonempty && dsrclen) break; // take the first non-empty set 		
	}

	return dsrclen;

}
	

Network* Network::addrandreticulation(string retid, int networkclass, bool uniform)
{
	int len = nn+rt; 
	NODEID esrc[len];        // source 

	len = _getsrcconnectors(networkclass, esrc);	

	if (!len)
	{
		cerr << "Src edge does not exist" << endl;
		return NULL; // no source edges 
	}
	
	// shuffle nodes for randomness
	shuffle(esrc, len);	


	NODEID dsrc[len*len*2]; // cand pairs	
	auto dsrclen = _getdestconnectors(networkclass, len, esrc, dsrc, !uniform);

	if (!dsrclen) return NULL;

	int pair = rand()%(dsrclen/2);
	 
	NODEID p, v = dsrc[pair*2];
	NODEID q, w = dsrc[pair*2+1];

	getvenc(v,p);		
	getvenc(w,q);


#ifdef RNDDEBUG		
		cout << " dlen=" << dlen << endl;
		cout << " v=" << v << " p=" << p << endl;

		cout << " escr=";
		for (NODEID i = 0; i < len; i++) cout << " " << esrc[i];
		cout << " dscr=";
		for (NODEID i = 0; i < dlen; i++) cout << " " << dsrc[i][0] << "," << dsrc[i][1];

		cout << endl;
		cout << v << " " << p << " -> " << w << " " << q << endl;
#endif	

	// yeah, connect (v,p) --> (w,q)
	return new Network(this, v, p, w, q, retid);

}


#if MAXRTNODES > 32
ostream& operator<<(ostream& os, const RETUSAGE& r)
{ 
	return os << r.lft <<"|" << r.rgh; 
}
#endif


void Network::initdid()
{
	if (rtcount() > 8*sizeof(DISPLAYTREEID)) 
    {
      cout << "Network has too many reticulation nodes (" << rt << "). The limit is " << 8*sizeof(DISPLAYTREEID) << "." << endl;    
      exit(-1);
    }	
}

// Only for DCE
// costfun ignored
COSTT Network::approxmindce(RootedTree &genetree, CostFun &costfun)
{
	RETUSAGE _;
	return approxmindceusage(genetree, _, costfun);
}

// Only for DCE
// costfun ignored
COSTT Network::approxmindceusage(RootedTree &genetree, RETUSAGE &retusage, CostFun &costfun)
{

#ifdef _DPDEBUG_	

	cout << " GENE TREE  -------------------------------------- " << endl;
	genetree.printdeb(cout,2);
	cout << " NETWORK ----------------------------------------- " << endl;
	printdeb(cout,2);
	cout << " ----------------------------------------- " << endl;

#endif	
	        
    DP_DCE dpdce(genetree, *this);

    dpdce.preprocess();    
    return dpdce.mindeltaroot(retusage);
}



