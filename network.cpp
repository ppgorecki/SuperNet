

#include "rtree.h"
#include "network.h"

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
	
	SPID freeint = lf;
	// cerr << "==============" << id << " " << rtcount() << endl; 
	t->root = _gendisplaytree(id, t, getroot(), MAXSP, freeint); 
	t->parent[t->root] = MAXSP;
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


bool Network::_skiprtedge(SPID i, SPID iparent, DISPLAYTREEID id)
{
	return ((iparent == retparent[i]) == bool(id & bitmask[i-rtstartid]));
}

SPID Network::_gendisplaytree(DISPLAYTREEID id, RootedTree *t, SPID i, SPID iparent, SPID &freeint) 
{ 
	// leaf
	if (i < lf) return i; 

	// if (i>=rtstartid)
	// {
	// 	cout << bitmask[i-rtstartid] << " __ " << (id & bitmask[i-rtstartid]) << endl;
	// }

	// skip current edge
	if (i>=rtstartid)
	{		
		//cout << "RT check!" << i << " " << rtstartid << endl;
		if (_skiprtedge(i, iparent, id)) return MAXSP;		
	}
		

	SPID cleft = MAXSP;
	getchild(i,cleft);
	SPID cleftgen = _gendisplaytree(id, t, cleft, i, freeint);

	// cout << (int)i << " cleft" << (int)cleft << " GD-LEFT " << (int)cleftgen << endl;

	SPID cright = cleft;
	if (!getchild(i,cright))
	{
		// reticulation only
		return cleftgen;  // maybe MAXSP (ignore)
	}

	SPID crightgen = _gendisplaytree(id, t, cright, i, freeint);

	if (cleftgen==MAXSP) return crightgen;  // maybe MAXSP
	if (crightgen==MAXSP) return cleftgen;  // single node

	// two nodes to be connected

	// cerr << freeint << " " << cleftgen << " " << crightgen << endl;
	t->parent[cleftgen] = t->parent[crightgen] = freeint;
	t->leftchild[freeint] = cleftgen;
	t->rightchild[freeint] = crightgen;

	return freeint++;
}

// Marks nodes reachable from v (including v)	
void Network::getreachablefrom(SPID v, bool *reachable)
{
	bool visited[nn];
	for (SPID i=0; i<nn; i++) reachable[i] = visited[i] = false;
	_getreachablefrom(v,reachable, visited);
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
	for (SPID i=0; i<nn; i++)		
	{	
		getreachablefrom(i, vreachable);
		int cnt = 0;
		for (SPID j=0; j<last; j++)		
			if (vreachable[j]) cnt++;
		s << cnt << " ";
	}
	return s;

}

void Network::_getreachablefrom(SPID v, bool *reachable, bool *visited)
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

// Marks nodes w such that v is reachble from w (including v)	
void Network::getreachableto(SPID v, bool *reachable)
{
	bool visited[nn];
	for (SPID i=0; i<nn; i++) reachable[i] = visited[i] = false;
	_getreachableto(v,reachable, visited);
}

void Network::_getreachableto(SPID v, bool *reachable, bool *visited)
{
	if (visited[v]) return;
	reachable[v] = true;
	if (v==root) return;
	_getreachableto(parent[v], reachable, visited);
	if (v>=rtstartid)  
		_getreachableto(retparent[v],reachable, visited);		
}

double Network::odtcost(vector<RootedTree*> &genetrees, int costfunc, bool usenaive, int runnaiveleqrt)
{
	if (usenaive) return odtcostnaive(genetrees, costfunc);
	return odtcostdpbb(genetrees, costfunc, runnaiveleqrt);
}

double Network::odtcostdpbb(vector<RootedTree*> &genetrees, int costfunc, int runnaiveleqrt)
{

	if (costfunc!=COSTDEEPCOAL)
	{
		cerr << "DP&BB cost computation only for DC" << endl;
		exit(-1);
	}

	double cost = 0;
	for (int gt=0; gt<genetrees.size(); gt++)
		cost+=mindc(*genetrees[gt], runnaiveleqrt);
    
    return cost;        
}

double Network::odtcostnaive(RootedTree *genetree, int costfunc)
{
	vector<RootedTree*> genetrees;
	genetrees.push_back(genetree);
	return odtcostnaive(genetrees,costfunc);
}

extern bool print_repr_inodtnaive;

double Network::odtcostnaive(vector<RootedTree*> &genetrees, int costfunc)
{
    RootedTree *t = NULL;
    DISPLAYTREEID tid = 0; // id of display tree
    double mincost;
    SPID *lcamaps[genetrees.size()];
    double gtcost[genetrees.size()];
        
    while ((t=gendisplaytree(tid,t))!=NULL)       
    {               	

    	t->initlca();
    	t->initdepth();    	 	

    	// Compute cost of the current tree vs all gene trees
    	for (int gt=0; gt<genetrees.size(); gt++)
    	{    		
    		RootedTree *genetree = genetrees[gt];
    		SPID *lcamap = lcamaps[gt];
    		if (tid) genetree->getlcamapping(*t,lcamaps[gt]);    // overwrite previous    			    		
    		else lcamaps[gt] = genetrees[gt]->getlcamapping(*t); // init lcamap

    		double curcost = genetree->_cost(*t, lcamaps[gt], costfunc);

    		if (!tid || (gtcost[gt] > curcost)) gtcost[gt] = curcost;  

    		if (print_repr_inodtnaive)
    		{
    			cout << "!";
    			t->printrepr(cout) << " " << curcost << endl;
    		}
    	}

    	tid++;

    }

    mincost = 0;
	for (int gt=0; gt<genetrees.size(); gt++)
		mincost += gtcost[gt];

    return mincost;        
}

// #define RNDDEBUG 

/*
	Add random reticulation	

*/
Network* Network::addrandreticulation(string retid, int networktype, bool uniform)
{
	int len=nn+rt-1; 
	SPID esrc[len];        // source 
	SPID dsrc[len*len][2]; // cand pairs
	
	len=0;
	for (SPID i = 0; i<rtstartid; i++) 		
		esrc[len++]=i;

	if (networktype==NT_GENERAL)
		for (SPID i = rtstartid; i < nn; i++) 
		{	
			esrc[len++]=i;
			esrc[len++]=-i;
		}

	if (!len)
	{
		cerr << "Src edge does not exist" << endl;
		return NULL; // no source edges 
	}

	
	// shuffle nodes for randomness
	shuffle(esrc,len);	
	bool reachable[size()];


// set v and parent
#define getvenc(nod,par) if (nod==root) { par=MAXSP; } else  { if (nod<0) { nod=-nod; par = retparent[nod]; } else { par = parent[nod]; } }

	for (SPID i = 0; i<len; i++)
	{		
		SPID v = esrc[i];
		SPID vsrc = v;
		SPID p = MAXSP;
		SPID dlen = 0;
		getvenc(v,p);
		
		// src edge (v,p)
		// gen reachable 	
		getreachableto(v, reachable);

		// cout << "CAND v=" << v << " ::";

		for (SPID i = 0; i<nn; i++)			
			if (!reachable[i]) 
			{
				SPID w = i; // candidate
				if (networktype==NT_TREECHLD)
				{
					// w, parent of w and sib w must be tree nodes/leaves
					if (w<rtstartid && parent[w]<rtstartid)
						if (sibling(w)<rtstartid)
						{
							dsrc[dlen][0]=vsrc;			
							dsrc[dlen++][1]=w;						
						}
				}
				else if (networktype==NT_CLASS1)
				{										
					// w is a tree node or ret
					// par tree node ->  sibling not ret. 
					if (parent[w]>=rtstartid || sibling(w)<rtstartid)
					{
						
						dsrc[dlen][0]=vsrc;					
						dsrc[dlen++][1]=w;						
					}
					

					if (w>=rtstartid) 
					{
						// check second parent
						if (retparent[w]>=rtstartid || retsibling(w)<rtstartid)
						{
							dsrc[dlen][0]=vsrc;						
							dsrc[dlen++][1]=-w;						
						}
					}										
				}
				else 
					{					
						dsrc[dlen][0]=vsrc;		
						dsrc[dlen++][1]=w;						
						if (w>=rtstartid) 
						{
							dsrc[dlen][0]=vsrc;		
							dsrc[dlen++][1]=-w;	
						}
					}


			}		
		
		if (!dlen) continue;

		v = dsrc[rand()%dlen][0];
		SPID w = dsrc[rand()%dlen][1];
		SPID q;
		getvenc(w,q);
	

		// yeah, connect (v,p) --> (w,q)

#ifdef RNDDEBUG		
		cout << " v=" << v << " p=" << p << endl;
		cout << " escr=";
		for (SPID i = 0; i < len; i++) cout << " " << esrc[i] ;
		cout << " dscr=";
		for (SPID i = 0; i < dlen; i++) cout << " " << dsrc[i] ;
		cout << endl;
		cout << v << " " << p << " -> " << w << " " << q << endl;
#endif	

		
		return new Network(this, v, p, w, q, retid);

	}
	return NULL;
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
