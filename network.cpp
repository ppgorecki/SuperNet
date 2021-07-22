

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
	return Dag::printdebstats(s) << " displaytreemaxid=" << displaytreemaxid;
}


// Generate id'th display tree
// id encodes reticulation schema -> for long size is 8 -> 64 bits
RootedTree* Network::gendisplaytree(DISPLAYTREEID id, RootedTree *t)
{
	if (id>=displaytreemaxid) return NULL;
	if (!t)	
		t = new RootedTree(lf, lab);				
	
	SPID freeint = lf;
	//cout << "==============" << id << endl; 
	t->root = _gendisplaytree(id, t, root, MAXSP, freeint); 

	if (t->root!=t->size()-1)
	{
		cerr << "Root Error?" << t->root << endl;
		return NULL;
	}	

	t->parent[t->root] = MAXSP;
	t->depthinitialized = false;
	
	return t; 
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
	bool skipedge = (i>=rtstartid) && ((iparent == retparent[i]) == bool(id & bitmask[i-rtstartid]));

	if (skipedge) return MAXSP; // ignore

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

double Network::odtnaivecost(vector<RootedTree*> &genetrees, int costfunc, DISPLAYTREEID &optid)
{
    RootedTree *t = NULL;
    DISPLAYTREEID tid = 0; // id of display tree
    double mincost;
    SPID *lcamaps[genetrees.size()];
    while ((t=gendisplaytree(tid,t))!=NULL)       
    {       
    	double ccost = 0;    
    	t->initlca();
    	t->initdepth();
    	

    	// Compute cost of the current tree vs all gene trees
    	for (int gt=0; gt<genetrees.size(); gt++)
    	{    		
    		RootedTree *genetree = genetrees[gt];
    		SPID *lcamap = lcamaps[gt];
    		if (tid) genetree->getlcamapping(*t,lcamaps[gt]);    // overwrite previous
    		else lcamaps[gt] = genetrees[gt]->getlcamapping(*t); // init lcamap
    		ccost += genetree->_cost(*t, lcamaps[gt], costfunc);
    	}
        
        if (!tid || mincost>ccost) 
        { 
        	// Store min
            mincost = ccost;
            optid = tid;
        }
        tid++;
    }
    return mincost;        
}