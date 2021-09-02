

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

 

Network* Network::addrandreticulation(string retid, int networktype)
{
	int len=nn+rt-1; // v<->parent.v, addtionally v->retparent; exclude v==root
	SPID esrc[len];
	SPID dsrc[len];
	
	len=0;
	for (SPID i = 0; i<rtstartid; i++) 		
		if (i!=root) esrc[len++]=i;

	if (networktype==NT_GENERAL)
		for (SPID i = rtstartid; i < nn; i++) 
		{	
			esrc[len++]=i;
			esrc[len++]=-i;
		}

	// for (SPID i = 0; i<len; i++)
	// 	cout << esrc[i] << " ";
	// cout << endl;

	if (!len)
	{
		cerr << "Src edge does not exist" << endl;
		return NULL; // no source edges 
	}

	// shuffle nodes for randomness
	shuffle(esrc,len);
	bool reachable[size()];	
	for (SPID i = 0; i<len; i++)
	{		
		SPID v = esrc[i];
		SPID p;
		SPID dlen = 0;
		if (v<0)
		{
			v=-v;
			p = retparent[v];			
		}
		else
			p = parent[v];			
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
							dsrc[dlen++]=w;						
				}
				else if (networktype==NT_CLASS1)
				{										
					// w is a tree node or ret
					// par tree node ->  sibling not ret. 
					if (parent[w]>=rtstartid || sibling(w)<rtstartid)
						dsrc[dlen++]=w;						
					

					if (w>=rtstartid) 
					{
						// check second parent
						if (retparent[w]>=rtstartid || retsibling(w)<rtstartid)
							dsrc[dlen++]=-w;						
					}										
				}
				else 
					{
						dsrc[dlen++]=w;						
						if (w>=rtstartid) 
							dsrc[dlen++]=-w;	
					}


			}		
		

		if (!dlen) continue;

		SPID w = dsrc[rand()%dlen];
		SPID q;


		if (w<0) 
		{ 
			w = -w;
			q = retparent[w];
		}
		else q = parent[w];

		// yeah, connect (v,p) --> (w,q)

#ifdef RNDDEBUG		
		cout << " v=" << v << " p=" << p << endl;
		cout << " dscr=";
		for (SPID i = 0; i < dlen; i++) cout << " " << dsrc[i] ;
		cout << endl;
		cout << v << " " << p << " -> " << w << " " << q << endl;
#endif	

		
		return new Network(this, v, p, w, q, retid);

	}






	return NULL;
}