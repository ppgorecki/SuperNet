#include "treespace.h"

// void TreeSpace::addedge(stcluster p, stcluster c)
// {
// 	bool fnd = 0;
// 	for (auto & cge : gtclu) 	
// 		if (EQSETS(p,cge.parent) && EQSETS(c,cge.child))
// 		{
// 			cge.usagecnt++;		
// 			fnd = true;
// 			break;
// 		}
    		
//     if (!fnd)    
//     	gtclu.push_back(CGEdge{p,c,1});
// }

// #define CACHE_STATS

#define HASHINT(a,b,c) (a>>1)^(b*c)

GNode *copytree(NODEID n, RootedTree *gt)
{
	if (n<gt->sizelf()) 		
		return new GNode({NULL,NULL,bcsingleton[gt->getlabel(n)]});
	
	GNode *a = copytree(gt->getleftchild(n),gt);
	GNode *b = copytree(gt->getrightchild(n),gt);
	return new GNode({a,b,UNION(a->cluster,b->cluster)});
}

void TreeSpace::addgenetree(RootedTree* gt)
{		 
	gtrees.push_back(copytree(gt->getroot(),gt));
	gtreesclusters.push_back(gt->getbitclusterrepr());
	gtreessrc.push_back(gt);
	if (gt->sizelf() > maxgenetreesize) maxgenetreesize = gt->sizelf();
}

// TreeSpace::TreeSpace(vector<RootedTree*> &genetrees)
// {
// 	vector<RootedTree*>::iterator gt;
// 	for (gt = genetrees.begin(); gt != genetrees.end(); ++gt)    
// 		addgenetree(*gt);

// 	for (auto & cge : gtclu) 	
// 	{
// 		ppstclusterspecies(cout, cge.parent) << " --> ";
// 		ppstclusterspecies(cout, cge.child);
// 		cout << " #" << cge.usagecnt;
// 		cout << endl;
// 	}

// }

// COSTT TreeSpace::dcerepr(RootedTree *s, bitcluster *stclu, int genetreeid)
// {

//     NODEID i;
//    	COSTT c = 0;
//    	RootedTree* gt = gtreessrc[genetreeid];
//    	bitcluster* gtclu = gtreesclusters[genetreeid];
//    	NODEID root = gt->getroot();
//    	// cout << "SPECIESTREE " << *s << endl;
 
//    	for (i=s->sizelf(); i<s->size(); i++)        
//    	{   		
//     	NODEID x = s->getleftchild(i);
//     	NODEID y = s->getrightchild(i);	
    
// 	 	for (int j=0; j<gt->size(); j++) 	
// 	 	{
// 	 		if (j==root) continue;

// 	 		NODEID p = gt->getparent(j);

// 	 		// cout << "G:";
// 	 		// ppbitclusterspecies(cout,gtclu[p]) << "->"; 
// 	 		// ppbitclusterspecies(cout,gtclu[j]) << " S:";
// 	 		// ppbitclusterspecies(cout,stclu[i]) << "->"; 
// 	 		// ppbitclusterspecies(cout,stclu[x]) << "|";
// 	 		// ppbitclusterspecies(cout,stclu[y]);
    
   
//    SPECIESTREE ((c,b),a)

// G:abc->a S:bc->c|b
// G:bc->b S:bc->c|b +Y
// G:bc->c S:bc->c|b +X
// G:abc->bc S:bc->c|b

// OK

// G:abc->a S:abc->bc|a +Y

// G:bc->b S:abc->bc|a +X ZLE!
// G:bc->c S:abc->bc|a +X ZLE

// G:abc->bc S:abc->bc|a +X

// 0.0000  6.0000  2.0000

   
// 		 	if (ISSUBSET(gtclu[j],stclu[x]))
//  			{
//  				//cout << " inx" << bool(INTERSECTION(gtclu[p],stclu[y])) <<  " " << DIFFERENCE(gtclu[p],stclu[i]) << " " ;
//  				if (INTERSECTION(gtclu[p],stclu[y]) || DIFFERENCE(gtclu[p],stclu[i])) { 
//  					c++;
//  					//cout << " +X ";
//  					//ppbitclusterspecies(cout,DIFFERENCE(gtclu[p],stclu[i]));
//  				}
// 	 		}	
// 	 		else if (ISSUBSET(gtclu[j],stclu[y]))
// 	 		{
// 	 			if (INTERSECTION(gtclu[p],stclu[x]) || DIFFERENCE(gtclu[p],stclu[i])) { 
// 	 				c++;
// 	 				//cout << " +Y ";
// 	 				//ppbitclusterspecies(cout,DIFFERENCE(gtclu[p],stclu[i]));
// 	 			}
// 	 		}
// 	 		//cout << endl;
// 		}
// 	}	
// 	return c;

// }

// #define CACHE_STATS

void TreeSpace::initleaves()
{
	COSTT *leafcost = new COSTT[gtrees.size()];

	for (int i=0; i<gtrees.size(); i++)
		leafcost[i]=0;

	for (int i=0; i<specnames.size(); i++)	
		leaves.push_back(new SNode(i, leafcost));
	
}


ostream& ppSNode(ostream& os, SNode *s)
{
   	return os;
}


// #define _PRINTTSPSTATS_

#ifdef _PRINTTSPSTATS_
clock_t lrgen =0, smpgen =0;
long sncnt = 0;
#endif

SNode *TreeSpace::find(SNode *l, SNode *r)
{

#ifdef _PRINTTSPSTATS_
	clock_t start, end;

	#define STARTCLOCK start = clock();
	#define TERMCLOCK(var) var += clock() - start;
#else
	#define STARTCLOCK
	#define TERMCLOCK(var)
#endif

#define REPRT(l) l->repr.t
#define REPRMM(l) REPRT(l)[1]
#define REPRLEN(l) REPRT(l)[0]

	if (REPRMM(l) > REPRMM(r))
	{
		//swap
		SNode *t = l;
		l = r;
		r = t;
	}

	NODEID l_lf = (REPRLEN(l)+1)/3;
	NODEID r_lf = (REPRLEN(r)+1)/3;
	NODEID lf = l_lf+r_lf;

	NODEID tr[lf*3];	
	tr[0] = lf*3-1;	// last index
	tr[1] = REPRMM(l); // mm 
	tr[2] = REPROPEN;

	memcpy ( tr+3, REPRT(l)+2, sizeof(NODEID)*(REPRLEN(l)-1));	
	memcpy ( tr+2+REPRLEN(l), REPRT(r)+2, sizeof(NODEID)*(REPRLEN(r)-1));	
	tr[tr[0]] = REPRCLOSE;

	SNode sk(tr, HASHINT(l->repr.hsh, r->repr.hsh, REPRMM(r)));

#ifdef CACHE_STATS	

	if ((n_missed+n_present)%CACHE_STATS_FREQ_REPORTING==0) 
	{
		cout << -n_missed << "+" << n_present << " " << 1.0*n_missed/(n_missed+n_present) << " nodes=" << repr2node.size();			
		cout << " Mem:" << get_memory_size() << "MB" << endl;
	}

#endif
	
	SNode *ret = NULL;

	unordered_set<SNode*, SNode::Hash, SNodeEq >::iterator it = repr2node.find(&sk);
 	if (it != repr2node.end()) 
   	{   		
		n_present++;		
#ifdef DEBUG_TSP   		

		// PPTT(l->repr.t);
		//  printf(" || ");
		//  PPTT(r->repr.t);
		//  printf(" ==> ");
		//  PPTT(res->repr.t);
		ppSNode(cout, *it); 
		cout << endl; 
#endif
		 
		ret = *it;
	}
	else
	{
	  //add new cluster 
#ifdef DEBUG_TSP 
		cout << " NEW! ";
#endif
    	n_missed++;	     	
    	ret = new SNode(gtrees.size());
    	ret->repr.t = new NODEID[lf*3];
		memcpy(ret->repr.t, tr, sizeof(NODEID)*(lf*3));    	
		ret->repr.hsh = sk.repr.hsh;
		ret->computecost2(l,r,this);	
		repr2node.insert(ret);       	

#ifdef DEBUG_TSP     	
    	cout << res->repr.hsh << " ";    	
    	cout << endl;
#endif 

   	}


#ifdef _PRINTTSPSTATS_
    	sncnt++;
    if (sncnt%1000000==0)
    {
    double sum = lrgen+smpgen;
	
    	cout << "TSP"	
    		<< std::fixed << std::setw(11) << std::setprecision(4) 
    		<< sncnt++ 
    		<< " LRG=" << lrgen/sum
    		<< " SMP=" << smpgen/sum;    		    
   		cout << endl;
   }
#endif   


// #define DEBUG_TSPP

#ifdef DEBUG_TSPP
	cout << endl << "Search: " << sk.repr.hsh << endl;
	ppSNode(cout, l); printf(" || "); 
	ppSNode(cout, r); printf(" "); 		
	cout << "size=" << repr2node.size() << endl;

 	//for (auto elt : repr2node)
	// {
	// 	cout << ":: " << elt->repr.hsh << " ";
	// 	ppSNode(cout, elt);
	// 	cout << " cst=" << elt->cost[0];
 //    	cout << endl;  	
	// }	
#endif

	return ret;
}


void SNode::computecost2(SNode *snodeleft, SNode *snoderight, TreeSpace *treespace)
{

	int genetreecnt = treespace->gtrees.size();	
	// init cost 
	// cost = new COSTT[genetreecnt];
	memcpy(cost, snodeleft->cost, sizeof(COSTT)*genetreecnt);

	bitcluster leftcluster = EMPTYSET();
	bitcluster rightcluster = EMPTYSET();

	// infer clusters

    NODEID i=4;

	if (repr.t[3] != REPROPEN)
	{
		// left tree is a leaf
		leftcluster = bcsingleton[repr.t[3]];		
	}
	else
	{		
		int open = 1;			    
	    while (1)
	    {
	        NODEID c = repr.t[i];                
	        if (c==REPROPEN) open++;                  
	        else if (c==REPRCLOSE) 
    		{
    			if (!--open) break;
    		}        	
	        else SETINSERT(leftcluster, c);                    	       	   
	        i++;
	    }
	}

	// add right leaves
    while (i<repr.t[0])
    {
    	if (repr.t[i]<REPRCLOSE) 
    		SETINSERT(rightcluster, repr.t[i]);                    
    	i++;
    }

    // print elements in cluster
	// ppSNode(cout, this) << ":::: "; 
	// ppbitclusterspecies(cout, leftcluster) << "||";
	// ppbitclusterspecies(cout, rightcluster);
	// cout << endl;

    bitcluster cluster = UNION(rightcluster,leftcluster);

    // compute dc cost for each gene tree
	GNode* stack[treespace->maxgenetreesize];
	for (int i=0; i<genetreecnt; i++)	
	 {

	 	int last=0;
	 	GNode *g = treespace->gtrees[i];
	 	cost[i]+=snoderight->cost[i];
	 	
	 	
		if (!INTERSECTION(g->cluster,cluster)) continue; // disjoint clusters; 	 	
		if (ISSUBSET(g->cluster, leftcluster)) continue;  // g maps to a child or below
		if (ISSUBSET(g->cluster, rightcluster)) continue; // g maps to a child or below

	 	stack[last++]=g; // g maps to s or above

	 	COSTT c = 0;

	 	while (last)
	 	{
	 		g = stack[--last]; // g maps to s or above
	 				
			if (ISSUBSET(g->left->cluster,leftcluster)) c++; 
			else if (ISSUBSET(g->left->cluster,rightcluster)) c++; 
			else if (INTERSECTION(g->left->cluster,cluster)) stack[last++] = g->left; // left maps to s or above

			if (ISSUBSET(g->right->cluster,leftcluster)) c++; 
			else if (ISSUBSET(g->right->cluster,rightcluster)) c++; 
			else if (INTERSECTION(g->right->cluster,cluster)) stack[last++] = g->right; // right maps to s or above
	 	}
	 	cost[i]+=c;	 	
	 }
}

