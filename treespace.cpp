#include "treespace.h"
#include "costs.h"
#include <algorithm>

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
	if (used())
	{
		cerr << "Tree space already initialized. Set all gene trees before cost computations" << endl;
		exit(-1);
	}
	gtrees.push_back(copytree(gt->getroot(),gt));
	bitcluster *bc = gt->getbitclusterrepr();
	gtreesclusters.push_back(bc);

	// Snapshot internal-node clusters into a sorted vector for fast
	// "cluster ∈ gt clusters" lookups (DTCACHE RF).
	vector<bitcluster> sorted;
	NODEID glf = gt->sizelf();
	NODEID gnn = gt->size();
	sorted.reserve(gnn - glf);
	for (NODEID i = glf; i < gnn; i++) sorted.push_back(bc[i]);
	std::sort(sorted.begin(), sorted.end());
	gtreesclustersorted.push_back(std::move(sorted));

	gtreessrc.push_back(gt);
	if (gt->sizelf() > maxgenetreesize) maxgenetreesize = gt->sizelf();
}


void TreeSpace::initleaves()
{
	COSTT *leafcost    = new COSTT[gtrees.size()];
	COSTT *leafcost_rf = new COSTT[gtrees.size()];

	for (size_t i = 0; i < gtrees.size(); i++)
	{
		leafcost[i]    = 0;
		leafcost_rf[i] = 0;   // a singleton leaf is not an internal cluster
	}

	for (size_t i = 0; i < specnames.size(); i++)
		leaves.push_back(new SNode(i, leafcost, leafcost_rf));
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

	// Canonicalise: smaller min-leaf side on the left.
	if (l->minleaf > r->minleaf)
	{
		SNode *t = l;
		l = r;
		r = t;
	}

	// Lookup key — pointer-compare on (left, right). No tr[] reconstruction.
	SNode sk(l, r, l->minleaf,
	         HASHINT(l->repr.hsh, r->repr.hsh, (size_t)r->minleaf));

#ifdef CACHE_STATS

	if ((n_missed+n_present)%CACHE_STATS_FREQ_REPORTING==0)
	{
		cout << -n_missed << "+" << n_present << " " << 1.0*n_missed/(n_missed+n_present) << " nodes=" << repr2node.size();
		cout << " Mem:" << get_memory_size() << "MB" << endl;
	}

#endif

	SNode *ret = repr2node.find(&sk);
 	if (ret)
   	{
		n_present++;
	}
	else
	{
    	n_missed++;
    	ret = new SNode(gtrees.size(), l, r, l->minleaf);
		ret->repr.hsh = sk.repr.hsh;
		ret->computecost2(l, r, this);
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
	int act = treespace->active_costtype;
	bool need_dc = (act == 0 || act == COSTDEEPCOAL || act == COSTDEEPCOALEDGE);
	bool need_rf = (act == 0 || act == COSTROBINSONFOULDS);

	// init cost (DC) — only when needed
	if (need_dc)
		memcpy(cost, snodeleft->cost, sizeof(COSTT)*genetreecnt);

	// Clusters come straight from children — set in their constructors.
	// `this->cluster` is already set (UNION) by SNode's constructor.
	bitcluster leftcluster  = snodeleft->cluster;
	bitcluster rightcluster = snoderight->cluster;
	bitcluster cluster      = this->cluster;

    // Loop over gene trees, computing whichever halves are active.
    // DC: walk the gene-tree subtree and count edges mapped below s.
    // RF: O(log |G|) binary_search whether the cluster matches an
    // internal cluster in the gene tree.

	GNode* stack[treespace->maxgenetreesize];
	for (int i=0; i<genetreecnt; i++)
	 {
		// RF: inherit + bonus.
		if (need_rf)
		{
			cost_rf[i] = snodeleft->cost_rf[i] + snoderight->cost_rf[i];
			const std::vector<bitcluster> &gtclusters = treespace->gtreesclustersorted[i];
			if (std::binary_search(gtclusters.begin(), gtclusters.end(), cluster))
				cost_rf[i] += 1;
		}

		if (!need_dc) continue;

	 	int last=0;
	 	GNode *g = treespace->gtrees[i];
	 	cost[i]+=snoderight->cost[i];

		if (EMPTYINTERSECTION(g->cluster, cluster)) continue; // disjoint clusters;
		if (ISSUBSET(g->cluster, leftcluster)) continue;  // g maps to a child or below
		if (ISSUBSET(g->cluster, rightcluster)) continue; // g maps to a child or below

	 	stack[last++]=g; // g maps to s or above

	 	COSTT c = 0;

	 	while (last)
	 	{
	 		g = stack[--last]; // g maps to s or above

			if (ISSUBSET(g->left->cluster, leftcluster)) c++;
			else if (ISSUBSET(g->left->cluster,rightcluster)) c++;
			else if (NONEMPTYINTERSECTION(g->left->cluster, cluster)) stack[last++] = g->left; // left maps to s or above

			if (ISSUBSET(g->right->cluster,leftcluster)) c++;
			else if (ISSUBSET(g->right->cluster,rightcluster)) c++;
			else if (NONEMPTYINTERSECTION(g->right->cluster,cluster)) stack[last++] = g->right; // right maps to s or above
	 	}
	 	cost[i]+=c;
	 }
}

