
#ifndef _TREESPACE_
#define _TREESPACE_

#include "bitcluster.h"
#include "tools.h"
#include "rtree.h"
#include <unordered_set>

constexpr std::uint32_t log2(std::uint32_t n) noexcept
{
    return (n > 1) ? 1 + log2(n >> 1) : 0;
};


#define CACHE_STATS_FREQ_REPORTING 1000000

#define PPTT(r) cout << "["; for (int i=0;i<(r)[0];i++) cout << (NODEID)((r)[i]) << "."; cout << (NODEID)(r)[(r)[0]] << "]"

typedef struct NODEIDKey
{
 	NODEID* t;
 	size_t hsh;
   
    struct Hash 
	{
    	size_t operator() (NODEIDKey const& t) const 
    	{       
        	return t.hsh;
    	}
	};

} NODEIDKey;

// struct TreeRepr
// {
// 	NODEID lf;
//   	NODEID *spp;
//   	char *bvec;
  	
//   	// Repr:
//   	// 10 - (
//   	// 11 - )
//   	// 00 - species  

//   TreeRepr(int lfs):
//   {
//   		lf = lfs;
//   		spp = new NODEID[lf];
//   		t = new char[];
//   }

//   int void add(NODEID leaf, int bpos)
//   {
//   	int start = bpos/(8*sizeof(size_t));
//   	int off = bpos%(8*sizeof(size_t));
//   	t[bpos/(sizeof(size_t)*8)][]		
//   }

//   // bool operator==(const TreeRepr& a) const
//   // {  
//   //     if (t[0]!=a.t[0]) return false;
//   //     for (int i=1; i<=a.t[0];i++)
//   //     	if (t[i]!=a.t[i]) return false;
//   //     return true;

//   // }
// };

// namespace std {

//   template <>
//   struct hash<NODEIDKey>
//   {
//     std::size_t operator()(const NODEIDKey& k) const
//     {
//       using std::size_t;
//       using std::hash;
//       using std::string;

//       // Compute individual hash values for first,
//       // second and third and combine them using XOR
//       // and bit shifting:
//  	  std::size_t seed = 0;
//       for (int i=2;i<=k.t[0];i++) seed = (seed << 1)^k.t[i];
//       return seed;
//     }
//   };

// }

typedef struct GNode 
{
	GNode *left, *right;		
	bitcluster cluster;	
	NODEID *repr; // node repr

	~GNode()
	{
		if (left) delete left;
		if (right) delete right;
	}

} GNode;

class TreeSpace;

typedef struct SNode
{
	NODEIDKey repr;    // hash key (repr.hsh used; repr.t now lazy)
	COSTT *cost;       // DC-flavoured cost vs each gene tree (DCE-encoded)
	COSTT *cost_rf;    // RF matches count vs each gene tree (#internal nodes
	                   // of this subgraph whose cluster appears in gt's clusters)
	SNode *left;       // canonical-form left child (NULL for leaves and tmp keys)
	SNode *right;      // canonical-form right child (NULL for leaves and tmp keys)
	bitcluster cluster;// union of leaves under this subgraph
	NODEID minleaf;    // min leaf label under this subgraph (used by find canonicalisation)


	// internal node — l and r must already be canonical (from find/leaves).
	SNode(int gts, SNode *l, SNode *r, NODEID minlf)
	{
		cost = new COSTT[gts];
		cost_rf = new COSTT[gts];
		left = l;
		right = r;
		minleaf = minlf;
		cluster = UNION(l->cluster, r->cluster);
		repr.t = NULL;     // not used by lookup/computecost2 anymore
	}

	void computecost2(SNode *snodeleft, SNode *snoderight, TreeSpace *treespace);

	// tmp struct for find — pure lookup key, no cost arrays.
	SNode(SNode *l, SNode *r, NODEID minlf, size_t phsh)
	{
		repr.t = NULL;
		repr.hsh = phsh;
		cost = NULL;
		cost_rf = NULL;
		left = l;
		right = r;
		minleaf = minlf;
	}

	// leaf
	SNode(NODEID l, COSTT *cst, COSTT *cst_rf)
	{

		cost = cst;
		cost_rf = cst_rf;
		repr.t = NULL;
		repr.hsh = (size_t)rand64(); // rand hash for leaves
		left = NULL;
		right = NULL;
		minleaf = l;
		cluster = bcsingleton[l];
	}

	~SNode()
	{
		if (cost) // clean internal only
		{
			delete[] cost;
			delete[] cost_rf;
		}
	}

	struct Hash 
	{
    	size_t operator() (SNode* const& t) const 
    	{       
        	return t->repr.hsh;
    	}
	};

} SNode;



// struct SNodeHash 
// {
//     size_t operator() (SNode* const& t) const 
//     {
//         // std::size_t seed = 0;
//         // for (int i=2;i<=t->repr[0];i++) seed ^= t->repr[i];
//         return t->repr.hsh;
//     }
// };


struct SNodeEq {
    // SNodes are canonicalised by find() so each subgraph has a unique
    // representative. Two internal SNodes are equal iff they share the same
    // (canonical) child pointers. Pointer compare is O(1).
    bool operator()( SNode const *a, SNode const *b) const
    {
      return a->left == b->left && a->right == b->right;
    }
};


// Open-addressing flat hash set for SNode* lookups. Replaces the STL
// `std::unordered_set<SNode*, SNode::Hash, SNodeEq>` used for the DTCACHE
// dedup table. Three improvements:
//   - Open addressing with linear probing → no separate-chaining allocations,
//     much better cache locality.
//   - Cached hash per bucket → skip hash recompute and skip pointer compare
//     when hashes differ.
//   - Power-of-2 size with bitmask → cheaper modulo.
// Custom because the use case is narrow (we control hash + eq, single value
// type) and a single-header dependency would carry far more code.
class SNodeFlatSet
{
	struct Bucket
	{
		SNode  *node;   // NULL = empty slot
		size_t  hsh;
	};
	std::vector<Bucket> buckets;
	size_t mask;
	size_t count_;

	void resize_to(size_t newcap)
	{
		std::vector<Bucket> old(std::move(buckets));
		buckets.assign(newcap, Bucket{NULL, 0});
		mask = newcap - 1;
		count_ = 0;
		for (auto &b : old)
			if (b.node) insert_known_unique(b.node, b.hsh);
	}

	inline void insert_known_unique(SNode *node, size_t hsh)
	{
		size_t i = hsh & mask;
		while (buckets[i].node) i = (i + 1) & mask;
		buckets[i].node = node;
		buckets[i].hsh  = hsh;
		count_++;
	}

public:
	SNodeFlatSet(size_t initial_cap = 1024)
	{
		size_t n = 1; while (n < initial_cap) n <<= 1;
		buckets.assign(n, Bucket{NULL, 0});
		mask = n - 1;
		count_ = 0;
	}

	// Lookup by a temporary key. Returns NULL on miss.
	SNode *find(SNode *key) const
	{
		size_t hsh = key->repr.hsh;
		size_t i = hsh & mask;
		while (buckets[i].node)
		{
			if (buckets[i].hsh == hsh
			    && buckets[i].node->left  == key->left
			    && buckets[i].node->right == key->right)
				return buckets[i].node;
			i = (i + 1) & mask;
		}
		return NULL;
	}

	// Insert assuming the key is not present. Caller must have checked find().
	void insert(SNode *node)
	{
		// Resize at ~87.5% load to keep probe sequences short.
		if ((count_ + 1) * 8 >= buckets.size() * 7)
			resize_to(buckets.size() * 2);
		insert_known_unique(node, node->repr.hsh);
	}

	void clear()
	{
		for (auto &b : buckets) b.node = NULL;
		count_ = 0;
	}

	size_t size() const { return count_; }

	template <typename F>
	void for_each(F f) const
	{
		for (auto &b : buckets) if (b.node) f(b.node);
	}
};


// Stores all data gene tree nodes
class TreeSpace
{
	// SNode dedup table — open-addressing flat hash. ~30% faster than
	// std::unordered_set on the find-heavy hot path.
	SNodeFlatSet repr2node;

	vector<RootedTree*> gtreessrc;
	vector<bitcluster*> gtreesclusters;
	vector<SNode*> leaves;

public:
	// Per gene tree: sorted vector of bitclusters of internal nodes (incl.
	// root). Used by DTCACHE RF lookups (binary_search to test whether a
	// display-tree subgraph cluster is present in the gene tree).
	vector<vector<bitcluster>> gtreesclustersorted;
private:

	long n_missed;
	long n_present;
	int maxdisplaytreecachesize;

	void initleaves();

	public: 
		
		// All nodes from gene trees
		vector<GNode*> gtrees;

		int maxgenetreesize;
	
		TreeSpace(RootedTree *g, int _maxdisplaytreecachesize)
		 {
			n_missed = 0;
			n_present = 0;
			maxgenetreesize = 0;
			maxdisplaytreecachesize = _maxdisplaytreecachesize;
			active_costtype = 0;
			addgenetree(g);
			initleaves();
		}

		TreeSpace(vector<RootedTree*> &genetrees, int _maxdisplaytreecachesize)
		{
			n_missed = 0;
			n_present = 0;
			maxgenetreesize = 0;
			maxdisplaytreecachesize = _maxdisplaytreecachesize;
			active_costtype = 0;
			for (auto g: genetrees) addgenetree(g);
			initleaves();
		}

		// Set the cost type that will be queried; computecost2 will then
		// skip building the unused half (e.g., for DC runs we can skip the
		// RF cluster lookups). Default 0 = compute both (safe fallback).
		// Should be called before any cost-driven cache miss.
		int active_costtype;
		void setactivecost(int ct) { active_costtype = ct; }

		~TreeSpace()
		{
			clearcache();			
			for (auto v: gtrees) delete v;
		}
		
		void addgenetree(RootedTree *g);

		// Returns a parent with l r kids; l, r are assumed to be present  
		// Creates the parent if not present
		SNode *find(SNode *l, SNode *r);

		// COSTT dce(RootedTree *s);
		// COSTT dcerepr(RootedTree *s, bitcluster *bitclu, int genetreeid);
		SNode *leaf(NODEID i) { return leaves[i]; }
	
		/*
		Maps: TrRepr -> Node, fast lookup for parents
		*/
		void computenodecost(SNode *s);

		void treecompleted()
		{
			if (repr2node.size()>maxdisplaytreecachesize)			
			{
#ifdef CACHE_STATS				
			 	cout << "CLEAR" << get_memory_size() << "MB nodes=" << repr2node.size() << endl;
#endif			 	
			 	clearcache();			 	
			}
		}

		void clearcache()
		{
			// repr2node contains only internal SNodes (leaves live in
			// `leaves[]`), so we can just delete everything here.
			repr2node.for_each([](SNode *s) { delete s; });
			repr2node.clear();
		}

		bool used()
		{	
			return n_missed || n_present;
		}
};

ostream& ppSNode(ostream& os, SNode *s);

#endif
