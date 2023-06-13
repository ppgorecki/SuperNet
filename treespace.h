
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


#define PPTT(r) cout << "["; for (int i=0;i<(r)[0];i++) cout << (NODEID)((r)[i]) << "."; cout << (NODEID)(r)[(r)[0]] << "]"

typedef struct NODEIDKey
{
 	NODEID* t;
 	size_t hsh;

    // bool operator==( const NODEIDKey &a) const 
    // {
    // 	cout << "#"; PPTT(t); 
    // 	cout << "vs"; 
    // 	PPTT(a.t);
    	
    //   	if (t[0] != a.t[0]) 
    //   	{
    //   		return false;
    //   	}

	// 	for (int i=1; i<=a.t[0]; i++)
	// 	{
	// 	   if (t[i]!=a.t[i]) return false;
	// 	}

	//     return true;
    // }

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
} GNode;

class TreeSpace;

typedef struct SNode 
{
	NODEIDKey repr; // node repr		
	COSTT *cost;	

	// internal node
	SNode(int gts)
	{
		cost = new COSTT[gts];		
	}
	 		   
   /*
   SPECIESTREE ((c,b),a)

G:abc->a S:bc->c|b
G:bc->b S:bc->c|b +Y
G:bc->c S:bc->c|b +X
G:abc->bc S:bc->c|b

OK

G:abc->a S:abc->bc|a +Y

G:bc->b S:abc->bc|a +X ZLE!
G:bc->c S:abc->bc|a +X ZLE

G:abc->bc S:abc->bc|a +X

0.0000  6.0000  2.0000

   */


		
		//  	if (ISSUBSET(gtclu[j],stclu[x]))
 	// 		{
 	// 			//cout << " inx" << bool(INTERSECTION(gtclu[p],stclu[y])) <<  " " << DIFFERENCE(gtclu[p],stclu[i]) << " " ;
 	// 			if (INTERSECTION(gtclu[p],stclu[y]) || DIFFERENCE(gtclu[p],stclu[i])) { 
 	// 				c++;
 	// 				//cout << " +X ";
 	// 				//ppbitclusterspecies(cout,DIFFERENCE(gtclu[p],stclu[i]));
 	// 			}
	 // 		}	
	 // 		else if (ISSUBSET(gtclu[j],stclu[y]))
	 // 		{
	 // 			if (INTERSECTION(gtclu[p],stclu[x]) || DIFFERENCE(gtclu[p],stclu[i])) { 
	 // 				c++;
	 // 				//cout << " +Y ";
	 // 				//ppbitclusterspecies(cout,DIFFERENCE(gtclu[p],stclu[i]));
	 // 			}
	 // 		}
	 // 		//cout << endl;
		// }

	
	void computecost2(SNode *snodeleft, SNode *snoderight, TreeSpace *treespace);
		

	// tmp struct for find
	SNode(NODEID *pt, size_t phsh)
	{
		repr.t = pt;
		repr.hsh = phsh;
		cost = NULL;
	}

	// leaf
	SNode(NODEID l, COSTT *cst)
	{		
		
		cost = cst;
		repr.t = new NODEID[3];
		repr.t[0] = 2;		
		repr.t[1] = l;		
		repr.t[2] = l;
		repr.hsh = (size_t)rand64(); // rand hash for leaves
	}

	~SNode()
	{
		if (cost && left) 
		{ 
			delete[] repr.t;
			delete[] cost; // clean internal only
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

    bool operator()( SNode const *a, SNode const *b) const 
    {
    	// cout << "#" << a->repr.hsh << "?" << b->repr.hsh << endl;
      	if (a->repr.t[0]!=b->repr.t[0]) return false;
		for (int i=1; i<=a->repr.t[0]; i++)
		   if (a->repr.t[i]!=b->repr.t[i]) return false;
	    return true;
    }
};


// Stores all data gene tree nodes
class TreeSpace
{

	// needed structs:
	// SpID -> LeafNode 
	// NodeRepr -> Node
	// vector<CGEdge> gtedges;
 
	// std::unordered_map<NODEID*, SNode*, comparespids> repr2node;	
	// std::unordered_map<NODEID*, SNode*> repr2node;	
	std::unordered_set<SNode*, SNode::Hash, SNodeEq > repr2node;	
	
	vector<RootedTree*> gtreessrc;
	vector<bitcluster*> gtreesclusters;
	vector<SNode*> leaves;

	long n_missed;
	long n_present;

	void initleaves();

	public: 
		
		// All nodes from gene trees
		vector<GNode*> gtrees;

		int maxgenetreesize;
	
		TreeSpace(RootedTree *g)
		 { 
			n_missed = 0;
			n_present = 0;
			maxgenetreesize = 0;
			addgenetree(g); 
			initleaves();			
		}

		TreeSpace(vector<RootedTree*> &genetrees) 
		{
			n_missed = 0;
			n_present = 0;
			maxgenetreesize = 0;
			for (auto g: genetrees) addgenetree(g);
			initleaves();
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
			if (get_memory_size()>10) 
			{

				//cout << "CLEAR" << get_memory_size() << "MB nodes=" <<  repr2node.size() << endl;
				//clearcache();
				//cout << "After" << get_memory_size() << endl;
			}
		}

		void clearcache()
		{
			// remove all internal SNodes (of the size >3)
			// std::unordered_set<SNode*, SNode::Hash, SNodeEq > repr2node;
			for (const auto& element : repr2node) 
			{
				//delete[] element->repr.t;
    		    delete element;
    		}
    		repr2node.clear();
		}
		

};

ostream& ppSNode(ostream& os, SNode *s);

#endif
