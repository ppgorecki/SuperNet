#ifndef _DAG_H
#define _DAG_H

#include <cstring>
#include "tools.h"
#include "clusters.h"


#define NET_GENERAL 2   // no limits
#define NET_CLASS1RELAXED 1    // int node has >=1 tree node/leaf child
#define NET_TREECHILD 0  // int node has <=1 ret. child

#define TIMECONSISTENT 1
#define NOTIMECONSISTENT 2

// TODO: cycles are allowed; add top. sort.

class Dag // New name: RootedDirectedGraph?
{

protected:

  // Nodes of a dag are represented using NODEID numeric type (usually short int)  
  NODEID *lab;     // Leaf labels (id of species from specnames)
  NODEID *parent;  // Parent array
  NODEID *leftchild, *rightchild; // Left and right child arrays
  NODEID *retchild, *retparent;   
  // In implementation retchild == leftchild 
  // Reticulation child and parent (only for ret. nodes)

  string* spid2retlabel;        // Dict. of reticulation labels (spid -> string)
  NODEID lf, rt, nn;  // The number of leaves, reticulations and all nodes
  NODEID root, rtstartid; // root id, start id of reticulations
  
  double weight;     // Dag weight (unused)
  bool exactspecies; // True iff leaf ids == species ids
  long int usagecount;    // for counting networks

  bool shallow; // whether the pointer arrays are allocated elsewhere (do not clean if true)

/* Example of tree structures:

> supnet -g '(a,((b,a),c))' -eD 

 Nodes(nn)=7 Leaves(lf)=4
 Leaves:0..3 TreeNodes:4..6
 root=6 exactspecies=0
0  Leaf p=6 $0 a subtree=a
1  Leaf p=4 $1 b subtree=b
2  Leaf p=4 $0 a subtree=a
3  Leaf p=5 $2 c subtree=c
4  Tree c=1 c=2 p=5 subtree=(b,a)
5  Tree c=4 c=3 p=6 subtree=((b,a),c)
6  Root c=0 c=5 subtree=(a,((b,a),c))
 parent=  0:6 1:4 2:4 3:5 4:5 5:6 6:32760
 leftchild=  4:1 5:4 6:0
 rightchild=  4:2 5:3 6:5
 depth array is uninitialized

Check also visualization:

supnet -g '(a,((b,a),c))' -eD | dot -Tpdf > n.pdf


---------------------------------------------------------
Example of network structures:

> supnet -n '((((c)#B,b))#A,(#A,(#B,a)))' -eD

 Nodes(nn)=9 Leaves(lf)=3 Reticulations(rt)=2
 Leaves:0..2 TreeNodes:3..6 RetNodes:7..8
 rtstartid=7 root=6 exactspecies=1 displaytreemaxid=4
0  Leaf p=7 $0 c subtree=c
1  Leaf p=3 $1 b subtree=b
2  Leaf p=4 $2 a subtree=a
3  Tree c=7 c=1 p=8 subtree=((c)#B,b)
4  Tree c=7 c=2 p=5 subtree=(#B,a)
5  Tree c=8 c=4 p=6 subtree=(#A,(#B,a))
6  Root c=8 c=5 subtree=((((c)#B,b))#A,(#A,(#B,a)))
7  Reti c=0 p=3 p=4 retlabel=#B subtree=(c)#B
8  Reti c=3 p=6 p=5 retlabel=#A subtree=(((c)#B,b))#A
 parent=  0:7 1:3 2:4 3:8 4:5 5:6 6:32760 7:3 8:6
 leftchild=  3:7 4:7 5:8 6:8 
 rightchild=  3:1 4:2 5:4 6:5
 lab=  0:0 1:1 2:2
 retchild=  7:0 8:3
 retparent=  7:4 8:5
 spid2retlabel=  7:#B 8:#A

Check also visualization:

supnet -n '((((c)#B,b))#A,(#A,(#B,a)))' -d | dot -Tpdf > n.pdf


*/

  // compute the value of exactspecies
  void setexactspecies();

  // Parse from char*; allocates all needed arrays
  void parse(char *s);

  NODEID _parse(char *s, int &p, int num, NODEID &freeleaf, NODEID &freeint, NODEID &freeret, NODEID* &parentset, map<string, NODEID> &retlabel2spid);

  // allocates arrays based on lf and rf sizes
  void init(int _lf, int _rt);

  void _dagrtreplace(NODEID s, NODEID d);

  virtual ostream& printdebstats(ostream&s);
  virtual ostream& printdebarrays(ostream&s);

  // iso dag tests
  inline bool _addisomap(NODEID src, NODEID dest, Dag *d, NODEID *isomap, NODEIDPair* cands, int& first, int& last);
  bool _eqdagsbypermutations(Dag *d, NODEID *isomap); 
  bool validisomap(Dag *d, NODEID *isomap);
  bool _eqdags(Dag *d, NODEID *isomap);  
  void _getclusters(NODEID v, GTCluster **res, Clusters *clusters);

  bool _marksubtrees(NODEID v, short *subtree);

public:

  // Parse Dag from a string
  Dag(const char *s, double dagweight=1.0);
  Dag(string s, double dagweight=1.0) : Dag(s.c_str(),dagweight) {}

  // Build caterpillar tree; first two leaves make cherry
  Dag(int _lf, NODEID *labels, double dagweight=1.0);

  // Make a dag by inserting new edge 
  // (v,p) --> (w,q); p,q are (ret)parents of v,w resp.
  // v may be root
  Dag(Dag *d, NODEID v, NODEID p, NODEID w, NODEID q, string retid, double dagweight=1.0);

  // Copy constructor
  Dag(const Dag &d) : Dag(d, false) {}

  // Copy constructor (shallow)
  Dag(const Dag &d, bool shallowcopy);

  virtual ~Dag();
  
  // Find leaf with a given label
  // returns MAXSP if such a leaf does not exist
  NODEID findlab(NODEID slab, int stoponerr=1);

  // Returns true if the leaf labelling is bijective
  bool bijectiveleaflabelling();


  // Sort reticulation nodes such that the first is the lowest
  void sortrtnodes();

  int _height(NODEID n, int heightarr[]);

  virtual NODEID rtcount() { return rt; }

  // Returns the parents; to get all parents use:
  // NODEID p=MAXNODEID;
  // while (getparentiter(i,p)) { .. p is the parent ... }
  // Double edges are colapsed
  bool getparentiter(NODEID i, NODEID &rparent);

  // Returns nodes; to get all nodes use:
  // NODEID i=MAXNODEID;
  // while (getnodeiter(i)) { .. i is a node ... }  
  virtual bool getnodeiter(NODEID &i);
  
  // As above but for children
  // Double edges are colapsed
  bool getchild(NODEID i, NODEID &rchild);

  // Print node and all edges to its parents in dot format
  ostream& printdot(ostream&s, int dagnum=0);

  // Print subtree rooted at i
  ostream& printsubtree(ostream&s, NODEID i, NODEID iparent = MAXSPECIES, int level = 0);

  // Print debug info 
  virtual ostream& printdeb(ostream&s, int gse, string tn="");

  // The number of nodes
  NODEID size() { return nn; }  

  // The number of nodes
  NODEID sizelf() { return lf; }  

  // Sibling of u; where parent(u) is a tree node u
  // MAXSP otherwise
  NODEID sibling(NODEID u) 
  { 
    NODEID p = parent[u];
    if (lf<=p<rtstartid) return leftchild[p]==u?rightchild[p]:leftchild[p];
    return MAXSPECIES;
  }

  // Sibling of retnode u where w.r.t tree node retparent(u)
  // MAXSP otherwise
  NODEID retsibling(NODEID u) { 
    if (u>=rtstartid)
  {
      NODEID p = retparent[u];
      if (lf<=p<rtstartid) return leftchild[p]==u?rightchild[p]:leftchild[p];
    }
    return MAXSPECIES;
  }

  friend ostream& operator<<(ostream&s, Dag &p)  { return p.print(s); }


  // Print 
  virtual ostream& print(ostream&s) 
  {
    printsubtree(s, root);
    return s;
  }

  // Check correctness of the child-parent relation in arrays
  // Return 0 if the dag is OK
  // 
  int verifychildparent();

  // -----------------------------------------------------
  // TODO: Topological sort; optionally with returning top. ordering
  bool isdag() { return true; } 

  bool istreechild() { return belongtoclass(NET_TREECHILD); }
  bool isrelaxed()   { return belongtoclass(NET_CLASS1RELAXED); }
  bool belongtoclass(int netclass);

  bool istimeconsistent();

  bool eqdags(Dag *d, bool maplabels=true);
  bool eqdagstc1(Dag *d); // only tree child and class 1 testing (to be proved)


  int compressedreprsize() { return COMPRSIZE(lf, rtcount()); } // +2 (prefix lf, rt )
  SPID *compressedrepr(SPID *r=NULL);
  ostream& printrepr(ostream&s = cout) 
  {
    SPID t[compressedreprsize()];    
    return pprepr(compressedrepr(t),s);          
  }

  bool eqdagsbypermutations(Dag *d); // brute force by permutations

  friend class DagSet;
  friend class DP_DCE;

  long int getusagecount() { return usagecount; }

  virtual NODEID getparent(NODEID v) { return parent[v]; }
  virtual NODEID getleftchild(NODEID v) { return leftchild[v]; }
  virtual NODEID getrightchild(NODEID v) { return rightchild[v]; }
  virtual NODEID getretchild(NODEID v) { return retchild[v]; }
  virtual NODEID getroot() { return root; }
  virtual NODEID getlabel(NODEID v) { return lab[v]; }  
  
  // Insert node clusters into clusters, limit to subtrees if subtreesonly=true
  void getclusters(Clusters *clusters, bool subtreesonly=false);
  
  // Check if the network has all clusters from given collection
  bool hasclusters(Clusters *clusters);

  // Returns true if arg clusters are present in network as clusters of subtrees
  bool hastreeclusters(Clusters *clusters);

  void printclusters(ostream &os);
  
};


#endif