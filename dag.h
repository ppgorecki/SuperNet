#ifndef _DAG_H
#define _DAG_H


#include "tools.h"

// TODO: cycles are allowed; add top. sort.

class Dag
{

protected:

  // Nodes of a dag are represented using SPID numeric type (usually short int)  
  SPID *lab;     // Leaf labels (id of species from specnames)
  SPID *parent;  // Parent array
  SPID *leftchild, *rightchild; // Left and right child arrays
  SPID *retchild, *retparent;   // Reticulation child and parent (only for ret. nodes)
  string* spid2retlabel;        // Dict. of reticulation labels (spid -> string)
  SPID lf, rt, nn;  // The number of leaves, reticulations and all nodes
  SPID root, rtstartid; // root id, start id of reticulations
  
  double weight;     // Dag weight (unused)
  bool exactspecies; // True iff leaf ids == species ids

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
 rightchild=  4:1 5:4 6:0
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
 rightchild=  3:7 4:7 5:8 6:8
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

  SPID _parse(char *s, int &p, int num, SPID &freeleaf, SPID &freeint, SPID &freeret, SPID* &parentset, map<string, SPID> &retlabel2spid);

  // allocates arrays based on lf and rf sizes
  void init(int _lf, int _rt);

  virtual ostream& printdebstats(ostream&s);
  virtual ostream& printdebarrays(ostream&s);

public:

  Dag(char *s, double dagweight=1.0);

  // Build caterpillar tree; first two leaves makes cherry
  Dag(int _lf, SPID *labels, double dagweight=1.0);

  // Find leaf with a given label
  // returns MAXSP if such a leaf does not exist
  SPID findlab(SPID slab, int stoponerr=1);

  // Returns true if the leaf labelling is bijective
  bool bijectiveleaflabelling();

  // Returns the parents; to get all parents use:
  // SPID p=MAXSP;
  // while (getparent(i,p)) { .. p is the parent ... }
  bool getparent(SPID i, SPID &rparent);
  
  // As above but for children
  bool getchild(SPID i, SPID &rchild);

  // Print node and all edges to its parents in dot format
  ostream& printdot(ostream&s);

  // Print subtree rooted at i
  ostream& printsubtree(ostream&s, SPID i, SPID iparent = MAXSP);

  // Print debug info 
  virtual ostream& printdeb(ostream&s, int gse, string tn="");

  // The number of nodes
  SPID size() { return nn; }

  // Return the root
  SPID getroot() { return root; }


  friend ostream& operator<<(ostream&s, Dag &p)  { return p.print(s); }


  // Print 
  virtual ostream& print(ostream&s) 
  {
    printsubtree(s, root);
    return s;
  }
    
};


#endif