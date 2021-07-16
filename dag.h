#ifndef _DAG_H
#define _DAG_H


#include "tools.h"

// TODO: cycles are allowed; add top. sort.

class Dag
{

protected:
  SPID *lab;  
  SPID *retparent, *parent;
  SPID *leftchild, *rightchild, *retchild;  
  string* spid2retlabel;
  SPID lf, rt, nn, root, rtstartid;
  
  double weight;     // Dag weight (unused)
  bool exactspecies; // True iff leaf ids == species ids

  // compute the value of exactspecies
  void setexactspecies();

  // Parse from char*; allocates all needed arrays
  void parse(char *s);

  SPID _parse(char *s, int &p, int num, SPID &freeleaf, SPID &freeint, SPID &freeret, SPID* &parentset, map<string, SPID> &retlabel2spid);

  // allocates arrays based on lf and rf sizes
  void init(int _lf, int _rt);

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