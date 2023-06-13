
#ifndef _CLUSTERS_H_
#define _CLUSTERS_H_

#include "tools.h"

class RootedTree;

#define deletespcluster(s) if (s!=topspcluster && ((s)[0]!=1)) delete[] s;

typedef struct GTCluster
{
  GTCluster *l, *r; // child nodes composition; needed to compute lca
  NODEID type;        // >=0 species id, -1 - all species, -2 - internal l-r
  NODEID lcamap;      // map to species tree
  NODEID* spcluster;  // cluster
  int usagecnt;     
  
  GTCluster(GTCluster *lc, GTCluster *rc, NODEID *s) : l(lc), r(rc), type(-2), lcamap(-1), spcluster(s), usagecnt(0)  
  { if (s[0] == 1) type = s[1]; }

  friend ostream& operator<<(ostream&s, GTCluster &p) {
    printspcluster(s, p.spcluster);
    s << " #"<< p.usagecnt << " T" << p.type;
    return s;
  }

  ~GTCluster() {  deletespcluster(spcluster);  }

} GTCluster;

/*
 Collection of GTC clusters
*/
class TreeClusters
{
  int ssize;
  int _usagecnt;  
  map<NODEID*, GTCluster*, comparespids> t;
  
  vector<GTCluster*> internal;
  vector<GTCluster*> leaves;

public:
  TreeClusters()
  {
    ssize = specnames.size();
    _usagecnt = 0;
    // initialize singleton clusters

    for (NODEID i = 0; i < ssize; i++)
    {
      GTCluster *gc = new GTCluster(0, 0, spec2spcluster[i]);
      t[gc->spcluster] = gc;
      leaves.push_back(gc);
    }
  }

  ~TreeClusters()
  {
    for (int i = 0; i < ssize; i++)
        delete leaves[i];

    for (size_t i = 0; i < internal.size(); i++)
        delete internal[i];      
  }

    
  GTCluster* leafcluster(NODEID n) 
  { 
    leaves[n]->usagecnt++; 
    _usagecnt++; 
    return leaves[n]; 
  }

  GTCluster* get(GTCluster *l, GTCluster *r);
  int size() { return t.size(); }
  int usagecnt() { return _usagecnt; }

  void addtree(RootedTree *g);

  string genrootedquasiconsensus(RootedTree *preserveroottree);

  friend ostream& operator<<(ostream&s, TreeClusters &c);

};

typedef GTCluster* GTClusterArr;


#endif