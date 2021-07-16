

#ifndef _CLUSTERS_H_
#define _CLUSTERS_H_

#include "tools.h"

class RootedTree;

#define deletespcluster(s) if (s!=topspcluster && ((s)[0]!=1)) delete[] s;

typedef struct GTCluster
{
  GTCluster *l, *r; // child nodes composition; needed to compute lca
  SPID type;        // >=0 species id, -1 - all species, -2 - l-r
  SPID lcamap;      // map to species tree
  SPID* spcluster;  // cluster
  int usagecnt;     
  
  GTCluster(GTCluster *lc, GTCluster *rc, SPID *s) : l(lc), r(rc), type(-2), lcamap(-1), spcluster(s), usagecnt(0)  
  { if (s[0] == 1) type = s[1]; }

  friend ostream& operator<<(ostream&s, GTCluster &p) {
    printspcluster(s, p.spcluster);
    return s;
  }

  ~GTCluster() {  deletespcluster(spcluster);  }

} GTCluster;

class TreeClusters
{
  int ssize;
  int _usagecnt;  
  map<SPID*, GTCluster*, comparespids> t;
  
  vector<GTCluster*> internal;
  vector<GTCluster*> leaves;

public:
  TreeClusters()
  {
    ssize = specnames.size();
    _usagecnt = 0;
    // initialize singleton clusters

    for (SPID i = 0; i < ssize; i++)
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

    
  GTCluster* leafcluster(SPID n) 
  { 
    leaves[n]->usagecnt++; 
    _usagecnt++; 
    return leaves[n]; 
  }

  GTCluster* get(GTCluster *l, GTCluster *r);
  int size() { return t.size(); }
  int usagecnt() { return _usagecnt; }

  void addtree(RootedTree *g);

  RootedTree* genrootedquasiconsensus(RootedTree *preserveroottree);

};

typedef GTCluster* GTClusterArr;


#endif