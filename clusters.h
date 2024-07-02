
#ifndef _CLUSTERS_H_
#define _CLUSTERS_H_

#include "tools.h"


class RootedTree;
class Dag;
class Network;

#define deletetablespcluster(s) (((s)[0])>1)
#define deletespcluster(s) if (deletetablespcluster(s)) delete[] s;

typedef struct GTCluster
{
  NODEID type;        // >=0 species id, -1 - all species, -2 - internal l-r
  NODEID* spcluster;  // cluster
  int usagecnt;     
  
  GTCluster(NODEID *s) : type(-2),  spcluster(s), usagecnt(0)  
  { if (s[0] == 1) type = s[1]; }

  friend ostream& operator<<(ostream&s, GTCluster &p) {
    pprintspcluster(s, p.spcluster);
    // s << " #"<< p.usagecnt << " T" << p.type;
    return s;
  }

  int size() const  
  { 
    return ((spcluster[0])==0)?specnames.size():spcluster[0]; 
  }

  ~GTCluster() {  if (deletetablespcluster(spcluster)) delete[] spcluster;  } 

} GTCluster;

/*
 Collection of GTC clusters
*/
class Clusters
{
  int _usagecnt;  
  map<NODEID*, GTCluster*, comparespids> t;
  
public:
  Clusters()
  {    
    _usagecnt = 0;        
  }

 
  ~Clusters()
  {    
      for (auto& pair: t)
      {
        if (deletetablespcluster(pair.second->spcluster))
          delete pair.second;
      }
  }

  GTCluster* has(NODEID *s);
  GTCluster* add(NODEID *s);
  GTCluster* get(NODEID *s); // has + add

  int size() { return t.size(); }
  int usagecnt() { return _usagecnt; }

  // Add clusters from a rooted binary tree
  void addtree(RootedTree *g);

  // Add clusters from a rooted binary dag
  void adddag(Dag *n);
  
  // Add clusters from a newick representation (may be )
  void addtree(char *s);

  GTCluster *_parse(char *s, int &p, int num);

  string genrootedquasiconsensus(RootedTree *preserveroottree, Clusters *guideclusters, Clusters *guidetree, float genetreessimilarity);

  friend ostream& operator<<(ostream&s, Clusters &c);
  GTCluster *addleafcluster(NODEID i);

  bool hasall(Clusters* clusters);

  Clusters(vector<RootedTree*> &genetreesv);
  


};

typedef GTCluster* GTClusterArr;


#endif