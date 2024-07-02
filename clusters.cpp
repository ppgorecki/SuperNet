
#include <string.h>
#include <algorithm>
#include "clusters.h"
#include "rtree.h"
#include "tools.h"

extern double opt_genetreessimilarity;

bool compcnt(const GTCluster *a, const GTCluster *b)
{ 
    return a->usagecnt > b->usagecnt; 
} 

bool compsize(const GTCluster *a, const GTCluster *b)
{ 
  return a->size() < b->size(); 
} 

// Only for compatible clusters collections
class GTCC {
public:
  NODEID *c;
  string s;
  int clean;
  GTCC(NODEID *_c, string _s, bool _clean=0) : c(_c), s(_s), clean(_clean) {}
};

bool propercompatiblecluster(GTCluster *candcluster, GTCluster *basecluster)
{
    // top cluster 
    if (basecluster->spcluster==topspcluster) return true;

    // check if equal (redundant)
    if (spsubseteq(candcluster->spcluster, basecluster->spcluster) && spsubseteq(basecluster->spcluster, candcluster->spcluster)) return false; 

    NODEID* sum = joinspclusters(basecluster->spcluster, candcluster->spcluster);  

    int sumsize = spsize(sum);
    int candclustersize = spsize(candcluster->spcluster);
    int baseclustersize = spsize(basecluster->spcluster);
        
    deletespcluster(sum);
    
    return sumsize==candclustersize || sumsize==baseclustersize || sumsize==(candclustersize+baseclustersize);
      // inclusion/disjont

}


// Generate quasi consensus tree
// If preserveroottree is given, the root split will be taken from the tree
string Clusters::genrootedquasiconsensus(RootedTree *preserveroottree, Clusters *guideclusters, 
  Clusters *guidetree, float genetreessimilarity)
{
  
  vector<GTCluster*> candidates, compclusters;  

  if (!specnames.size())
  {
      cerr << ERR_NOSPECIESDEFINED << endl;
      exit(-1);
  }

  for (NODEID i = 0; i < specnames.size(); i++)
  {    
    get(spec2gtcluster[i]->spcluster);
  }
  
  for (auto& pairObj: t)
  {
      candidates.push_back(pairObj.second);  
  }
    
  sort(candidates.begin(), candidates.end(), compcnt);

  if (preserveroottree)
  {
    NODEID** t=preserveroottree->getspclusterrepr(); 
    NODEID prroot = preserveroottree->getroot(); 
    NODEID rc = MAXNODEID;

    preserveroottree->getchild(prroot,rc);
    NODEID lc=rc;
    
    preserveroottree->getchild(prroot,rc);        
    
    compclusters.push_back(new GTCluster(t[lc])); 
    compclusters.push_back(new GTCluster(t[rc]));    
    // TODO: Clean array t
  }


  int maxcnt=1;
  if (candidates.size()>0) 
    maxcnt = candidates[0]->usagecnt;

  float minusage=0.01*maxcnt;

#ifdef _CLUDEBUG_
  for (int i=0;i<candidates.size();i++)
    { 
      cout << i << " cr#" << candidates[i]->usagecnt << " cluster=" << *candidates[i] << " ";      
      if (candidates[i]->spcluster[0]<specnames.size()/2 && candidates[i]->spcluster[0]>1)
        if (candidates[i]->spcluster[0]==1)
          cout << "+";
        cout << endl;
    }
#endif 

  if (guideclusters)
  {
    // add all guide clusters (check if OK)

    for (auto &candpair: guideclusters->t)
    {    
      GTCluster *candcluster = candpair.second;
      if (candcluster->size()==1) continue;
      if (candcluster->size()==specnames.size()) continue;

      bool skip = false;
      // check compatiblity
      for (auto &comp: compclusters)
      {
        if (spsubseteq(candcluster->spcluster, comp->spcluster) && spsubseteq(comp->spcluster, candcluster->spcluster))
        {
          skip = true; // equal, already present, skip
          break;
        }

        if (!propercompatiblecluster(candcluster, comp))
            {
              cerr << "Ooops. Incompatible clusters in a guide tree and/or in preserveroottree opt? " << endl;
              exit(-1);
            } 
      }
      if (!skip)
      {
          // add cluster
          compclusters.push_back(candcluster); 
      }
    }
  }


  // Ugly repeat!
  if (guidetree)
  {
    // add all guide tree clusters (check if OK)

    for (auto &candpair: guidetree->t)
    {    
      GTCluster *candcluster = candpair.second;
      if (candcluster->size()==1) continue;
      if (candcluster->size()==specnames.size()) continue;

      bool skip = false;
      // check compatiblity
      for (auto &comp: compclusters)
      {
        if (spsubseteq(candcluster->spcluster, comp->spcluster) && spsubseteq(comp->spcluster, candcluster->spcluster))
        {
          skip = true; // equal, already present, skip
          break;
        }

        if (!propercompatiblecluster(candcluster, comp))
            {
              cerr << "Ooops. Incompatible clusters in a guide tree/clusters and /or preserveroottree opt? " << endl;
              exit(-1);
            } 
      }
      if (!skip)
      {
          // add cluster
          compclusters.push_back(candcluster); 
      }
    }
  }

  for (auto &cand: candidates)
  {       
        // randomize     

        if (cand->spcluster[0]>1 && (cand->spcluster[0]>genetreessimilarity*specnames.size() || 
            cand->usagecnt<minusage)) 
        {
            // ignore
            continue;
        }

        int ok=1;        
        
        for (auto &comp: compclusters)
        {
            if (!propercompatiblecluster(cand, comp))
            {
              ok=0;
              break;
            }                                
        }

        if (ok)       
            compclusters.push_back(cand); 
                 
    }

  // Compatible clusters 
  sort(compclusters.begin(),compclusters.end(),compsize);


  vector<GTCC> vs;
  for (size_t i=0; i<compclusters.size(); i++)   
  {
    GTCluster *c = compclusters[i];

    if (c->spcluster[0]==1) // leaf
    {
       vs.push_back(GTCC(c->spcluster, species(c->spcluster[1]))); 
       continue;
    }

    vector<GTCC> res; 
    for (int j = 0; j < (int)vs.size(); j++)
      if (spsubseteq(vs[j].c, c->spcluster))
      {
        res.push_back(vs[j]);
        vs.erase(vs.begin() + j);
        j--;
      }

    
    //for (int j=0;j<res.size();j++ ) { cout << "z "; printspcluster(cout,res[j].c); cout << endl; }

    if (res.size() < 2)
    {
      cerr << "At least two subsets in genrootedquasiconsensus expected..." << endl;
      return NULL;
    }

    while (res.size()>1)
    {
        int a = rand()%res.size();
        GTCC l = res[a];
        res.erase(res.begin() + a);

        a = rand() % res.size();
        GTCC r = res[a];
        res.erase(res.begin() + a);

        res.push_back(GTCC(joinspclusters(l.c,r.c), string("(" + l.s + "," + r.s + ")"))); // genrootedquasiconsensus

        deletespcluster(l.c);
        deletespcluster(r.c);

    }    

    vs.push_back(res[0]);      
  }

  // Ugly
  while (vs.size()>1)
    {
        int a = rand()%vs.size();
        GTCC l = vs[a];
        vs.erase(vs.begin() + a);

        a = rand() % vs.size();
        GTCC r = vs[a];
        vs.erase(vs.begin() + a);

        vs.push_back(GTCC(joinspclusters(l.c,r.c), string("(" + l.s + "," + r.s + ")"))); // genrootedquasiconsensus

        deletespcluster(l.c);
        deletespcluster(r.c);
    }    

  return vs[0].s;
}

void Clusters::adddag(Dag *dag)
{      
    dag->getclusters(this);    
}

GTCluster* Clusters::has(NODEID *s)
{
    map<NODEID*, GTCluster*, comparespids>::iterator it = t.find(s);
    if (it == t.end()) 
        return NULL;
    return t[s];
}

GTCluster* Clusters::add(NODEID *s)
{
  GTCluster *gc;
  if (s[0]==1)
    gc = spec2gtcluster[s[1]];
  else gc = new GTCluster(s);
  t[s] = gc;  
  return gc;  
}


GTCluster* Clusters::get(NODEID *s)
{

    GTCluster *gc = has(s);

    if (!gc)
    {
      gc = add(s);
    }

    gc->usagecnt++;
    _usagecnt++;

    return gc;
}

ostream& operator<<(ostream&s, Clusters &c)
{
  for (auto& pairObj: c.t)
  {      
      //s << std::cout<< *pairObj.first << " -> " << *pairObj.second << std::endl;
      s << *pairObj.second << endl;
  }
  return s;
}

void Clusters::addtree(char *s)
{
  int p=0;  
  while (s[p])
  {
    _parse(s, p, 0);    
    if (s[p]==';') p++;
    if (!s[p]) break;    
  }
}


GTCluster *Clusters::_parse(char *s, int &p, int num)
{
  char *basecluster = getTok(s, p, num);
  char *token; 

  NODEID *spclu = NULL;

  char *st = s+p;

  if (basecluster[0] == '(')
  {

    while (1)
    {
      GTCluster *a = _parse(s, p, num);

      if (!spclu) spclu = a->spcluster;
      else
      {
        // TODO: clean
        spclu = joinspclusters(spclu, a->spcluster, NULL, false); 
      }    
      
      token = getTok(s, p, num);

      if (token[0] != ',')
        break;                    
    }

    if (!has(spclu))
      return add(spclu); 
  }

  NODEID spid = getspecies(basecluster, s + p - basecluster);        
  return get(spec2gtcluster[spid]->spcluster);
}

bool Clusters::hasall(Clusters* clusters)
{
    if (clusters->size()>size()) return false;
    for (auto& pairObj: clusters->t)
    {
        if (!has(pairObj.second->spcluster)) return false;
    }
    return true;
}


Clusters::Clusters(vector<RootedTree*> &genetreesv)
{
  for (auto gtree: genetreesv) 
  {
     adddag(gtree);    
  }
}
