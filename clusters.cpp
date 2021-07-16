
#include <string.h>
#include <algorithm>
#include "clusters.h"
#include "rtree.h"

bool compcnt(const GTCluster *a, const GTCluster *b)
{ 
    return a->usagecnt > b->usagecnt; 
} 

bool compsize(const GTCluster *a, const GTCluster *b)
{ 
  return (spsize(a->spcluster) < spsize(b->spcluster)); 
} 


// Only for compatible clusters collections
class GTCC {
public:
  SPID *c;
  string s;
  int clean;
  GTCC(SPID *_c, string _s, bool _clean=0) : c(_c), s(_s), clean(_clean) {}
};


// Generate quasi consensus tree
// If preserveroottree is given, the root split will be takes from the tree
RootedTree* TreeClusters::genrootedquasiconsensus(RootedTree *preserveroottree)
{
  
  vector<GTCluster*> sc,compclusters;
  for (size_t i=0; i<internal.size(); i++) sc.push_back(internal[i]);
  for (size_t i=0;i<leaves.size(); i++) sc.push_back(leaves[i]);
  sort(sc.begin(),sc.end(),compcnt);

  // for (int i=0;i<specnames.size();i++) cout<< i << "-" << specnames[i] << " ";
  //   cout << endl;

  SPID lc=MAXSP,rc=MAXSP;
  if (preserveroottree)
  {
    SPID** t=preserveroottree->getspclusterrepr();
    SPID prroot = preserveroottree->getroot(); 
    preserveroottree->getchild(prroot,rc);
    lc=rc;
    preserveroottree->getchild(prroot,rc);        
    compclusters.push_back(new GTCluster(0,0,t[lc]));
    compclusters.push_back(new GTCluster(0,0,t[rc]));    

    // TODO: Clean array t
  }

  int maxcnt=sc[0]->usagecnt;
  float minusage=0.01*maxcnt;

  // for (int i=0;i<sc.size();i++)
  //   { cout << "cr#" << sc[i]->usagecnt << " cluster=" << *sc[i] << " ";      
  //     if (sc[i]->spcluster[0]<specnames.size()/2 && sc[i]->spcluster[0]>1)
  //       if (sc[i]->spcluster[0]==1)
  //         cout << "+";
  //       cout << endl;

  //     }

    for (size_t i=0; i<sc.size();i++)
    {   
        GTCluster *gc=sc[i];
        if (gc->spcluster[0]>1 && (gc->spcluster[0]>0.7*specnames.size() || 
            gc->usagecnt<minusage)) continue;

        int ok=1;        
        
        for (size_t j=0; j<compclusters.size();j++)
        {

            SPID *cur=compclusters[j]->spcluster;            
            if (cur==topspcluster) continue;
                                    
            
            if (spsubseteq(gc->spcluster,cur) && spsubseteq(cur, gc->spcluster)) { ok=0; break; }

            SPID* sum=joinspclusters(cur,gc->spcluster);  // genrootedquasiconsensus
            
            if (spsize(sum)==spsize(gc->spcluster) 
                || spsize(sum)==spsize(cur) || spsize(sum)==spsize(gc->spcluster)+spsize(cur)) { 
                deletespcluster(sum);
                continue; 
                } // inclusion/disjont

            deletespcluster(sum);

            ok=0;
            break;      
        }

        if (ok) { compclusters.push_back(gc); 
            //cout << "Inserting " << *gc << endl;
            } 
        //else cout << "Failure" << endl;
    }

  // Compatible clusters 
  sort(compclusters.begin(),compclusters.end(),compsize);


  vector<GTCC> vs;
  for (size_t i=0; i<compclusters.size(); i++)   
  {
    GTCluster *c = compclusters[i];

    // cout << "  CUR=" << *c << endl;

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

  
  char *s = strdup(vs[0].s.c_str());
  RootedTree *r = new RootedTree(s); 
  free(s);
  return r; 

}


void TreeClusters::addtree(RootedTree *t)
{  	
	t->setspclusters(this);
}


GTCluster* TreeClusters::get(GTCluster *l, GTCluster *r)
/*
Returns a cluster from two children clusters l and r.
*/
{
  GTCluster *gc;
  
  SPID *s = joinspclusters(l->spcluster, r->spcluster); // OK

  if (s[0]==1) return t[s]; // leaf cluster

  map<SPID*, GTCluster*, comparespids>::iterator it = t.find(s);

  if (it == t.end()) {
    // add new cluster 
    gc = new GTCluster(l, r, s);
    t[s] = gc;
    internal.push_back(gc);
  }
  else {
    //cluster present

    gc = t[s];
    deletespcluster(s); 
  }

  gc->usagecnt++;
  _usagecnt++;
  return gc;
}


