/************************************************************************
 SuperNetwork Inference - based on FastUrec v2.05 
(c) Copyright 2005-2020 by Pawel Gorecki
 Written by P.Gorecki.
 Permission is granted to copy and use this program provided no fee is
 charged for it and provided that this copyright notice is not removed.
 *************************************************************************/

const char *SUPNET="0.01";

#include <time.h> 
#include <sys/time.h>
 
#define TRNEW 0
#define TRPROCESSED 1

#define COSTDUPLICATIONLOSS 0   // dl
#define COSTDUPLICATION 1    // duplication
#define COSTDEEPCOAL 2   // deep coalescence
#define COSTROBINSONFOULDS 3   // robinson-foulds
#define COSTLOSS 4    // loss

int _gsetrees=0;

#define SPC int


#define LCATAB // If set lca queries are stored in O(n^2) array

#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>

#include <string.h>
#include <iostream>
#include <set>
#include <vector>
#include <algorithm>
using namespace std;

#include "tools.h"

#include <sstream>
#include <queue>


#define pcc(d,v)
#define pcct(d,v)

int randomtreescnt = 0;
int quasiconsensuscnt = 0;

char *projectname = strdup("byopt");

#define VERBOSE_QUIET 1
#define VERBOSE_CNTCOST 2
#define VERBOSE_SMP 3
#define VERBOSE_DETAILED 4

#define verb_cntcost if (verbose==VERBOSE_CNTCOST) cout
#define verb_smp if (verbose>=VERBOSE_SMP) cout
#define verb_det if (verbose>=VERBOSE_DETAILED) cout

#define deletespcluster(s) if (s!=topspcluster && ((s)[0]!=1)) delete[] s;

int verbose = VERBOSE_SMP; 

template <class T>
inline std::string to_string (const T& t)
{
  std::stringstream ss;
  ss << t;
  return ss.str();
}

double weightloss = 1.0;
double weightdup = 1.0;

int gsid = GSFULL;
char *gsdelim = NULL;
int gspos = 0;

int lfnum(char *str)
{  
  int cnt = 0;
  
  cnt=1;
  while (*str) 
  {
    if (*str == ',') cnt++;
    str++;
  }
  return cnt;

}

#define LEAFLOOP for (SPID i=0; i<lf; i++)
#define INTERNALNODESLOOP for (SPID i=lf; i<nn; i++)
#define NODESLOOP for (SPID i=0; i<nn; i++)


#define isinternalnode(i) ((i)>=lf)
#define isleaf(i) ((i)<lf)

#define DEFBL -100.0   // def. branch length in a gene tree

class RootedTree;

typedef vector<RootedTree*> VecRootedTree;

struct comparespids
{
  bool operator()(const SPID* a, const SPID* b) const
{
      int i = 1, al = a[0] + 1, bl = b[0] + 1;

      if (al<bl) { //cout << "Tx";
	return true; }
      if (al>bl) { //cout << "Fx";
	return false; }
    for ( ; i < al; i++ ) {
      if (a[i] < b[i]) { //cout << "T";
	return true;}
      if (a[i] > b[i]) { //cout << "F";
	return false; }
    }
      //    cout <<"F";
    return false;
	
  }
};

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

  ~GTCluster()
  {
    deletespcluster(spcluster);
  }

} GTCluster;

typedef GTCluster* GTClusterArr;

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

    
  GTCluster* leafcluster(SPID n) { leaves[n]->usagecnt++; _usagecnt++; return leaves[n]; }

  GTCluster* get(GTCluster *l, GTCluster *r);
  int size() { return t.size(); }
  int usagecnt() { return _usagecnt; }

  void addtree(RootedTree *g);



  RootedTree* genrootedquasiconsensus(RootedTree *preserveroottree);

};

TreeClusters *gtc;

typedef set<SPID> cluster;
typedef map<cluster, int> cntmap;

bool compcnt(const GTCluster *a, const GTCluster *b)
  { return a->usagecnt > b->usagecnt; } 

bool compsize(const GTCluster *a, const GTCluster *b)
  { return (spsize(a->spcluster) < spsize(b->spcluster)); } 


// Only for compatible clusters collections
class GTCC {
public:
  SPID *c;
  string s;
  int clean;
  GTCC(SPID *_c, string _s, bool _clean=0) : c(_c), s(_s), clean(_clean) {}
};



////////////////////////////////////////////////////////
// generic tree
class Tree
{
protected:

  bool exactspecies; // true iff for each leaf lab[leaf]==leaf

  void init(int _lf, int _it, char *s = NULL)
  {
    lf = _lf;
    it = _it;
    nn = lf + it;
    //cerr << lf << " " << it << " " << nn << " " << MAXSP <<  endl;
    if (2*lf >= MAXSP)
    {
      cout << "Tree is too large. " << lf << " leaves." << endl;
      cout << "Compile with SPLARGE or SPMED macro." << endl;
      exit(-1);
    }

    
 
    if (nn>0) par = new SPID[nn];
    else par=NULL;
    lab = new SPID[lf];    
  }

  virtual ostream& _printdebnode(ostream&s, SPID i)
  {
    s << setw(3) << (SPC)i;
    if (par[i] == MAXSP) s << "-  ";
    else s << setw(3) << (SPC)par[i];

    if (i < lf)
      s << " $" << (SPC)lab[i] << " " << species(lab[i]);
    else s << " *   ";
    return s;
  }

  void setexactspecies()
  {
    exactspecies = 1;
    LEAFLOOP { if (lab[i] != i) { exactspecies = 0; break; } }
  }


public: // GenericTree
  string name;
  double weight;
  int lf;
  int it;
  int nn;
  SPID *par; 
  SPID *lab;  
  

  virtual ~Tree()
  {
    if (par) delete[] par;
    delete[] lab;
  }

  virtual ostream& printdebnode(ostream&s, SPID i)
  {
    s << setiosflags(ios::left);
    _printdebnode(s, i) << resetiosflags(ios::left);
    return s;
  }



  virtual ostream& printsubtree(ostream&s, SPID i) { return s; }

  Tree(double treeweight = 1.0) : name(""), weight(treeweight), lf(0), it(0), nn(0)
  {}

  // Tree
  virtual ostream& printdeb(ostream&s, int gse, string tn="") {
    if (gse&2)
    {
      s << "leaves=" << lf << " internal=" << it << " total=" << nn << " exactspecies=" << exactspecies << endl;
      NODESLOOP printdebnode(s, i) << endl; 


    }
    if (gse&1)
    {
      s << "&s (";
      char cnt='A';
      NODESLOOP if (par[i]==MAXSP)
      {    
           if (cnt!='A') s << ",";          
           printsubtree(s,i) << " edgelab=\"" << cnt << "\"" ;
           cnt++;
      }
      s << ")" << " treename=\"" << _gsetrees++ << " " << tn << "\"" ;
    }
    return s;
  }

  virtual ostream& print(ostream&s)  {
    s << lf << it << nn << endl;
    NODESLOOP printdebnode(s, i) << endl;
    return s;
  }
}; // Tree


////////////////////////////////////////////////////////
// rooted binary tree
class RootedTree : public Tree
{
protected:

  // Init rooted tree structures
  void init(int _lf, int _it, char *s = NULL)
  {
    Tree::init(_lf, _it, s);
    rightchild = new SPID[nn];
    leftchild = rightchild - lf; // proper indexing 

    level = new SPID[nn];
    depth = NULL;
  }


  // Set depth of nodes in a subtree (recursive)
  void _setdepth(SPID i, int dpt)
  {
    depth[i] = dpt;
    if (isleaf(i)) return;
    _setdepth(leftchild[i], dpt + 1);
    _setdepth(rightchild[i], dpt + 1);
  }

  // Insert clusters of children node n into TreeClusters
  GTCluster* _getgtccluster(SPID n, TreeClusters *gtc)
  {
    if (gtcclusters[n]) return gtcclusters[n];    
    GTCluster *r1 = _getgtccluster(leftchild[n], gtc);
    GTCluster *r2 = _getgtccluster(rightchild[n], gtc);
    return gtcclusters[n] = gtc->get(r1, r2);
  }


/* 

Parse tree from a char*

Returns id of the root of the subtree that starts at s+p

Tree is represented in prefix order:
> supnet -g '(a,((b,a),c))'  -eD
leaves=4 internal=3 total=7 exactspecies=1

i   par[i]  lab[i]  LAB      LFT RGH child  subtree
0     6       $0    a                       a
1     4       $1    b                       b
2     4       $0    a                       a
3     5       $2    c                       c
4     5       *               1   2         (b,a)
5     6       *               4   3         ((b,a),c)
6     MAXSP   *         Root  0   5         (a,((b,a),c))

 
  leftchild=[ 1 4 0 ]
 rightchild=[ 2 3 5 ]
        par=[ 6 4 4 5 5 6 126 ]
      level=[ 0 0 0 0 1 2 3 ]


Variables:
  lf - the number of leaves
  nn - the number of all nodes
  0..lf-1 - ids of leaves
  lf..nn-1 - ids of internal nodes
  root id == nn-1

Arrays:
  lab[0..lf-1] - id of the species (see getspecies)
  par[0..nn-1] - parent
  leftchild[lf..nn-1] - children maps; don't use indexes <lf
  rightchild[lf..nn-1] - children maps; don't use indexes <lf

  leftchild[i]  - left child of i
  rightchild[i] - right child of i

*/


  SPID parse(char *s, int &p, int num, SPID &freeleaf, SPID &freeint)
  {
    char *cur = getTok(s, p, num);

    if (cur[0] == '(')
    {
      int a = parse(s, p, num, freeleaf, freeint);
      getTok(s, p, num);
      int b = parse(s, p, num, freeleaf, freeint);
      getTok(s, p, num);
      if (freeint >= nn)
      {
        cerr << "Is it a binary tree? Too many nodes in rooted tree" << endl;
        exit(-1);
      }
      rightchild[freeint] = b;
      leftchild[freeint] = a;
      par[a] = freeint;
      par[b] = freeint;
      level[freeint] = max(level[a], level[b]) + 1;
      return freeint++;
    }
    // leaf
    if (freeleaf >= lf)
    {
      cerr << "Is it a binary tree? Too many leaves in rooted tree" << endl;
      exit(-1);
    }
    lab[freeleaf] = getspecies(cur, s + p - cur);
    level[freeleaf] = 0;
    return freeleaf++;
  }

  // Print node info
  virtual ostream& _printdebnode(ostream&s, SPID i)
  {
    Tree::_printdebnode(s, i);
    if (root == i) s << "Rt"; else s << "  ";
    if (isinternalnode(i))
      s << " " << setw(3) << (SPC)leftchild[i] << " " << setw(3) << (SPC)rightchild[i];    
    printsubtree(s, i);
    return s;
  }

  // Set levels of nodes
  void setlevel(int i)
  {
    if (isleaf(i)) level[i] = 0;
    else
    {
      setlevel(leftchild[i]);
      setlevel(rightchild[i]);
      level[i] = 1 + max(level[leftchild[i]], level[rightchild[i]]);
    }
  }

public :
  SPID root, *level, *leftchild, *rightchild, *depth;  
  
  GTCluster **gtcclusters; // node clusters

#ifdef LCATAB
  SPID **lcatab;
#endif

RootedTree(char *s, double weight = 1.0): Tree(weight)
{
  
    int lfn = lfnum(s);
    init(lfn, lfn - 1, s);
    int p = 0;
    SPID freeleaf = 0;
    SPID freeint = lf;
    while (s[0]>='0' && s[0]<='9') s++; // skip optional first number
    root = parse(s, p, 0, freeleaf, freeint);
    checkparsing(s + p);

    setexactspecies();

    depth = NULL;
#ifdef LCATAB
  lcatab = NULL;
#endif
    par[root] = MAXSP;
    depth = new SPID[nn];
  }

  virtual ~RootedTree()
  {
    if (depth) delete[] depth;
    delete[] level;
    delete[] rightchild; // don't delete leftchild!
  }

  // Subtree nodes generator
  void subtreenodes(SPID i, SPID t[], int &cnt)
  {
    t[cnt++] = i;
    if (isleaf(i)) return;
    subtreenodes(leftchild[i], t, cnt);
    subtreenodes(rightchild[i], t, cnt);
  }
  
  // Print tree
  virtual ostream& print(ostream&s) {
    printsubtree(s, root);
    return s;
  }

  // Get children ids
  virtual void getchildren(SPID i, SPID &l, SPID &r)
  {
    l = leftchild[i];
    r = rightchild[i];
  }

  // Insert clusters of a tree into TreeClusters
  void setspclusters(TreeClusters *gtc)
  {
    if (nn<1) return;
    gtcclusters = new GTClusterArr[nn];
    NODESLOOP gtcclusters[i] = NULL;
    LEAFLOOP gtcclusters[i] = gtc->leafcluster(lab[i]);
    INTERNALNODESLOOP _getgtccluster(i, gtc);        
  }


  double cost(RootedTree &speciestree, int costfunc)
  {
    SPID *lcamap = getlcamapping(speciestree);
    speciestree.setdepth();
    if (costfunc == COSTLOSS) return costloss(speciestree,lcamap);
    if (costfunc == COSTDUPLICATIONLOSS) return costduplicationloss(speciestree,lcamap);
    if (costfunc == COSTDUPLICATION) return costduplication(speciestree,lcamap);
    if (costfunc == COSTDEEPCOAL) return costdeepcoal(speciestree,lcamap);
    if (costfunc == COSTROBINSONFOULDS) return costrobinsonfoulds(speciestree,lcamap);

    cerr << "Unknown costfunc type?" << endl;
    exit(-1);

  }

  long costloss(RootedTree &speciestree, SPID *lcamap)
  {
      return costdeepcoal(speciestree,lcamap)+2*costduplication(speciestree,lcamap);
  }

  long costrobinsonfoulds(RootedTree &speciestree, SPID *lcamap)
  {
      //TODO
      return 0;
  }

  long costduplicationloss(RootedTree &speciestree, SPID *lcamap)
  {
      return costdeepcoal(speciestree,lcamap)+3*costduplication(speciestree,lcamap);
  }

  long costduplication(RootedTree &speciestree, SPID *lcamap)
  {      
      long s = 0;
      INTERNALNODESLOOP { if (lcamap[leftchild[i]]==lcamap[i] || lcamap[rightchild[i]]==lcamap[i]) s++; }
      return s;
  }

  long costdeepcoal(RootedTree &speciestree, SPID *lcamap)
  {      
      long s = 0;            
      NODESLOOP { if (i!=root) s+=speciestree.depth[lcamap[i]]-speciestree.depth[lcamap[par[i]]]-1; }
      return s;
  }


  // Returns true if the labelling is bijective
  bool isspeciestree()
  {            
      SPID occ[specnames.size()];
      for (size_t i=0; i<specnames.size(); i++) occ[i]=0;
      for (int i=0; i<lf; i++) 
        if (occ[lab[i]]++==1) return false;
      return true;
  }


  virtual int subtreesize(SPID n)
  {
    if (isleaf(n)) return 1;
    return subtreesize(leftchild[n]+subtreesize(rightchild[n]));
  }


  // Get all clusters of nodes
  // Now used only to compute the root cluster in quasiconsensus (todo: replace it by GTClusters)
  SPID** getspclusterrepr()
  {
        SPID **res=new SPID*[nn];
        LEAFLOOP res[i]=spec2spcluster[lab[i]];
        INTERNALNODESLOOP res[i]=joinspclusters(res[leftchild[i]],res[rightchild[i]]);
        return res;
  }

  void saveraw(string fn = "st.raw")
  {
    ofstream s;
    s.open(fn.c_str());
    s << specnames.size() << endl;
    for (int i = 0; i < (int)specnames.size(); i++) s << species(i) << " ";
    s << endl;
    s << nn << " " << lf << " ";
    s << endl;
    for (int i = 0; i < nn; i++) s << par[i] << " ";    
    s << endl;
    for (int i = 0; i < lf; i++) s << lab[i] << " ";     
    s << endl;
    for (int i = lf; i < nn; i++) s << leftchild[i] << " ";     
    s << endl;
    for (int i = lf; i < nn; i++) s << rightchild[i] << " ";     
    s << endl;
    s << root;
    s << endl;
    s.close();
  }

  // size of the cluster including duplicates
  SPID multiclustersize(SPID n) 
  {
    if (isleaf(n)) return 1;
    return multiclustersize(leftchild[n]) + multiclustersize(rightchild[n]);
  }



  int _lca(SPID a, SPID b)
  {
    int lev = min(level[a],level[b]); 
    while (1)
    {
      if (a == b) return a;
      if (a == root) return root;
      if (b == root) return root;
      if (level[a] == lev) a = par[a];
      if (level[b] == lev) b = par[b];
      lev++;
    }
  }

  friend ostream& operator<<(ostream&s, RootedTree &p)  { return p.print(s); }

  void initlca()
  {
    setlevel(root);
#ifdef LCATAB
    initlcatab();
#endif
  }

#ifdef LCATAB
  void initlcatab()
  {
    if (!lcatab) {
      lcatab = new SPID*[nn];
      for (int i = 1; i < nn; i++) lcatab[i] = new SPID[i];
    }
    for (int i = 1; i < nn; i++)
      for (int j = 0; j < i; j++)
        lcatab[i][j] = _lca(i, j);
  }
#endif

  int lca(SPID a, SPID b)
  {
#ifdef LCATAB
    if (a == b) return a;
    if (!lcatab) initlca();
    if (a < b) return lcatab[b][a];
    return lcatab[a][b];
#else
    return _lca(a, b);
#endif
  }

  virtual ostream&  printsubtree(ostream &s, SPID i)
  {

#ifdef GSESUBTREES

    if (i == MAXSP) s << "---";
    else    
    if (isleaf(i)) s << species(lab[i]);
    else
    {
        s << "(";
        printsubtree(s, leftchild[i]) << ",";
        printsubtree(s, rightchild[i]) << ")";
    }
    return s<< " nodeinfo=\"" << i << "\" edgelab=\"" << par[i] << "\"";

#else
 if (i == MAXSP) return s << "---";
    if (isleaf(i)) return s << species(lab[i]);
    s << "(";
    printsubtree(s, leftchild[i]) << ",";
    return printsubtree(s, rightchild[i]) << ")";
#endif
  }

  // Requirement: a <= b
  inline virtual SPID dist(SPID a, SPID b)
  {
    return depth[a]-depth[b];    
  }

  inline virtual SPID nodedist(SPID a, SPID b)
  {
    SPID lcaab=lca(a,b);
    return dist(a,lcaab)+dist(b,lcaab);
  }

  // Used in sorted representation of tree
  virtual void _repr(SPID *t, SPID i, int &pos, SPID* mm)
  {
    if (isleaf(i))
    {
      t[pos++] = lab[i];
      return;
    }
    t[pos++] = MAXSP; // (

    if (mm[leftchild[i]] == mm[i])
    {
      _repr(t, leftchild[i], pos, mm);
      t[pos - 1] = -t[pos - 1]; //,
      _repr(t, rightchild[i], pos, mm);
    }
    else
    {
      _repr(t, rightchild[i], pos, mm);
      t[pos - 1] = -t[pos - 1]; //,
      _repr(t, leftchild[i], pos, mm);
    }
    t[pos++] = MAXSP + 1; // )
  }

  // Used in sorted representation of tree to compute 
  // a node in a subtree with min spid
  virtual SPID _setmm_repr(SPID i, SPID* mm)
  {
    if (isleaf(i)) mm[i] = lab[i];
    else
    {
      SPID lm = _setmm_repr(leftchild[i], mm);
      SPID rm = _setmm_repr(rightchild[i], mm);
      mm[i] = min(lm, rm);
    }
    return mm[i];
  }

  // The length of an array in sorted representation 
  int reprlen() 
  {
    return  lf * 3 - 2;
  }

  // Generate unique representation of a tree as a vector
  // (sorted representation)
  virtual SPID* repr(SPID *t = NULL)
  {
    if (!t) t = new SPID[reprlen()];
    SPID mm[nn];
    _setmm_repr(root, mm);
    int cpos = 0;
    _repr(t, root, cpos, mm);
    return t;
  }


  // Find leaf label
  SPID findlab(SPID slab, int stoponerr=1)
  {
    if (exactspecies && slab < lf) return slab;
      
    LEAFLOOP if (slab == lab[i]) return i;
    if (stoponerr)
    {
      cerr << "Species lab " << specnames[slab] << " ($" << (int)slab << ") not found in the species tree: ";
      print(cerr);
      cerr << endl;
      exit(-1);
    }
    return MAXSP;
  }


  void setdepth(SPID i = -1, SPID dpt = 0)
  {    
    _setdepth(root, 0);
  }

  int getlosses(SPID s, SPID s1, SPID s2)
  {
    if (s != s2)
    {
      if (s == s1)
      {
        return (int)depth[s2] - (int)depth[s];
      }
      else
      {
        return (int)depth[s1] + (int)depth[s2] - 2 * (int)depth[s] - 2;
      }
    }
    return (int)depth[s1] - (int)depth[s];
  }

  SPID getparent(SPID n) { return par[n]; }


  SPID* getlcamapping(RootedTree &speciestree)
  {
    if (nn<1) return NULL;
    SPID* lcamap = new SPID[nn];
    NODESLOOP lcamap[i] = -1;
    LEAFLOOP lcamap[i] = speciestree.findlab(lab[i]);
    INTERNALNODESLOOP lcamap[i] = speciestree.lca(lcamap[leftchild[i]], lcamap[rightchild[i]]);
    return lcamap;
  }

  virtual ostream& printdeb(ostream&s, int gse, string tn="")
  {
    Tree::printdeb(s, gse, tn);
    s << "  leftchild=[ ";
    INTERNALNODESLOOP s << (int)leftchild[i] << " ";
    s << "]" << endl;
    s << " rightchild=[ ";
    INTERNALNODESLOOP s << (int)rightchild[i] << " ";
    s << "]" << endl;        
    s << "        par=[ ";
    NODESLOOP s << (int)par[i] << " ";
    s << "]" << endl;        
    s << "      level=[ ";
    NODESLOOP s << (int)level[i] << " ";
    s << "]" << endl;    
    
    return s;
  }


}; // RootedTree
////////////////////////////////////////////////////////

#define BUFSIZE 100000

// Read trees from a given file
void readtrees(char *fn, vector<char*> &stvec)
{
  FILE *f;
  f = fopen(fn, "r");
  if (!f)
  {
    cerr << "Cannot open file " << fn << endl;
    exit(-1);
  }
  char buf[BUFSIZE];
  while (1)
  {
    if (!fgets(buf, BUFSIZE, f)) break;
    stvec.push_back(strdup(buf));
  }
  fclose(f);
}


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

  SPID lc,rc;
  if (preserveroottree)
  {
    SPID** t=preserveroottree->getspclusterrepr();
    preserveroottree->getchildren(preserveroottree->root,lc,rc);    
    
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

void TreeClusters::addtree(RootedTree *rt)
{  
    rt->setspclusters(this);
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


// inline bool isInteger(const std::string & s)
// {
//    if(s.empty() || ((!isdigit(s[0])) && (s[0] != '-') && (s[0] != '+'))) return false ;

//    char * p ;
//    strtol(s.c_str(), &p, 10);

//    return (*p == 0) ;
// }


RootedTree *randspeciestree()
{
  SPID tr[specnames.size()];
  for (int i = 0; i < (int)specnames.size(); i++) tr[i] = i;
  return new RootedTree(strdup(genrandomtree(tr, specnames.size()).c_str()));
}

// Read trees from a string; trees are separated by semicolon
void insertrtrees(vector<char*> &v,const char *t)
{
  char *x, *y = strdup(t);
  x=y;

  while (strlen(x)>0)
  {
    char *p=strchr(x,';');
    if (p)
    {      
      *p=0;      
      v.push_back(strdup(x));
      x=p+1;
    }
    else
    {
      v.push_back(strdup(x));      
      break;
    }
  }
  free(y);
}

int main(int argc, char **argv) 
{
  
  int opt;

  if (argc < 2) usage(argc, argv);

  VecRootedTree stvec, gtvec;
  
  struct timeval time; 
  gettimeofday(&time,NULL);
     
  srand((time.tv_sec * 1000) + (time.tv_usec / 1000));
  
  int costfunc = COSTDUPLICATIONLOSS;
  
  int OPT_PRINTQUASICONSTREES = 0;
  int OPT_PRINTRANDSPECIESTREES = 0;
  int OPT_PRINTGENE = 0;
  int OPT_PRINTSPECIES = 0;
  int OPT_PRINTCOST = 0;
  int OPT_PRINTDETAILED = 0;
  int OPT_PRINTINFO=0;  
  
  int OPT_PRESERVEROOT=0;

  
  
  const char* optstring = "e:g:s:G:l:w:L:D:C:r:S:";
  vector<char*> sgtvec, sstvec;
  
  while ((opt = getopt(argc, argv, optstring))   != -1)
  {    
    switch (opt) {

      case 'e':        
        if (strchr(optarg,'i')) OPT_PRINTINFO=1;
        if (strchr(optarg,'g')) OPT_PRINTGENE = 1;
        if (strchr(optarg,'s')) OPT_PRINTSPECIES = 1;
        if (strchr(optarg,'r')) OPT_PRESERVEROOT = 1;
        if (strchr(optarg,'c')) OPT_PRINTCOST = 1; // just numbers
        if (strchr(optarg,'C')) OPT_PRINTCOST = 2; // with trees
        if (strchr(optarg,'D')) OPT_PRINTDETAILED = 1;

        //if (strchr(optarg,'i')) OPT_USERSPTREEISSTARTING=0;        
        //if (strchr(optarg,'S')) OPT_SHOWSTARTINGTREES=1;
        break;
    
    case 'v': 
      verbose=atoi(optarg);
      break;

    case 'g':    
      insertrtrees(sgtvec,optarg);                
      break;
    
    case 's':
      insertrtrees(sstvec,optarg);                    
      break;


    case 'G': {      
      readtrees(optarg, sgtvec);
      string oa(optarg);      
      string bf = oa.substr(oa.find_last_of("/\\") + 1);      
      if (bf.find_last_of(".")!=string::npos)
        if (bf.substr(0,bf.find_last_of("."))==".txt")
          bf=bf.substr(0,bf.find_last_of("."));
      projectname=strdup(bf.c_str());
      
      break;
    }

    // Matching gene tree and species tree labels
    case 'l':
    {
      switch (optarg[0])
      {
      case 'p':
        gspos = atoi(optarg + 1);
        gsid = GSPOS;
        break;
      case 'a':
        gsid = GSAFTER;
        gsdelim = strdup(optarg + 1);
        break;
      case 'b':
        gsid = GSBEFORE;
        gsdelim = strdup(optarg + 1);
        break;
      default:
        cerr << "Invalid code for -l option expected aDELIM,bDELIM or p[-]number" << endl;
        exit(-1);
      }
      break;
    }

    // Print quasi-consensus trees
    case 'w':
      OPT_PRINTQUASICONSTREES = 1;
      if (sscanf(optarg, "%d", &quasiconsensuscnt) != 1) {
        cerr << "Number expected in -w" << endl;
        exit(-1);
      }
      
      break;

    case 'r':
      OPT_PRINTRANDSPECIESTREES = 1;      
      if (sscanf(optarg, "%d", &randomtreescnt) != 1) {
        cerr << "Number expected in -r" << endl;
        exit(-1);
      }
      
      break;

    // define species a,b,c,...
    case 'S':
    {
      int spcnt=0;
      if (sscanf(optarg, "%d", &spcnt) != 1) {
        cerr << "Number expected in -S" << endl;
        exit(-1);
      }
      if (spcnt>'z'-'a'+1)
      {
        cerr << "Too many species in -S" << endl;
        exit(-1);
      }

      char buf[2]={'a',0};
      for (int i=0;i<spcnt;i++)          
      {
        getspecies(buf,0);        
        buf[0]++;
      }          
      break;      
    }

    // Loss weight
    case 'L':      
      if (sscanf(optarg, "%lf", &weightloss) != 1) {
        cerr << "Number expected in -L" << endl;
        exit(-1);
      }
      break;

    // Duplication weight
    case 'D':
      
      if (sscanf(optarg, "%lf", &weightdup) != 1) {
        cerr << "Number expected in -D" << endl;
        exit(-1);
      }
      break;

    // Set cost function 
    case 'C':
    {      
      
      costfunc = -1;
      if (!strcmp(optarg, "D")) costfunc = COSTDUPLICATION;
      if (!strcmp(optarg, "L")) costfunc = COSTLOSS;
      if (!strcmp(optarg, "DC")) costfunc = COSTDEEPCOAL;
      if (!strcmp(optarg, "DL")) costfunc = COSTDUPLICATIONLOSS;      
      if (!strcmp(optarg, "RF")) costfunc = COSTROBINSONFOULDS;
      
      if (costfunc < 0)
      {
        cerr << "Unknown cost function" << endl;
        exit(-1);
      }
      break;
    }

    default:
      cerr << "Unknown option: " << ((char) opt) << endl;
      exit(-1);
    }

  } // for



  for (size_t i = 0; i < sstvec.size(); i++)
    stvec.push_back(new RootedTree(sstvec[i]));

  for (size_t i = 0; i < sgtvec.size(); i++)
    gtvec.push_back(new RootedTree(sgtvec[i]));

      
  // verb_smp << "Number of gene trees: " << sgtvec.size() << endl;

  // int toosmall=0;

  // if (toosmall)
  //   verb_smp << "Warning: " << toosmall << " gene tree(s) of size 1 or 2 - ignored" << endl;

  if (OPT_PRINTINFO)
    {
      for (size_t i=0; i<specnames2id.size(); i++)
	       cout << i << " " << specnames[i] << endl;
    }
  
  VecRootedTree::iterator stpos;
  VecRootedTree::iterator gtpos;

  if (OPT_PRINTSPECIES)
    for (stpos = stvec.begin(); stpos != stvec.end(); ++stpos)
    { 
      (*stpos)->print(cout) << " " << "ispeciestree=" << (*stpos)->isspeciestree() << endl;
    } 

  if (OPT_PRINTGENE)
    for (gtpos = gtvec.begin(); gtpos != gtvec.end(); ++gtpos)    
      (*gtpos)->print(cout) << endl;

  if (OPT_PRINTDETAILED)    
  {

    cout << "Species trees:" << endl;
    for (stpos = stvec.begin(); stpos != stvec.end(); ++stpos)    
      (*stpos)->printdeb(cout,2) << endl;
    cout << "Gene trees:" << endl;
    for (gtpos = gtvec.begin(); gtpos != gtvec.end(); ++gtpos)    
      (*gtpos)->printdeb(cout,2) << endl;     
  }


  if (OPT_PRINTCOST)
  {
    
    for (stpos = stvec.begin(); stpos != stvec.end(); ++stpos)
    {
      
      for (gtpos = gtvec.begin(); gtpos != gtvec.end(); ++gtpos)        
      {
        if (OPT_PRINTCOST==2)
        {
          (*stpos)->print(cout) << " ";           
          (*gtpos)->print(cout) << " ";
        }
        cout << (*gtpos)->cost(*(*stpos),costfunc)  << endl;
      }
    }

  }


  gtc = new TreeClusters();
  for (gtpos = gtvec.begin(); gtpos != gtvec.end(); ++gtpos)    
      gtc->addtree(*gtpos);
  
  // Print quasi-consensus trees
  if (OPT_PRINTQUASICONSTREES)
  {
    RootedTree *sr=NULL;
    if (OPT_PRESERVEROOT)
    {
      if (stvec.size()>0)
          sr=stvec[0];
      else 
      {
          cout << "Exactly one species tree should be defined for quasi-consensus trees with preserve-root option" << endl;
          exit(2);        
      }
    }
      
    if (!gtvec.size())
    {
      cout << "Empty collection of gene trees - quasi consensus not found" << endl;
      return -1;
    }
    RootedTree *r;
    for (int i = 0; i < quasiconsensuscnt; i++)
    {
      
      r = gtc->genrootedquasiconsensus(sr);
      if (!r)
      {
        cerr << "Cannot create initial species tree" << endl;
        exit(-1);
      }
      cout << *r << endl;
    }
    return 0;
  }


  if (OPT_PRINTRANDSPECIESTREES)
  {      

      for (int i = 0; i < randomtreescnt; i++)
      {
        RootedTree *r = randspeciestree();
        if (!r)
        {
          cerr << "Cannot create initial random species tree" << endl;
          exit(-1);
        }
        cout << *r << endl;
      }   
  }     

        

  // Clean    
  for (size_t i = 0; i < sgtvec.size(); i++) 
    free(sgtvec[i]);

  for (size_t i = 0; i < sstvec.size(); i++) 
    free(sstvec[i]);

  for (size_t i = 0; i < stvec.size(); i++) 
    delete stvec[i];

  delete gtc; // Clean TreeClusters

  for (gtpos = gtvec.begin(); gtpos != gtvec.end(); ++gtpos)    
      delete *gtpos;

  cleanspecies();

  if (projectname) free(projectname);
}

