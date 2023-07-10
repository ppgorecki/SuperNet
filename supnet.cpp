/************************************************************************
 SuperNetwork Inference - based partially on FastUrec v2.05 
(c) Copyright 2005-2022 by Pawel Gorecki
 Written by P.Gorecki.
 Permission is granted to copy and use this program provided no fee is
 charged for it and provided that this copyright notice is not removed.
 *************************************************************************/

const char *SUPNET="0.03";

#include <time.h> 
#include <sys/time.h>
 
#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>

#include <string.h>
#include <iostream>
#include <set>
#include <vector>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <queue>
using namespace std;

#include "tools.h"
#include "bb.h"
#include "clusters.h"
#include "rtree.h"
#include "network.h"
#include "contrnet.h"
#include "hillclimb.h"
#include "costs.h"
#include "bitcluster.h"
#include "treespace.h"

double weightloss = 1.0;  // Unused
double weightdup = 1.0;   // Unused

int gsid = GSFULL;
char *gsdelim = NULL;
int gspos = 0;
int verbosehccost = 1; // 0=quiet
int verbosealg = 3; // <4 quiet

bool print_repr_inodtnaive = false;

typedef vector<RootedTree*> VecRootedTree;
typedef vector<Network*> VecNetwork;

#define BUFSIZE 100000

// Read trees from a given file
void readtrees(char *fn, vector<char*> &stvec)
{
  FILE *f;
  if (!strcmp(fn,"-")) f=stdin;
  else f = fopen(fn, "r");
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


// Read trees/networks from a string; trees/networks are separated by semicolon
void insertstr(vector<char*> &v,const char *t)
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


Network *addrandreticulations(int reticulationcnt_R, Network *n, int networktype)
{      
    for (int i=0; i<reticulationcnt_R; i++)
    {
      Network *prev = n;
      n = n->addrandreticulation("",networktype, uniform);
      if (!n)
      {
        cerr << "Cannot insert random " << (i+1) << "-th reticulation into " << *prev << endl;
        exit(-1);
      }        
      delete prev;      
    }
    return n;
}


Network *randnetwork(int reticulationcnt_R, int networktype)
{
    string r = randspeciestreestr();
    if (!r.length())
    {
      cerr << "Cannot create initial random species tree" << endl;
      exit(-1);
    }        

    return addrandreticulations(reticulationcnt_R,new Network(r),networktype);        

}

Network *randquasiconsnetwork(int reticulationcnt_R, int networktype, TreeClusters *gtc, RootedTree *preserverootst)
{
    string r = gtc->genrootedquasiconsensus(preserverootst);
    if (!r.length())
    {
      cerr << "Cannot create initial quasi consensus species tree" << endl;
      exit(-1);
    }      


    return addrandreticulations(reticulationcnt_R, new Network(r), networktype);

}



// interator over networks
// i should be initialized with -1
Network* netiterator(long int &i, VecNetwork &netvec, int &randomnetworkscnt, int &quasiconsensuscnt, 
  TreeClusters *gtc,
  RootedTree *preserverootst,

  int reticulationcnt_R, int networktype)

{
  
  i++;

  if (netvec.size()>0 && i<netvec.size())  
      return netvec[i];    
  
  
  if (randomnetworkscnt!=0)   // with -1 infitite 
  { 
    if (randomnetworkscnt>0)  randomnetworkscnt--;             

    return randnetwork(reticulationcnt_R, networktype);

  }


  if (quasiconsensuscnt!=0)  // with -1 infitite 
  { 
      if (quasiconsensuscnt>0) quasiconsensuscnt--;              
      return randquasiconsnetwork(reticulationcnt_R, networktype, gtc, preserverootst);                     
  }

  return NULL;

}

TreeSpace *globaltreespace;

int main(int argc, char **argv) 
{
  
  int opt;

  if (argc < 2) usage(argc, argv);

  VecRootedTree stvec, gtvec;
  VecNetwork netvec;
  
  struct timeval time; 
  gettimeofday(&time,NULL);
     
  srand((time.tv_sec * 1000) + (time.tv_usec / 1000));

  initbitmask(); // network
  
  CostFun *costfun = NULL;
  
  int OPT_QUASICONSTREES = 0;
  int OPT_RANDNETWORKS = 0;

  int OPT_PRINTGENE = 0;
  int OPT_PRINTNETWORK = 0;
  int OPT_PRINTSPECIES = 0;
  
  int OPT_PRINTCOST = 0;
  int OPT_PRINTDETAILED = 0;
  int OPT_PRINTINFO=0;  
  int OPT_DOT=0;  
  int OPT_PRINTDISPLAYTREES = 0;
  int OPT_COMPAREDAGS_BFTEST = 0;
  int OPT_COMPAREDAGS = 0;
  int OPT_DAGSHAPES = 0;
  int OPT_UNIQUEDAGS = 0;
  int OPT_PRINTSUBTREES = 0; // of s
  int OPT_UNIQUEDAGS_CNTS = 0;
  char *odt_option = NULL;

  int OPT_PRESERVEROOT=0;
  int OPT_EDITOPERATIONTEST = 0;
  int OPT_PRINTNODESTATS = 0;
  int OPT_DP = 0;
  int OPT_CONTRACTTEST = 0;
  int OPT_BB = 0;
  int OPT_BBSTATS = 0;

  int OPT_ODTNAIVE=0;

  int randomnetworkscnt = 0;
  int quasiconsensuscnt = 0;

  int reticulationcnt_R = 0;
  int networktype = NET_TREECHILD; // default network type
  int improvementthreshoold = 0;  

  int runnaiveleqrt_t = 13;  // default for DC, based on experiments

  bool odtlabelled = false;

  bool randnetuniform = false;

  string odtfile = "odt.log";
  string datfile = "odt.dat";

  const char* optstring = "e:g:s:G:S:N:l:q:L:D:C:r:A:n:do:O:R:K:t:b:z:E:aT:v:m:c:";
  vector<char*> sgtvec, sstvec, snetvec;

  COSTT bbstartscore = 0;
  bool bbstartscoredefined = false;
  float odtnaivesampling = 0.0;

  bool treereprtesting = 0;

  int maximprovements = 0;

  int maxdisplaytreecachesize = 1000000;
  
  while ((opt = getopt(argc, argv, optstring))   != -1)
  {    
    switch (opt) {

      case 'b':    
        bbstartscore = atof(optarg);
        bbstartscoredefined = true;
        break;

      case 'z':
        srand((unsigned int)atoi(optarg));
        break;

      case 'c':
        maxdisplaytreecachesize = atoi(optarg);
        break;

      case 'E':
        odtnaivesampling = atof(optarg);
        break;

      case 'a': 
        treereprtesting = 1;
        break;

      case 'e':        
        if (strchr(optarg,'i')) OPT_PRINTINFO=1;
        if (strchr(optarg,'a')) odtlabelled=true;

        if (strchr(optarg,'g')) OPT_PRINTGENE = 1;
        if (strchr(optarg,'t')) OPT_PRINTDISPLAYTREES = 1;
        if (strchr(optarg,'T')) OPT_PRINTDISPLAYTREES = 2;
        if (strchr(optarg,'n')) OPT_PRINTNETWORK = 1;
        if (strchr(optarg,'o')) OPT_ODTNAIVE = 1;

        if (strchr(optarg,'s')) OPT_PRINTSPECIES = 1;
        if (strchr(optarg,'S')) OPT_PRINTSPECIES = 2; // via printrepr
        if (strchr(optarg,'_')) OPT_PRINTSUBTREES = 1;
        if (strchr(optarg,'r')) OPT_PRESERVEROOT = 1;
        if (strchr(optarg,'c')) OPT_PRINTCOST = 1; // just numbers
        if (strchr(optarg,'C')) OPT_PRINTCOST = 2; // with trees
        if (strchr(optarg,'D')) OPT_PRINTDETAILED = 1;        
        if (strchr(optarg,'l')) OPT_PRINTNODESTATS = 1; // visible leaves stats
        if (strchr(optarg,'L')) OPT_PRINTNODESTATS = 2; // visible nodes stats

        if (strchr(optarg,'p')) OPT_COMPAREDAGS = 1;  // allvsall
        if (strchr(optarg,'x')) OPT_DAGSHAPES = 1;  
        if (strchr(optarg,'u')) OPT_UNIQUEDAGS = 1;   
        if (strchr(optarg,'U')) OPT_UNIQUEDAGS_CNTS = 1;  // TODO: unified approach to rand test and generator


        if (strchr(optarg,'R')) print_repr_inodtnaive = 1; 

        if (strchr(optarg,'1')) networktype = NET_CLASS1RELAXED; 
        if (strchr(optarg,'2')) networktype = NET_GENERAL;         
        if (strchr(optarg,'d')) OPT_DP = 1;
        if (strchr(optarg,'b')) OPT_BB = 1;


        if (strchr(optarg,'j')) OPT_BBSTATS |= 1;
        if (strchr(optarg,'J')) OPT_BBSTATS |= 2;
        if (strchr(optarg,'k')) OPT_BBSTATS |= 4;

        if (strchr(optarg,'N')) OPT_EDITOPERATIONTEST = 1; // nni test
        if (strchr(optarg,'M')) OPT_EDITOPERATIONTEST = 2; // tailmove
      
        // 
        // one net -r1 
        // 3 species -A3
        // 5 
        // ./supnet -R5 -r1 -A3 -e1X 6 1 && dot -Tpdf contr.dot -o c.pdf && evince c.pdf
        if (strchr(optarg,'X')) OPT_CONTRACTTEST = 1;
        
        if (strchr(optarg,'.')) OPT_COMPAREDAGS_BFTEST = 1; // hidden
        break;
    
    // case 'v': 
    //   verbose=atoi(optarg);
    //   break;

    case 'g':    
      insertstr(sgtvec,optarg);                
      break;

    case 'R':    
      reticulationcnt_R = atoi(optarg);
      break;

    case 't':    
      runnaiveleqrt_t = atoi(optarg);
      break;

    case 'm':    
      maximprovements = atoi(optarg);
      break;
    
    case 's':
      insertstr(sstvec,optarg);                    
      break;

    case 'n':
      insertstr(snetvec,optarg);                    
      break;


    case 'G': 
    {      
      readtrees(optarg, sgtvec);      
      break;
    }

    case 'S': 
    {      
      readtrees(optarg, sstvec);      
      break;
    }

    case 'N': 
    {      
      readtrees(optarg, snetvec);      
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

    // Gen quasi-consensus trees
    case 'q':
      OPT_QUASICONSTREES = 1;
      if (sscanf(optarg, "%d", &quasiconsensuscnt) != 1) {
        cerr << "Number expected in -w" << endl;
        exit(-1);
      }
      
      break;

    // Gen random trees/networks
    case 'r':
      OPT_RANDNETWORKS = 1;      
      if (sscanf(optarg, "%d", &randomnetworkscnt) != 1) {
        cerr << "Number expected in -r" << endl;
        exit(-1);
      }
      
      break;

    case 'K':      
      if (sscanf(optarg, "%d", &improvementthreshoold) != 1) {
        cerr << "Number expected in -K" << endl;
        exit(-1);
      }
      break;

    // define letter based species a,b,c,...
    case 'A':
    {
      int spcnt=0;
      if (sscanf(optarg, "%d", &spcnt) != 1) {
        cerr << "Number expected in -A" << endl;
        exit(-1);
      }
      if (spcnt>'z'-'a'+1)
      {
        cerr << "Too many species in -A" << endl;
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

    case 'd':
        OPT_DOT = 1;
        break;

    case 'O':
        odtfile = optarg;
        break;

    case 'T':
        datfile = optarg;
        break;

    case 'o':        
        odt_option = strdup(optarg);
        break;

    case 'v':

      // ugly
      if (strchr(optarg,'0')) verbosehccost = verbosealg = 0;
      if (strchr(optarg,'1')) verbosehccost = 1;
      if (strchr(optarg,'2')) verbosehccost = 2;
      if (strchr(optarg,'3')) verbosehccost = 3;
      if (strchr(optarg,'4')) verbosealg = 4;
      if (strchr(optarg,'5')) verbosealg = 5;
      if (strchr(optarg,'6')) verbosealg = 6;
      break;

    // Set cost function 
    case 'C':
    {      
      
      if (!strcmp(optarg, "D")) costfun = new CFDuplication();
      if (!strcmp(optarg, "L")) costfun = new CFLoss();
      if (!strcmp(optarg, "DC")) costfun = new CFDeepCoalescence();
      if (!strcmp(optarg, "DCE")) costfun = new CFDeepCoalescenceEdge();
      if (!strcmp(optarg, "DL")) costfun = new CFDuplicationLoss();      
      if (!strcmp(optarg, "RF")) costfun = new CFRobinsonFoulds();
      
      if (!costfun)
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
  {
      RootedTree *s = new RootedTree(sstvec[i]);
      if (!s->bijectiveleaflabelling())
      {
        cerr << "Bijective leaf labelling expected in a species tree: " << *s << endl;
        exit(-1);
      }
      stvec.push_back(s);
  } 

  for (size_t i = 0; i < sgtvec.size(); i++)
  {
    RootedTree *gtree = new RootedTree(sgtvec[i]);
    gtree->setid(i);
    gtvec.push_back(gtree);
    
  }

  globaltreespace = new TreeSpace(gtvec,maxdisplaytreecachesize);


 
  for (size_t i = 0; i < snetvec.size(); i++)
  {
    Network *n = new Network(snetvec[i]);
    if (!n->bijectiveleaflabelling())
    {
      cerr << "Bijective leaf labelling expected in a network: " << *n << endl;
      exit(-1);
    }    


    netvec.push_back(addrandreticulations(reticulationcnt_R,n,networktype));

  }
    
  // Print species names    
  if (OPT_PRINTINFO)    
      for (size_t i=0; i<specnames2id.size(); i++)
	       cout << i << " " << specnames[i] << endl;
      
  VecRootedTree::iterator stpos, gtpos;
  VecNetwork::iterator ntpos;


  TreeClusters *gtc = NULL;

  // Prepare clusters 
  if (OPT_QUASICONSTREES || OPT_RANDNETWORKS)
  {
      gtc = new TreeClusters();    
      for (gtpos = gtvec.begin(); gtpos != gtvec.end(); ++gtpos)    
        gtc->addtree(*gtpos);
      // cout << *gtc;
  }

  // Gen quasi consensus trees and insert into netvec as networks 
  // Add reticulations if -R is set 
  RootedTree *preserverootst = NULL;
  if (OPT_QUASICONSTREES)
  {
  
    if (OPT_PRESERVEROOT)
    {
      if (stvec.size()>0) preserverootst=stvec[0];
      else 
      {
          cout << "Exactly one species tree should be defined for quasi-consensus trees with preserve-root option" << endl;
          exit(2);        
      }
    }
      
    if (!gtvec.size())
    {
      cout << "Gene trees are required to infer quasi consensus" << endl;
      return -1;
    }
    
    if (!odt_option) // odt generated separately
       for (int i = 0; i < quasiconsensuscnt; i++)      
          netvec.push_back(randquasiconsnetwork(reticulationcnt_R, networktype, gtc, preserverootst));
  }



  // Random networks generated on the fly
  if (OPT_UNIQUEDAGS || OPT_UNIQUEDAGS_CNTS)
  {
      DagSet dagset(OPT_UNIQUEDAGS_CNTS,OPT_DAGSHAPES);               
      // get next network 
      Network *n; 
      long int i = -1;

      while  ((n = netiterator(i, netvec, randomnetworkscnt, quasiconsensuscnt, gtc, preserverootst, reticulationcnt_R, networktype))!=NULL)              

         dagset.add(n);        
      
      cout << dagset;     
      cerr <<  "unique=" << dagset.size() << " all=" << netvec.size() << endl;
      exit(0); // ignore rest opt
  }

  // Just printing
  if (OPT_PRINTSPECIES)
  {
    for (stpos = stvec.begin(); stpos != stvec.end(); ++stpos)    
    {
        if (OPT_PRINTSPECIES==2)
        {
          (*stpos)->printrepr() << endl;      
        }
        else 
        {
          cout << **stpos << endl;      
        }
    }
  }



  if (OPT_PRINTSUBTREES)
  {
    for (stpos = stvec.begin(); stpos != stvec.end(); ++stpos)    
    {        
        (*stpos)->printsubtrees(cout);              
    }    
  }
     
  if (OPT_PRINTGENE)
  {
    for (gtpos = gtvec.begin(); gtpos != gtvec.end(); ++gtpos)     
    {     
      cout << **gtpos << endl;  
    }
  }


  if (OPT_PRINTCOST)
  {    
    for (stpos = stvec.begin(); stpos != stvec.end(); ++stpos)        
    {
      for (gtpos = gtvec.begin(); gtpos != gtvec.end(); ++gtpos)        
      {
     
        if (OPT_PRINTCOST==2)        
        {
          cout << **stpos << " " << **gtpos << " ";
        }

        cout << costfun->computegt(**gtpos, **stpos)  << endl;
      }
    }
  }   


  // Gen random trees and store in netvec
  // Add reticulations if -R is set 
  if (OPT_RANDNETWORKS && !odt_option)
  {      
      for (int i = 0; i < randomnetworkscnt; i++)      

        netvec.push_back(randnetwork(reticulationcnt,networktype,randnetuniform));              

  }     

  if (OPT_PRINTNETWORK)
  {
    for (ntpos = netvec.begin(); ntpos != netvec.end(); ++ntpos)      
    {
      //(*ntpos)->printdeb(cout,2) << endl;
      cout << **ntpos << endl;     
    }
  }
      

  if (OPT_PRINTDETAILED)    
  {

    if (stvec.size()) cout << "Species trees:" << endl;
    for (stpos = stvec.begin(); stpos != stvec.end(); ++stpos)    
      (*stpos)->printdeb(cout,2) << endl;

    if (gtvec.size()) cout << "Gene trees:" << endl;
    for (gtpos = gtvec.begin(); gtpos != gtvec.end(); ++gtpos)    
      (*gtpos)->printdeb(cout,2) << endl;     

    if (netvec.size()) cout << "Networks:" << endl;
    for (ntpos = netvec.begin(); ntpos != netvec.end(); ++ntpos)    
      (*ntpos)->printdeb(cout,2) << endl;     
  }

  if (OPT_PRINTNODESTATS)    
  {    
    for (ntpos = netvec.begin(); ntpos != netvec.end(); ++ntpos)    
      (*ntpos)->visibilenodestats(OPT_PRINTNODESTATS,cout) << endl;     
  }

  // Generate dot file with networks and trees
  if (OPT_DOT)
  {
    ostream &s = cout; // todo files
    int dagnum = 0;
    s << "digraph SN {" << endl;


    for (stpos = stvec.begin(); stpos != stvec.end(); ++stpos)    
      (*stpos)->printdot(s,dagnum++) << endl;
  
    for (gtpos = gtvec.begin(); gtpos != gtvec.end(); ++gtpos)    
      (*gtpos)->printdot(s,dagnum++) << endl;     

    for (ntpos = netvec.begin(); ntpos != netvec.end(); ++ntpos)    
      (*ntpos)->printdot(s,dagnum++) << endl;     

    s << "}" << endl;
  }


  // Cost prepration 
  if (!costfun)
      costfun = new CFDeepCoalescence(); // default

  // Compute ODT cost by naive enumeration of display trees
  if (OPT_ODTNAIVE)
  {
    DISPLAYTREEID optid;
    ODTStats odtstats;
    for (ntpos = netvec.begin(); ntpos != netvec.end(); ++ntpos)            
      {
        if (OPT_ODTNAIVE==2)        
        {
          cout << **ntpos << " ";        
        }
        cout << (*ntpos)->odtcostnaive(gtvec, *costfun, odtstats, odtnaivesampling) << endl;        
    }
  }

  // Run DP algorithm to compute approx DCE 
  if (OPT_DP)
  {    
    for (ntpos = netvec.begin(); ntpos != netvec.end(); ++ntpos)            
    {
      for (gtpos = gtvec.begin(); gtpos != gtvec.end(); ++gtpos)          
      {
          cout << ((*ntpos)->approxmindce(**gtpos, *costfun)) << endl; 
      }
    }
    exit(0);
  }

  // Run BB algorithm to compute DCE
  if (OPT_BB)
  {
      BBTreeStats bbtreestats;
      ODTStats odtstats;

      for (ntpos = netvec.begin(); ntpos != netvec.end(); ++ntpos)            
        for (gtpos = gtvec.begin(); gtpos != gtvec.end(); ++gtpos)                        
        {
          double tm = gettime();
          COSTT dce = (*ntpos)->mindce(**gtpos, runnaiveleqrt_t, *costfun, odtstats, &bbtreestats, 
            bbstartscore, bbstartscoredefined);               
          cout << dce - (**gtpos).sizelf()*2 - 2;
          if (OPT_BBSTATS&4) 
            cout   << " " << dce 
                   << " " << (gettime() - tm) 
                   << " " << bbtreestats.minrtnumber 
                   << " " << bbtreestats.stats.naivecnt
                   << " " << bbtreestats.stats.naivetime
                   << " " << bbtreestats.stats.dpcnt
                   << " " << bbtreestats.stats.dptime;

          cout << endl;
        }

      if (OPT_BBSTATS&1) bbtreestats.savedot();
      if (OPT_BBSTATS&2) bbtreestats.savetsv();
      exit(0);   
  }
  
  // ODT heuristic using HC, BB and DP
  if (odt_option)
  {
    DISPLAYTREEID optid;    
    bool usenaive_oe = false;     

    if (strchr(odt_option,'q')) 
    {   
        // do not generate odt/dat files
        odtfile = ""; 
        datfile = ""; 
    }

    if (strchr(odt_option,'e')) 
    { 
        // use exhaustive enumeration
        usenaive_oe = true; 
    }

    int printstats = 0;
    if (strchr(odt_option,'s')) printstats=1;
    if (strchr(odt_option,'S')) printstats=2;

    HillClimb hc(gtvec); 
    
    EditOp *op; 
    if (strchr(odt_option,'N')) op = new NNI();
    else op = new TailMove(networktype);        

    NetworkHCStats *globalstats = new NetworkHCStats(networktype);
    
    long int i = -1;
    int lastimprovement=0;
    if (verbosealg>=4) 
        {
            cout << "HC" << " start:"               
              << " usenaive_oe=" << usenaive_oe
              << " runnaiveleqrt_t=" << runnaiveleqrt_t
              << endl;
        }

    while (1)    
    {

        // stopping criterion
        if (improvementthreshoold && (i-lastimprovement)>improvementthreshoold) break; // stop


        // get next initial network 
        Network *n = netiterator(i, netvec, randomnetworkscnt, quasiconsensuscnt, gtc, preserverootst, reticulationcnt_R, networktype);

        if (!n) break;        

        NetworkHCStats nhcstats(networktype);
        nhcstats.start();        

        // climb
        double cost = hc.climb(*op, n, *costfun, nhcstats, usenaive_oe, runnaiveleqrt_t, maximprovements);        

        nhcstats.finalize();
        
        if (globalstats->merge(nhcstats, printstats)) 
        {
            lastimprovement=i;
        }
          
        delete n; 

                  
    }

    globalstats->print(true);      
    delete op;


    if (odtfile.length()) 
    {

      // save odt file
      globalstats->save(odtfile);      
      
      // save dat file
      globalstats->savedat(datfile, odtlabelled); 
     

      if (verbosealg>=4)
      {
        cout << "Optimal networks saved: " << odtfile << endl;  
        cout << "Stats data save to: " << datfile << endl;
      }
    }



  }


  if (OPT_PRINTDISPLAYTREES)
  {
    for (ntpos = netvec.begin(); ntpos != netvec.end(); ++ntpos)    
    {
        DISPLAYTREEID tid = 0;
        Network *n = *ntpos;
        RootedTree *t = NULL;
        // cout << *n << endl;     
        while ((t=n->gendisplaytree(tid,t))!=NULL)   
        {
          if (OPT_PRINTDISPLAYTREES==2)
            cout << tid << " ";
          t->printrepr() << endl;       
          tid++;      
        }
    }
  }

  if (treereprtesting)
  {

    for (ntpos = netvec.begin(); ntpos != netvec.end(); ++ntpos)    
    {
        DISPLAYTREEID tid = 0;
        Network *n = *ntpos;
        SNode *t = NULL;
        cout << *n << endl;     
        while ((t=n->gendisplaytree2(tid, t, globaltreespace))!=NULL)   
        {
          //cout << tid << " <<";          
          ppSNode(cout, t) << endl;
          tid++;      
        }
    }

  }

  if (OPT_COMPAREDAGS)
  {
      int cnt=0, cntall=0;
      for (int i=0; i<netvec.size(); i++)    
      {
        // cout << *(netvec[i]) << endl;
        Network *n1 = netvec[i];
        cout << endl;
        for (int j=i+1; j<netvec.size(); j++)    
        { 
            Network *n2 = netvec[j];            
            bool e1 = n1->eqdags(n2,!OPT_DAGSHAPES);                        
            cout << *(n1) << "\t" << *(n2) << "\t" << e1 << endl;                        
        }
      }
      //cout <<  "eqcnt=" << cnt << " all=" << cntall << endl;
  }

  // ----------------- DEBUGs -----------------------

  if (OPT_EDITOPERATIONTEST==1)
  {
    for (ntpos = netvec.begin(); ntpos != netvec.end(); ++ntpos)            
      {
        std::ofstream dotf;
        dotf.open ("edit.dot", std::ofstream::out );

        int cnt=0;
        dotf << "digraph SN {" << endl;

        NNI nni;
        nni.init(*ntpos);

        cout << **ntpos << " verify==" << (*ntpos)->verifychildparent() << endl;
        (*ntpos)->printdeb(cout,2);
        (*ntpos)->printdot(dotf,cnt++);
        
        while (nni.next())
        {
          cout << "========" << endl << endl;
          int err = (*ntpos)->verifychildparent();

          cout << **ntpos << " verify==" << err << endl;          
          (*ntpos)->printdeb(cout,2); 
          (*ntpos)->printdot(dotf,cnt++);         

          if (err) break;
        }  

        dotf << "}" << endl;
        dotf.close();
        
      }
  }


  if (OPT_EDITOPERATIONTEST==2)
  {
    for (ntpos = netvec.begin(); ntpos != netvec.end(); ++ntpos)            
      {
        std::ofstream dotf;

        dotf.open ("tm.dot", std::ofstream::out );

        int cnt=0;
        dotf << "digraph SN {" << endl;

        TailMove tailmove(0);
        tailmove.init(*ntpos);

        cout << **ntpos << " verify==" << (*ntpos)->verifychildparent() << endl;
        // (*ntpos)->printdeb(cout,2);
        (*ntpos)->printdot(dotf,cnt++);
        
        while (tailmove.next())
        {
          cout << "========" << endl << endl;
          int err = (*ntpos)->verifychildparent();

          cout << **ntpos << " verify==" << err << endl;          
          // (*ntpos)->printdeb(cout,2); 
          (*ntpos)->printdot(dotf,cnt++);         

          if (err) break;
        }  

        dotf << "}" << endl;
        dotf.close();
        
      }
  }

 

  // Debug on
 if (OPT_COMPAREDAGS_BFTEST)
  {

      int cnt=0;
      while (1)
      {
        cnt++;
        string r1 = randspeciestreestr();
        string r2 = randspeciestreestr();
        int EE=0;
        int e1,e2,e3;
        if (!r1.length() || !r2.length())
        {
           cerr << "Cannot create initial random species tree" << endl;
           exit(-1);
        }        
        Network *n1 = addrandreticulations(reticulationcnt_R,new Network(r1),networktype);
        Network *n2 = addrandreticulations(reticulationcnt_R,new Network(r2),networktype);

        e1 = n1->eqdags(n2);                              
        e2 = n1->eqdagsbypermutations(n2);            

        cout << (*n1) << "\t" << (*n2) << "\tE1=" << e1 << "\tE2=" << e2 << "\tError=" << (e1!=e2) <<  endl;        

        cerr << (*n1) << "\t" << (*n2) << "\tE1=" << e1 << "\tE2=" << e2 << "\tError=" << (e1!=e2) <<  endl;        

        delete n1;
        delete n2;

        if (e1!=e2) break; // stop
    }
    exit(0);        
  }

  

  if (OPT_CONTRACTTEST)
  {
      // Last two args: leftretusage rightretusage
      // ./supnet -g "(a,(b,c))" -n "((b)#1,(((#2)#3,((#1)#2,(c,#3))),a))" -eXng 2 1
      // ./dot -Tpdf contr.dot -o c.pdf && evince c.pdf


      RETUSAGE retusage;
      emptyretusage(retusage);            

      //ugly
      if (optind+1<=argc)      
      {
        
        long us = atol(argv[optind]);
        
        int rid=0;
        while (us)
        {
          cout << "A" << us << " " << rid << endl;
          if (us&1) addleftretusage(retusage,rid);
          us=us>>1;
          rid++;
        }
      }


      if (optind+2<=argc)
      {        
        long us = atol(argv[optind+1]);
        int rid=0;
        while (us)
        {
          cout << "B" << us << " " << rid << endl;
          if (us&1) addrightretusage(retusage,rid);
          us=us>>1;
          rid++;
        }
      }

      cout << "R" << retusage << endl;
      cout << "Conflicted " << conflicted(retusage)  << endl;;

      for (int i=0; i<netvec.size(); i++)    
      {       

        Network *n1 = netvec[i];        
        std::stringstream ss;
        n1->print(ss);
        cout << ss.str() << endl;
        
        ContractedNetwork *c=new ContractedNetwork(ss.str());

        c->contract(retusage);        

        cout << "Rtcount " << c->rtcount() << " " << "rt=" << n1->rtcount() << endl;;

        for (gtpos = gtvec.begin(); gtpos != gtvec.end(); ++gtpos)          
          cout << "retmindc:" << (c->approxmindce(**gtpos, *costfun)) << endl;

        ofstream s("contr.dot");
        std::ofstream sf;
        sf.open ("contr.dot", std::ofstream::out );
        sf << "digraph SN {" << endl;
        sf << " inp [label=\"InRT=" <<  retusage << "\"]" <<endl;
        c->gendot(sf);
        c->gendotcontracted(sf);
        sf << "}" << endl;
        sf.close();        
        cout << c->newickrepr() << endl;


        RootedTree *t = NULL;
        DISPLAYTREEID tid = 0; // id of display tree
            
        while ((t=c->gendisplaytree(tid,t))!=NULL)       
        {  
          cout << tid << " " << (*t) << endl;
          tid++;
        }
      }
  }

  // clustergraphtester();
  
  // Cleaning


  delete globaltreespace;

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

}

