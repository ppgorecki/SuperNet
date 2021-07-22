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

 
#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>

#include <string.h>
#include <iostream>
#include <set>
#include <vector>
#include <algorithm>
#include <fstream>
using namespace std;

#include "tools.h"
#include "clusters.h"
#include "rtree.h"
#include "network.h"
#include "hillclimb.h"

#include <sstream>
#include <queue>

#define VERBOSE_QUIET 1
#define VERBOSE_CNTCOST 2
#define VERBOSE_SMP 3
#define VERBOSE_DETAILED 4

double weightloss = 1.0;  // Unused
double weightdup = 1.0;   // Unused

int gsid = GSFULL;
char *gsdelim = NULL;
int gspos = 0;

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
  
  int costfunc = COSTDUPLICATIONLOSS;
  
  int OPT_PRINTQUASICONSTREES = 0;
  int OPT_PRINTRANDSPECIESTREES = 0;
  int OPT_PRINTGENE = 0;
  int OPT_PRINTNETWORK = 0;
  int OPT_PRINTSPECIES = 0;
  int OPT_PRINTCOST = 0;
  int OPT_PRINTDETAILED = 0;
  int OPT_PRINTINFO=0;  
  int OPT_DOT=0;  
  int OPT_PRINTDISPLAYTREES = 0;
  char *odt = NULL;

  int OPT_PRESERVEROOT=0;
  int OPT_EDITOPERATIONTEST = 0;

  int OPT_ODTNAIVE=0;

  int randomtreescnt = 0;
  int quasiconsensuscnt = 0;

   
  const char* optstring = "e:g:s:G:S:N:l:w:L:D:C:r:A:n:do:";
  vector<char*> sgtvec, sstvec, snetvec;
  
  while ((opt = getopt(argc, argv, optstring))   != -1)
  {    
    switch (opt) {

      case 'e':        
        if (strchr(optarg,'i')) OPT_PRINTINFO=1;
        if (strchr(optarg,'g')) OPT_PRINTGENE = 1;
        if (strchr(optarg,'t')) OPT_PRINTDISPLAYTREES = 1;
        if (strchr(optarg,'T')) OPT_PRINTDISPLAYTREES = 2;
        if (strchr(optarg,'n')) OPT_PRINTNETWORK = 1;
        if (strchr(optarg,'o')) OPT_ODTNAIVE = 1;

        if (strchr(optarg,'s')) OPT_PRINTSPECIES = 1;
        if (strchr(optarg,'S')) OPT_PRINTSPECIES = 2; // via printrepr
        if (strchr(optarg,'r')) OPT_PRESERVEROOT = 1;
        if (strchr(optarg,'c')) OPT_PRINTCOST = 1; // just numbers
        if (strchr(optarg,'C')) OPT_PRINTCOST = 2; // with trees
        if (strchr(optarg,'D')) OPT_PRINTDETAILED = 1;

        if (strchr(optarg,'N')) OPT_EDITOPERATIONTEST = 1; // nni test
        if (strchr(optarg,'M')) OPT_EDITOPERATIONTEST = 2; // tailmove

        //if (strchr(optarg,'i')) OPT_USERSPTREEISSTARTING=0;        
        //if (strchr(optarg,'S')) OPT_SHOWSTARTINGTREES=1;
        break;
    
    // case 'v': 
    //   verbose=atoi(optarg);
    //   break;

    case 'g':    
      insertstr(sgtvec,optarg);                
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

    case 'o':        
        odt = strdup(optarg);
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
    gtvec.push_back(new RootedTree(sgtvec[i]));
 
  for (size_t i = 0; i < snetvec.size(); i++)
  {
    Network *n = new Network(snetvec[i]);
    if (!n->bijectiveleaflabelling())
    {
      cerr << "Bijective leaf labelling expected in a network: " << *n << endl;
      exit(-1);
    }
    netvec.push_back(n);
  }
    
  // Print species names    
  if (OPT_PRINTINFO)    
      for (size_t i=0; i<specnames2id.size(); i++)
	       cout << i << " " << specnames[i] << endl;
      
  VecRootedTree::iterator stpos, gtpos;
  VecNetwork::iterator ntpos;

  if (OPT_PRINTSPECIES)
    for (stpos = stvec.begin(); stpos != stvec.end(); ++stpos)    
    {
        if (OPT_PRINTSPECIES==2)
          (*stpos)->printrepr(cout) << endl;      
        else cout << **stpos << endl;      
    }
     
  if (OPT_PRINTGENE)
    for (gtpos = gtvec.begin(); gtpos != gtvec.end(); ++gtpos)          
      cout << **gtpos << endl;     

  if (OPT_PRINTNETWORK)
    for (ntpos = netvec.begin(); ntpos != netvec.end(); ++ntpos)      
      cout << **ntpos << endl;     
      

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

  if (OPT_PRINTCOST)
  {    
    for (stpos = stvec.begin(); stpos != stvec.end(); ++stpos)        
      for (gtpos = gtvec.begin(); gtpos != gtvec.end(); ++gtpos)        
      {
        if (OPT_PRINTCOST==2)        
          cout << **stpos << " " << **gtpos << " ";
        cout << (*gtpos)->cost(*(*stpos),costfunc)  << endl;
      }
  }


  if (OPT_ODTNAIVE)
  {
    DISPLAYTREEID optid;
    for (ntpos = netvec.begin(); ntpos != netvec.end(); ++ntpos)            
      {
        if (OPT_ODTNAIVE==2)        
          cout << **ntpos << " ";
        cout << (*ntpos)->odtnaivecost(gtvec,costfunc,optid) << " optid=" << optid << endl;
        
    }
  }


  // Run hill climbing
  if (odt)
  {
    DISPLAYTREEID optid;
    int verbose = 0;
    if (strchr(odt,'1')) verbose = 1;
    if (strchr(odt,'2')) verbose = 2;
    if (strchr(odt,'3')) verbose = 3;

    HillClimb hc(gtvec,verbose,strchr(odt,'s'));    

    for (ntpos = netvec.begin(); ntpos != netvec.end(); ++ntpos)            
      {
        EditOp *op; 

        if (strchr(odt,'N')) op = new NNI();
        else op = new TailMove(strchr(odt,'t'));        
        
        // climb
        double cost = hc.climb(*op,*ntpos, costfunc);
        
        // Print network:
        if (strchr(odt,'n')) 
          cout << **ntpos;

        // Optimal cost 
        cout << " " << cost << endl;        

        delete op;
    }
  }


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
          cout << *t << endl;       
          tid++;      
        }
    }
  }

 
  TreeClusters *gtc = NULL;

  // Prepare clusters 
  if (OPT_PRINTQUASICONSTREES || OPT_PRINTRANDSPECIESTREES)
  {
      gtc = new TreeClusters();    
      for (gtpos = gtvec.begin(); gtpos != gtvec.end(); ++gtpos)    
        gtc->addtree(*gtpos);
  }

  if (OPT_PRINTQUASICONSTREES)
  {
    RootedTree *sr=NULL;
    if (OPT_PRESERVEROOT)
    {
      if (stvec.size()>0) sr=stvec[0];
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

    RootedTree *r;
    for (int i = 0; i < quasiconsensuscnt; i++)
    {
      
      r = gtc->genrootedquasiconsensus(sr);
      if (!r)
      {
        cerr << "Cannot create initial species tree" << endl;
        exit(-1);
      }
      cout << *r;
    }    
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

}

