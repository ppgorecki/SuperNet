
/************************************************************************
 SuperNetwork Inference - based partially on FastUrec v2.05
(c) Copyright 2005-2023 by Pawel Gorecki
 Written by P.Gorecki.

 Permission is granted to copy and use this program provided no fee is
 charged for it and provided that this copyright notice is not removed.
 *************************************************************************/

const char *SUPNET = "0.10";

#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <time.h>
#include <unistd.h>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <queue>
#include <set>
#include <sstream>
#include <vector>
using namespace std;

#include <bits/stdc++.h>

#include "bb.h"
#include "bitcluster.h"
#include "clusters.h"
#include "contrnet.h"
#include "costs.h"
#include "hillclimb.h"
#include "network.h"
#include "randnets.h"
#include "rtree.h"
#include "tools.h"
#include "topsort.h"
#include "treespace.h"
using namespace std;

double weightloss = 1.0; // Unused
double weightdup = 1.0;  // Unused

int gsid = GSFULL;
char *gsdelim = NULL;
int gspos = 0;
int verbosehccost = 1; // 0=quiet
int verbosealg = 3;    // <4 quiet

string opt_outfiles = ODTBASENAME;
string opt_outdirectory = "";

int flag_globaldagcache = 0; 
int flag_saveextnewick = 0; 
int flag_autooutfiles = 0; 
int flag_hcsavewhenimproved = 0;
int flag_hcsamplerstats = 0;
int flag_hcsamplingmaxnetstonextlevel = 0; // unlimited
int flag_hcdetailedsummarydat = 0;
float opt_hcstoptime = 0; 
bool flag_hcignorecostgeq = 0; // if set 
float opt_hcignorecostgeq = 0;
int flag_hcsavefinalodt = 0; // ignore hcignorecostgeq for final files

int print_repr_inodtnaive = 0;

typedef vector<RootedTree *> VecRootedTree;
typedef vector<Network *> VecNetwork;

string flag_retidprefix = "";

#define BUFSIZE 100000

// Read trees from a given file
void readtrees(char *fn, vector<char *> &stvec) {
  FILE *f;
  if (!strcmp(fn, "-"))
    f = stdin;
  else
    f = fopen(fn, "r");
  if (!f) {
    cerr << "Cannot open file " << fn << endl;
    exit(-1);
  }
  char buf[BUFSIZE];
  while (1) {
    if (!fgets(buf, BUFSIZE, f))
      break;
    stvec.push_back(strdup(buf));
  }
  fclose(f);
}

// Read trees/networks from a string; trees/networks are separated by semicolon
void insertstr(vector<char *> &v, const char *t) {
  char *x, *y = strdup(t);
  x = y;

  while (strlen(x) > 0) {
    char *p = strchr(x, ';');
    if (p) {
      *p = 0;
      v.push_back(strdup(x));
      x = p + 1;
    } else {
      v.push_back(strdup(x));
      break;
    }
  }
  free(y);
}


TreeSpace *globaltreespace;

int main(int argc, char **argv) {
  int opt;

  if (argc < 2)
    usage(argc, argv);

  VecRootedTree stvec, gtvec;
  VecNetwork netvec;

  struct timeval time;
  gettimeofday(&time, NULL);

  unsigned int randseed = (time.tv_sec * 1000) + (time.tv_usec / 1000);
  
  initbitmask(); // network

  CostFun *costfun = NULL;

  int opt_quasiconsensus = 0;
  int opt_randnetworks = 0;

  int flag_print_gene_trees = 0;
  int flag_print_networks = 0;
  int flag_print_species_trees = 0;



  int flag_ptreecost = 0;
  int flag_extra_print_detailed = 0;
  int flag_extra_print_species_dictionary = 0;
  int flag_dot = 0;
  int flag_print_display_trees = 0;
  int flag_print_display_trees_with_ids = 0;
  int flag_testcomparedags = 0;
  int flag_comparedags = 0;
  int flag_dagshapes = 0;
  int flag_uniquedags = 0;
  int flag_printstsubtrees = 0; // of s
  int flag_uniquedagscnts = 0;

  int flag_preserveroot = 0;
  int flag_testeditnni = 0;

  int flag_extra_visible_node_stats = 0;
  int flag_extra_visible_leaf_stats = 0;
  int flag_dpalgorithm = 0;
  int flag_testcontract = 0;
  int flag_bbalgorithm = 0;
  int flag_bbalgorithmSTATS = 0;

  int flag_odt_naive_gtvsnet = 0;

  int randomnetworkscnt = 0;
  int quasiconsensuscnt = 0;

  int reticulationcnt_R = 0;
  int networkclass = NET_TREECHILD; // default network type

  int hcrunnaiveleqrt_t = 13; // default for DC, based on experiments

  int odtlabelled = 0;

  int randnetuniform = 0;

  
  
  vector<char *> sgtvec, sstvec, snetvec;

  COSTT bbstartscore = 0;
  int bbstartscoredefined = 0;
  string displaytreesampling = "";

  int maxdisplaytreecachesize = 1000000;

  int flag_hcusenaive = 0;
  

  int net_treechild = 1;
  int net_relaxed = 0;
  int net_general = 0;

  int flag_timeconsistent = 0;
  int flag_notimeconsistent = 0;
  int timeconsistency = 0;  // 0 - no restrictions, 1-timeconsistent, 2-notimeconsistent
  
  int flag_ptreecostext = 0;
  int flag_noodtfiles = 0;
  int flag_bbtreesearch = 0;
  int flag_bbtsvstats = 0;
  int flag_bbtimestats = 0;
  int flag_hcalgorithm = 0;
  int flag_runnaiveleqrt = 0;
  int flag_hcedit_nni = 0;
  int flag_detectclass = 0;
  int flag_testtreerepr = 0;
  int flag_hcrunstats = 0;
  int flag_hcrunstatsext = 0;
  int flag_hcrunstatsalways = 0;
  int flag_hcdetailedsummary = 0;
  int flag_match_after = 0;
  int flag_match_before = 0;
  int flag_match_atpos = 0;
  int flag_print_network_clusters = 0;
  int flag_cutwhendtimproved = 0;

  int testdisplaytreesampling = 0;

  int hcmaximprovements = 0;
  int hcstopinit = 0;
  int hcstopclimb = 0;
  int flag_odtcost = 0;
  int flag_bestneworks = 0;
  char *matchparam = NULL;

  char *opt_guideclusters = NULL;
  char *opt_guidetree = NULL;
  
  costfun = new CFDeepCoalescence();

  struct option longopts[] = {

      // -l, -r, -R, -A

      {"verbose", required_argument, NULL, 'v'},

      // OK -------------------------------------
      {"genetrees", required_argument, NULL, 'g'},
      {"networks", required_argument, NULL, 'n'},
      {"speciestrees", required_argument, NULL, 's'},
      {"genetreesfile", required_argument, NULL, 'G'},
      {"networksfile", required_argument, NULL, 'N'},
      {"speciestreesfile", required_argument, NULL, 'S'},
      {"noodtfiles", no_argument, &flag_noodtfiles, 1},
      {"dot", no_argument, &flag_dot, 1}, // former -d
      {"maxdisplaytreecachesize", required_argument, NULL, 'V'},
      {"randseed", required_argument, NULL, 'z'},      
      {"ptreesinodtnaive", no_argument, &print_repr_inodtnaive, 1},
      {"odtlabelled", no_argument, &odtlabelled, 1},
      {"pgenetrees", no_argument, &flag_print_gene_trees, 1},
      {"pspeciestrees", no_argument, &flag_print_species_trees, 1},
      {"pdisplaytrees", no_argument, &flag_print_display_trees, 1},
      {"pdisplaytreesext", no_argument, &flag_print_display_trees_with_ids, 1},
      {"pnetworks", no_argument, &flag_print_networks, 1},
      {"pnetworkclusters", no_argument, &flag_print_network_clusters, 1},
      {"ptreecost", no_argument, &flag_ptreecost, 1},
      {"ptreecostext", no_argument, &flag_ptreecostext, 1},
      {"comparedags", no_argument, &flag_comparedags, 1},       //'p'
      {"dagshapes", no_argument, &flag_dagshapes, 1},           //'x'
      {"uniquedags", no_argument, &flag_uniquedags, 1},         
      {"uniquedagscnts", no_argument, &flag_uniquedagscnts, 1}, //'U'

      {"preserveroot", no_argument, &flag_preserveroot, 1},
      {"treechild", no_argument, &net_treechild, 1},
      {"relaxed", no_argument, &net_relaxed, 1},
      {"general", no_argument, &net_general, 1},
      {"timeconsistent", no_argument, &flag_timeconsistent, 1},
      {"notimeconsistent", no_argument, &flag_notimeconsistent, 1},

      {"outfiles", required_argument, NULL, 'O'},
      {"outdirectory", required_argument, NULL, 'd'},

      {"detectclass", no_argument, &flag_detectclass, 1},
      {"uniformedgesampling", no_argument, &randnetuniform, 1},
      {"displaytreesampling", required_argument, NULL, 'U'},
      {"testdisplaytreesampling", required_argument, NULL, '2'},
      {"pspeciesdictionary", no_argument, &flag_extra_print_species_dictionary,
       1},
      {"pstsubtrees", no_argument, &flag_printstsubtrees, 1},
      {"reachablenodescnt", no_argument, &flag_extra_visible_node_stats, 1},
      {"reachableleavescnt", no_argument, &flag_extra_visible_leaf_stats, 1},
      {"pdetailed", no_argument, &flag_extra_print_detailed, 1},

      {"DP", no_argument, &flag_dpalgorithm, 1},

      {"odtnaivecost", no_argument, &flag_odt_naive_gtvsnet, 1},
      {"cost", required_argument, NULL, 'C'},
      {"odtcost", no_argument, &flag_odtcost, 1},

      {"BB", no_argument, &flag_bbalgorithm, 1},
      {"bbtreesearch", no_argument, &flag_bbtreesearch, 1},
      {"bbtsvstats", no_argument, &flag_bbtsvstats, 1},
      {"bbtimestats", no_argument, &flag_bbtimestats, 1},
      {"bbstartscore", required_argument, NULL, '1'},

      {"HC", optional_argument, &flag_hcalgorithm, 'o'},
      {"runnaiveleqrt", required_argument, &flag_runnaiveleqrt, 't'},
      {"hcusenaive", no_argument, &flag_hcusenaive, 1},
      {"hcnnimove", no_argument, &flag_hcedit_nni, 1},
      {"hcrunstats", no_argument, &flag_hcrunstats, 1},
      {"hcsamplerstats", no_argument, &flag_hcsamplerstats, 1},
      {"hcrunstatsext", no_argument, &flag_hcrunstatsext, 1},
      {"hcrunstatsalways", no_argument, &flag_hcrunstatsalways, 1},
      {"hcsavewhenimproved", no_argument, &flag_hcsavewhenimproved, 1},
      {"hcsamplingmaxnetstonextlevel", required_argument, NULL, '3'},
      {"hcdetailedsummary", no_argument, &flag_hcdetailedsummary, 1},
      {"hcdetailedsummarydat", no_argument, &flag_hcdetailedsummarydat, 1},
      {"autooutfiles", no_argument, &flag_autooutfiles, 1},
      {"cutwhendtimproved", no_argument, &flag_cutwhendtimproved, 1},

      {"bestnetworks", no_argument, &flag_bestneworks, 1},
      
      {"hcmaximprovements", required_argument, NULL, 'F'},
      {"hcignorecostgeq", required_argument, NULL, '4'},
      {"hcsavefinalodt", no_argument, &flag_hcsavefinalodt, 1},
      {"hcstopinit", required_argument, NULL, 'Y'},
      {"hcstopclimb", required_argument, NULL, 'Z'},
      {"hcstoptime", required_argument, NULL, 'P'},

      // TODO: implement rules
      {"matchafter", no_argument, &flag_match_after, 1},
      {"matchbefore", no_argument, &flag_match_before, 1},
      {"matchatpos", no_argument, &flag_match_atpos, 1},
      {"matchparam", required_argument, NULL, 'l'},

      // debug
      {"testtreerepr", no_argument, &flag_testtreerepr, 1},
      {"testcontract", no_argument, &flag_testcontract, 1},
      {"testcomparedags", no_argument, &flag_testcomparedags, 1},
      {"testeditnni", no_argument, &flag_testeditnni, 1},
      {"testedittailmove", no_argument, &flag_testeditnni, 2},

      {"guideclusters", required_argument, NULL, 'u'},
      {"guidetree", required_argument, NULL, 'e'},

      {"globaldagcache", no_argument, &flag_globaldagcache, 1 },

      {"saveextnewick", no_argument, &flag_saveextnewick, 1 },

      {"retidprefix", required_argument, NULL, 'Q'},

      {0, 0, 0, 0}

  };

  const char *optstring =
      "c:z:b:U:g:s:G:S:N:n:R:t:m:V:X:Y:Z:l:q:r:A:L:D:O:v:C:1:Wu:e:P:;";

  while ((opt = getopt_long(argc, argv, optstring, longopts, NULL)) != -1) {
    switch (opt) {
    case 'V':
      maxdisplaytreecachesize = atoi(optarg);
      break;

    case 'z':
      randseed = (unsigned int)atoi(optarg);      
      break;

    case '1':
      bbstartscore = atof(optarg);
      break;

    case 'U':
      displaytreesampling = string(optarg);
      break;

    case '3':
      flag_hcsamplingmaxnetstonextlevel=atoi(optarg);
      break;

    case '4': 
      opt_hcignorecostgeq = atof(optarg);
      flag_hcignorecostgeq = true;
      break;


    case '2':
      testdisplaytreesampling = atoi(optarg);
      break;

    case 'P':
      opt_hcstoptime = atof(optarg);
      break;

    case 'g':
      insertstr(sgtvec, optarg);
      break;

    case 'R':
      reticulationcnt_R = atoi(optarg);
      break;

    case 't':
      flag_runnaiveleqrt = atoi(optarg);
      break;

    case 'X': {
      hcmaximprovements = atoi(optarg);
      break;
    }

    case 'Y': {
      hcstopinit = atoi(optarg);
      break;
    }

    case 'Z': {
      hcstopclimb = atoi(optarg);
      break;
    }

    case 's':
      insertstr(sstvec, optarg);
      break;

    case 'n':
      insertstr(snetvec, optarg);
      break;

    case 'G': {
      readtrees(optarg, sgtvec);
      break;
    }

    case 'O': 
    {
      opt_outfiles = optarg;      
      break;
    }

    case 'd': 
    {
      opt_outdirectory = optarg;      
      break;
    }

    case 'S': {
      readtrees(optarg, sstvec);
      break;
    }

    case 'N': {
      readtrees(optarg, snetvec);
      break;
    }

    case 'u': 
    {
      opt_guideclusters = strdup(optarg);
      break;
    }

    case 'e': 
    {
      opt_guidetree = strdup(optarg);
      break;
    }

    case 'Q':
    {
      flag_retidprefix = string(optarg);
      break;
    }

    // Matching gene tree and species tree labels
    case 'l': {
      // TODO
      cerr << "Matching (-l) not implemented yet" << endl;
      matchparam = strdup(optarg);

      // switch (optarg[0])
      // {
      // case 'p':
      //   gspos = atoi(optarg + 1);
      //   gsid = GSPOS;
      //   break;
      // case 'a':
      //   gsid = GSAFTER;
      //   gsdelim = strdup(optarg + 1);
      //   break;
      // case 'b':
      //   gsid = GSBEFORE;
      //   gsdelim = strdup(optarg + 1);
      //   break;
      // default:
      //   cerr << "Invalid code for -l option expected aDELIM,bDELIM or
      //   p[-]number" << endl; exit(-1);
      // }
      break;
    }

    // Gen quasi-consensus trees
    case 'q':
      opt_quasiconsensus = 1;
      if (sscanf(optarg, "%d", &quasiconsensuscnt) != 1) {
        cerr << "Number expected in -w" << endl;
        exit(-1);
      }

      break;

    // Gen random trees/networks
    case 'r':
      opt_randnetworks = 1;
      if (sscanf(optarg, "%d", &randomnetworkscnt) != 1) {
        cerr << "Number expected in -r" << endl;
        exit(-1);
      }

      break;

    // define letter based species a,b,c,...
    case 'A': {
      int spcnt = 0;
      if (sscanf(optarg, "%d", &spcnt) != 1) {
        cerr << "Number expected in -A" << endl;
        exit(-1);
      }
      if (spcnt > 'z' - 'a' + 1) {
        cerr << "Too many species in -A" << endl;
        exit(-1);
      }

      char buf[2] = {'a', 0};
      for (int i = 0; i < spcnt; i++) {
        getspecies(buf, 0);
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

    case 'v':

      // ugly
      if (strchr(optarg, '0'))
        verbosehccost = verbosealg = 0;
      if (strchr(optarg, '1'))
        verbosehccost = 1;
      if (strchr(optarg, '2'))
        verbosehccost = 2;
      if (strchr(optarg, '3'))
        verbosehccost = 3;
      if (strchr(optarg, '4'))
        verbosealg = 4;
      if (strchr(optarg, '5'))
        verbosealg = 5;
      if (strchr(optarg, '6'))
        verbosealg = 6;
      break;

    // Set cost function
    case 'C': {
      if (!strcmp(optarg, "D"))
        costfun = new CFDuplication();
      if (!strcmp(optarg, "L"))
        costfun = new CFLoss();
      // if (!strcmp(optarg, "DC")) costfun = new CFDeepCoalescence();
      if (!strcmp(optarg, "DCE"))
        costfun = new CFDeepCoalescenceEdge();
      if (!strcmp(optarg, "DL"))
        costfun = new CFDuplicationLoss();
      if (!strcmp(optarg, "RF"))
        costfun = new CFRobinsonFoulds();

      if (!costfun) {
        cerr << "Unknown cost function" << endl;
        exit(-1);
      }
      break;
    }

    case 0: /* getflag_long() set a variable, just keep going */
      break;

    default:
      cerr << "Unknown option: " << ((char)opt) << endl;
      exit(-1);
    }

  } // for

  srand(randseed);

  // Set network classes
  if (net_general)
  {
    networkclass = NET_GENERAL;
  }
  if (net_relaxed)
  {
    networkclass = NET_CLASS1RELAXED;
  }

  // Odt files flags
  if (flag_noodtfiles) 
  {
    // do not generate odt/dat files
    opt_outfiles = "";    
  }

  // Apply time consistent flags
  if (flag_timeconsistent)
  {
    timeconsistency = TIMECONSISTENT;
  }
  else 
    if (flag_notimeconsistent) 
    {
      timeconsistency = NOTIMECONSISTENT;
    }

  // Parse species trees
  for (size_t i = 0; i < sstvec.size(); i++) 
  {
    RootedTree *s = new RootedTree(sstvec[i]);
    if (!s->bijectiveleaflabelling()) {
      cerr << "Bijective leaf labelling expected in a species tree: " << *s
           << endl;
      exit(-1);
    }
    stvec.push_back(s);
  }

  // Parse gene trees
  for (size_t i = 0; i < sgtvec.size(); i++) 
  {
    RootedTree *gtree = new RootedTree(sgtvec[i]);
    gtree->setid(i);
    gtvec.push_back(gtree);
  }

  globaltreespace = new TreeSpace(gtvec, maxdisplaytreecachesize);

  // Add random reticulations
  for (size_t i = 0; i < snetvec.size(); i++) 
  {
    Network *n = new Network(snetvec[i]);
    if (!n->bijectiveleaflabelling()) 
    {
      cerr << "Bijective leaf labelling expected in a network: " << *n << endl;
      exit(-1);
    }

    netvec.push_back(addrandreticulations(reticulationcnt_R, n, networkclass,
                                          timeconsistency, randnetuniform, NULL, NULL));
  }

  // Print species names
  if (flag_extra_print_species_dictionary)
    for (size_t i = 0; i < specnames2id.size(); i++)
      cout << i << " " << specnames[i] << endl;

  Clusters *genetreeclusters = NULL;

  // Prepare clusters
  if (opt_quasiconsensus || opt_randnetworks) 
  {
    genetreeclusters = new Clusters();
    for (auto & gtpos: gtvec) 
      genetreeclusters->adddag(gtpos);    
  }

  Clusters *guideclusters = NULL;
  if (opt_guideclusters)
  {
    guideclusters = new Clusters();
    guideclusters->addtree(opt_guideclusters);    
  }

  Clusters *guidetree = NULL;
  if (opt_guidetree)
  {
    guidetree = new Clusters();
    guidetree->addtree(opt_guidetree);    
  }

  // Gen quasi consensus trees and insert into netvec as networks
  // Add reticulations if -R is set
  RootedTree *preserverootst = NULL;
  if (opt_quasiconsensus) 
  {
    if (flag_preserveroot) 
    {
      if (stvec.size() > 0)
        preserverootst = stvec[0];
      else 
      {
        cout << "Exactly one species tree should be defined for "
                "quasi-consensus trees with preserve-root option"
             << endl;
        exit(2);
      }
    }

    if (!gtvec.size() && !guideclusters && !guidetree) 
    {
        cout << "Gene trees (-g/-G), a guide tree (--guidetree) or guide clusters (--guideclusters) are required to infer quasi consensus" << endl;
        return -1;
    }


    if (!flag_hcalgorithm) // odt generated separately
    {
      
      for (int i = 0; i < quasiconsensuscnt; i++)
        netvec.push_back(randquasiconsnetwork(reticulationcnt_R, networkclass,
                                              timeconsistency, genetreeclusters,
                                              preserverootst, 
                                              guideclusters,
                                              guidetree));
    }
  }

  // Random networks generated on the fly
  if (flag_uniquedags || flag_uniquedagscnts) 
  {
    DagSet dagset(flag_uniquedagscnts, flag_dagshapes);
    // get next network
    Network *n;
    long int i = -1;

    Dag *src;
    NetIterator netiterator(netvec, randomnetworkscnt, quasiconsensuscnt, genetreeclusters,
                          preserverootst, reticulationcnt_R, networkclass,
                          timeconsistency, randnetuniform, guideclusters, 
                           guidetree);
    while ((n = netiterator.next()) != NULL)
       dagset.add(n, &src);

    cout << dagset;
    cerr << "unique=" << dagset.size() << " all=" << netvec.size() << endl;
    exit(0); 
    // ignore rest opt
  }

  // Just printing
  if (flag_print_species_trees) 
  {
    for (auto & stpos : stvec) 
    {
      if (flag_print_species_trees == 2) 
      {
        stpos->printrepr() << endl;
      } else 
      {
        cout << *stpos << endl;
      }
    }
  }

  if (flag_printstsubtrees) 
  {
    for (auto &stpos : stvec) 
    {
      stpos->printsubtrees(cout);
    }
  }

  if (flag_print_gene_trees) 
  {
    for (auto gtpos : gtvec) 
    {
      cout << *gtpos << endl;
    }
  }


  if (flag_ptreecost || flag_ptreecostext) 
  {
    for (auto &stpos : stvec) 
    {
      for (auto &gtpos : gtvec) 
      {
        if (flag_ptreecostext) 
        {
          cout << *stpos << " " << *gtpos << " ";
        }

        cout << costfun->computegt(*gtpos, *stpos) << endl;
      }
    }
  }

  // Gen random trees and store in netvec
  // Add reticulations if -R is set
  if (opt_randnetworks && !flag_hcalgorithm) 
  {

    for (int i = 0; i < randomnetworkscnt; i++)
    {
      netvec.push_back(randnetwork(reticulationcnt_R, networkclass,
                                   timeconsistency, randnetuniform));
    }
  }

  if (flag_print_networks) 
  {
    for (auto &ntpos: netvec) 
    {
      //(*ntpos)->printdeb(cout,2) << endl;
      cout << *ntpos << endl;
      if (flag_print_network_clusters)
      {
        ntpos->printclusters(cout);
      }
    }
  }

  if (flag_extra_print_detailed) 
  {
    if (stvec.size())
    {
      cout << "Species trees:" << endl;
    }

    for (auto &stpos: stvec)
    {
      stpos->printdeb(cout, 2) << endl;
    }

    if (gtvec.size())
    {
      cout << "Gene trees:" << endl;
    }

    for (auto &gtpos: gtvec)
    {
      gtpos->printdeb(cout, 2) << endl;
    }

    if (netvec.size())
    {
      cout << "Networks:" << endl;
    }

    for (auto &ntpos: netvec)
    {
      ntpos->printdeb(cout, 2) << endl;
    }
  }

  if (flag_extra_visible_node_stats) 
  {
    for (auto &ntpos: netvec)
    {
      ntpos->visibilenodestats(1, cout) << endl;
    }
  }

  if (flag_extra_visible_leaf_stats) 
  {
    for (auto &ntpos: netvec)
      ntpos->visibilenodestats(2, cout) << endl;
  }

  // Generate dot file with networks and trees
  if (flag_dot) 
  {    
    ostream &s = cout; // todo files
    int dagnum = 0;
    s << "digraph SN {" << endl;

    for (auto &stpos: stvec)
    {
      stpos->printdot(s, dagnum++) << endl;
    }

    for (auto &gtpos: gtvec)
    {
      gtpos->printdot(s, dagnum++) << endl;
    }

    for (auto &ntpos: netvec)
    {
      ntpos->printdot(s, dagnum++) << endl;
    }

    s << "}" << endl;
  }

  // Cost prepration
  if (!costfun)
  {
    costfun = new CFDeepCoalescence(); // default
  }

  // Compute ODT cost by naive enumeration of display trees
  if (flag_odt_naive_gtvsnet) 
  {
    DISPLAYTREEID optid;
    ODTStats odtstats;
    for (auto &ntpos: netvec) 
    {
      if (flag_odt_naive_gtvsnet == 2) 
      {
        cout << *ntpos << " ";
      }
      cout << ntpos->odtcostnaive(gtvec, *costfun, odtstats, 0)
           << endl;
    }
  }

  // Compute ODT cost 
  if (flag_odtcost) 
  {
    DagSet visiteddags;
    NetworkHCStatsBase *stats =
        new NetworkHCStatsBase(networkclass, timeconsistency, visiteddags, randseed);
    for (auto &net: netvec) 
    {
      double cost = net->odtcost(gtvec, *costfun, flag_hcusenaive, hcrunnaiveleqrt_t, stats->getodtstats());
      cout << cost << " " << *net << endl;
    }
    exit(0);
  }

  // Run DP algorithm to compute approx DCE
  if (flag_dpalgorithm) 
  {
    for (auto &ntpos: netvec) 
    {
      for (auto &gtpos: gtvec) 
      {
        cout << ntpos->approxmindce(*gtpos, *costfun) << endl;
      }
    }
    exit(0);
  }

  // Run BB algorithm to compute DCE
  if (flag_bbalgorithm) 
  {
    BBTreeStats bbtreestats;
    ODTStats odtstats;

    for (auto &ntpos: netvec)
    {
      for (auto &gtpos: gtvec) 
      {
        double tm = gettime();
        COSTT dce =
            ntpos->mindce(*gtpos, flag_runnaiveleqrt, *costfun, odtstats,
                             &bbtreestats, bbstartscore, bbstartscoredefined);
        cout << dce - (*gtpos).sizelf() * 2 - 2;

        if (flag_bbtimestats)
          cout << " " << dce << " " << (gettime() - tm) << " "
               << bbtreestats.minrtnumber << " " << bbtreestats.stats.naivecnt
               << " " << bbtreestats.stats.naivetime << " "
               << bbtreestats.stats.dpcnt << " " << bbtreestats.stats.dptime;

        cout << endl;
      }
    }

    if (flag_bbtreesearch)
    {
      bbtreestats.savedot();
    }

    if (flag_bbtsvstats)
    {
      bbtreestats.savetsv();
    }

    exit(0);
  }

  // ODT heuristic using HC, BB and DP
  if (flag_hcalgorithm || flag_bestneworks) 
  {

      if (flag_bestneworks) 
      {
        opt_hcstoptime = -1; // no climb employed
      }

      DagSet visiteddags;

      NetworkHCStatsGlobal *globalstats =
        new NetworkHCStatsGlobal(networkclass, timeconsistency, visiteddags, randseed);

      EditOp *editop;
      if (flag_hcedit_nni)
        editop = new NNI();
      else
        editop = new TailMove(networkclass, guideclusters, guidetree);

      if (opt_outfiles.length()) 
      {
        globalstats->setoutfiles(opt_outdirectory, opt_outfiles, odtlabelled);
      }

      int printstats = 0;
      if (flag_hcrunstats) printstats = 1;
      else if (flag_hcrunstatsext) printstats = 2;
      else if (flag_hcrunstatsalways) printstats = 3;

      if (verbosealg >= 4) 
      {
        cout << "HC"
           << " start:"
           << " hcusenaive=" << flag_hcusenaive
           << " runnaiveleqrt=" << hcrunnaiveleqrt_t
           << " tailmove=" << !flag_hcedit_nni << endl;
       }

      NetIterator netiterator(netvec, randomnetworkscnt, 
                      quasiconsensuscnt, genetreeclusters,
                      preserverootst, reticulationcnt_R, 
                      networkclass,
                      timeconsistency, randnetuniform, 
                      guideclusters,
                      guidetree);


     vector<NetworkHCStatsGlobal*> globalstatsarr;

     DagSet visiteddags1;
     DagSet visiteddags2;
     
     if (displaytreesampling.length())
     {
        vector<float> v;
        istringstream iss(displaytreesampling);   
        copy(std::istream_iterator<float>(iss), istream_iterator<float>(), back_inserter(v));
        
        for (auto samplingvalue: v)
        {
          globalstatsarr.push_back(new NetworkHCStatsGlobalSampler(networkclass, timeconsistency, visiteddags2, randseed, samplingvalue));    
        }        
      }

     globalstatsarr.push_back(globalstats);    

     // Test effectiveness of display tree sampler
     if (testdisplaytreesampling)
        {
          for (auto nhc: globalstatsarr) 
          {
             cout << " Sampler " << nhc->getsampling() << " " << nhc->testsampling(testdisplaytreesampling) << endl;
          }
          exit(0);
        }
    
     supnetheuristic(   
        gtvec,       
        &netiterator,
        editop,
        costfun,
        printstats,       
        hcstopinit,
        hcstopclimb,
        flag_hcusenaive,
        hcrunnaiveleqrt_t,    
        hcmaximprovements,
        globalstatsarr,
        flag_cutwhendtimproved
    );

    delete editop;

    if (flag_hcdetailedsummary)
    {
       // print summary stats     
       for (auto nhc: globalstatsarr) 
       {        
          cout << "HC stats: "; 
          nhc->info(cout);
          cout << " ";        
          nhc->print();        
          cout << endl;
       }
    }

    globalstatsarr.pop_back(); // remove the last global

    // merge all data and save odt/dat file(s); optional     
    globalstats->savedatmerged(verbosealg >= 4, globalstatsarr, true);     

    // best dags to file
    globalstats->savebestdags(verbosealg >= 4, true);

    // print summary
    globalstats->print();
    globalstats->printnetworkinfo();     
    cout << endl;


  } 

  if (flag_print_display_trees || flag_print_display_trees_with_ids) 
  {
    for (auto &ntpos: netvec) 
    {
      DISPLAYTREEID tid = 0;      
      RootedTree *t = NULL;
      // cout << *n << endl;
      while ((t = ntpos->gendisplaytree(tid, t)) != NULL) 
      {
        if (flag_print_display_trees_with_ids)
          cout << tid << " ";
        t->printrepr() << endl;
        tid++;
      }
    }
  }

  if (flag_testtreerepr) 
  {
    for (auto &ntpos: netvec) 
    {
      DISPLAYTREEID tid = 0;
      SNode *t = NULL;
      cout << *ntpos << endl;
      while ((t = ntpos->gendisplaytree2(tid, t, globaltreespace)) != NULL) {
        ppSNode(cout, t) << endl;
        tid++;
      }
    }
  }

  if (flag_comparedags) 
  {
    int cnt = 0, cntall = 0;
    for (int i = 0; i < netvec.size(); i++) 
    {      
      Network *n1 = netvec[i];
      cout << endl;
      for (int j = i + 1; j < netvec.size(); j++) 
      {
        Network *n2 = netvec[j];
        bool e1 = n1->eqdags(n2, !flag_dagshapes);
        cout << *(n1) << "\t" << *(n2) << "\t" << e1 << endl;
      }
    }
  }

  if (flag_detectclass) 
  {    
    for (auto &ntpos: netvec) 
    {      
      cout << ntpos->istimeconsistent() <<  " " << ntpos->istreechild() << " " << ntpos->isrelaxed()  << " " << *ntpos << endl;
    }
  }

  // ----------------- DEBUGs -----------------------

  if (flag_testeditnni==1) 
  {
    for (auto &ntpos: netvec) 
    {
      std::ofstream dotf;
      dotf.open("edit.dot", std::ofstream::out);

      int cnt = 0;
      dotf << "digraph SN {" << endl;

      NNI nni;
      nni.init(ntpos);

      cout << *ntpos << " verify==" << ntpos->verifychildparent() << endl;
      ntpos->printdeb(cout, 2);
      ntpos->printdot(dotf, cnt++);

      while (nni.next()) 
      {
        cout << "========" << endl << endl;
        int err = ntpos->verifychildparent();

        cout << *ntpos << " verify==" << err << endl;
        ntpos->printdeb(cout, 2);
        ntpos->printdot(dotf, cnt++);

        if (err)
        {
          break;
        }
      }

      dotf << "}" << endl;
      dotf.close();
    }
  }

  if (flag_testeditnni==2) 
  {
    for (auto &ntpos: netvec) 
    {
      std::ofstream dotf;

      dotf.open("tm.dot", std::ofstream::out);

      int cnt = 0;
      dotf << "digraph SN {" << endl;

      TailMove tailmove(0);
      tailmove.init(ntpos);

      cout << *ntpos << " verify==" << ntpos->verifychildparent() << endl;
      // (*ntpos)->printdeb(cout,2);
      ntpos->printdot(dotf, cnt++);

      while (tailmove.next()) 
      {
        cout << "========" << endl << endl;
        int err = ntpos->verifychildparent();

        cout << ntpos << " verify==" << err << endl;
        // (*ntpos)->printdeb(cout,2);
        ntpos->printdot(dotf, cnt++);

        if (err)
        {
          break;
        }
      }

      dotf << "}" << endl;
      dotf.close();
    }
  }

  // Debug
  if (flag_testcomparedags) 
  {
    int cnt = 0;
    while (1) 
    {
      cnt++;
      string r1 = randspeciestreestr();
      string r2 = randspeciestreestr();
      int EE = 0;
      int e1, e2, e3;
      if (!r1.length() || !r2.length()) {
        cerr << "Cannot create initial random species tree" << endl;
        exit(-1);
      }
      Network *n1 =
          addrandreticulations(reticulationcnt_R, new Network(r1), networkclass,
                               timeconsistency, randnetuniform, NULL, NULL);
      Network *n2 =
          addrandreticulations(reticulationcnt_R, new Network(r2), networkclass,
                               timeconsistency, randnetuniform, NULL, NULL);

      e1 = n1->eqdags(n2);
      e2 = n1->eqdagsbypermutations(n2);

      cout << (*n1) << "\t" << (*n2) << "\tE1=" << e1 << "\tE2=" << e2
           << "\tError=" << (e1 != e2) << endl;

      cerr << (*n1) << "\t" << (*n2) << "\tE1=" << e1 << "\tE2=" << e2
           << "\tError=" << (e1 != e2) << endl;

      delete n1;
      delete n2;

      if (e1 != e2)
        break; // stop
    }
    exit(0);
  }

  if (flag_testcontract) 
  {
    // Last two args: leftretusage rightretusage
    // ./supnet -g "(a,(b,c))" -n "((b)#1,(((#2)#3,((#1)#2,(c,#3))),a))" -eXng 2
    // 1
    // ./dot -Tpdf contr.dot -o c.pdf && evince c.pdf

    RETUSAGE retusage;
    emptyretusage(retusage);

    // ugly
    if (optind + 1 <= argc) {
      long us = atol(argv[optind]);

      int rid = 0;
      while (us) {
        cout << "A" << us << " " << rid << endl;
        if (us & 1)
          addleftretusage(retusage, rid);
        us = us >> 1;
        rid++;
      }
    }

    if (optind + 2 <= argc) {
      long us = atol(argv[optind + 1]);
      int rid = 0;
      while (us) {
        cout << "B" << us << " " << rid << endl;
        if (us & 1)
          addrightretusage(retusage, rid);
        us = us >> 1;
        rid++;
      }
    }

    cout << "R" << retusage << endl;
    cout << "Conflicted " << conflicted(retusage) << endl;
    
    for (int i = 0; i < netvec.size(); i++) 
    {
      Network *n1 = netvec[i];
      std::stringstream ss;
      n1->print(ss);
      cout << ss.str() << endl;

      ContractedNetwork *c = new ContractedNetwork(ss.str());

      c->contract(retusage);

      cout << "Rtcount " << c->rtcount() << " "
           << "rt=" << n1->rtcount() << endl;
      ;

      for (auto &gtpos: gtvec)
      {
        cout << "retmindc:" << (c->approxmindce(*gtpos, *costfun)) << endl;
      }

      ofstream s("contr.dot");
      std::ofstream sf;
      sf.open("contr.dot", std::ofstream::out);
      sf << "digraph SN {" << endl;
      sf << " inp [label=\"InRT=" << retusage << "\"]" << endl;
      c->gendot(sf);
      c->gendotcontracted(sf);
      sf << "}" << endl;
      sf.close();
      cout << c->newickrepr() << endl;

      RootedTree *t = NULL;
      DISPLAYTREEID tid = 0; // id of display tree

      while ((t = c->gendisplaytree(tid, t)) != NULL) {
        cout << tid << " " << (*t) << endl;
        tid++;
      }
    }
  }

  delete globaltreespace;

  for (auto &stpos: sgtvec)  
    free(stpos);

  for (auto &sst: sstvec)  
    free(sst);  

  for (auto &stpos: stvec)  
    free(stpos);

  delete genetreeclusters; // Clean Clusters

  for (auto &gtpos: gtvec)
    delete gtpos;

  for (auto &n: netvec)
    delete n;

  cleanspecies();
}
