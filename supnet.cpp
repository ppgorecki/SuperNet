
/************************************************************************
 SuperNetwork Inference - based partially on FastUrec v2.05
(c) Copyright 2005-2023 by Pawel Gorecki
 Written by P.Gorecki.

 Permission is granted to copy and use this program provided no fee is
 charged for it and provided that this copyright notice is not removed.
 *************************************************************************/

const char *SUPNET = "0.13";

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
#include "netgen.h"
#include "rtree.h"
#include "tools.h"
#include "topsort.h"
#include "treespace.h"
#include "testers.h"
#include "algotok.h"
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
double opt_genetreessimilarity = 0.7;

int flag_globaldagcache = 0; 
int flag_saveextnewick = 0; 
int flag_autooutfiles = 0; 
int flag_savewhenimproved = 0;
int flag_hcsamplerstats = 0;
int flag_hcsamplingmaxnetstonextlevel = 3; 
int flag_hcdetailedsummarydat = 0;
float opt_hcstoptime = 0; 
bool flag_hcignorecostgeq = 0; // if set 
float opt_hcignorecostgeq = 0;
int flag_hcsavefinalodt = 0; // ignore hcignorecostgeq for final files
int maxdisplaytreecachesize = 1000000;

int randomnetworkscnt = 0;
int timeconsistency = NET_ANY;

int print_repr_inodtnaive = 0;

string flag_retidprefix = "";

TreeSpace *globaltreespace;

#define BUFSIZE 100000

// Read trees from a given file
void readtrees(char *fn, vector<char *> &speciestreesv) 
{ 
  FILE *f;
  if (!strcmp(fn, "-"))
    f = stdin;
  else
    f = fopen(fn, "r");
  if (!f) 
  {
    cerr << "Cannot open file " << fn << endl;
    exit(-1);
  }
  char buf[BUFSIZE];
  while (1) {
    if (!fgets(buf, BUFSIZE, f))
      break;
    speciestreesv.push_back(strdup(buf));
  }
  fclose(f);
}

string readscript(string fn) 
{    
  std::ifstream ifs(fn);
  return string((std::istreambuf_iterator<char>(ifs)),
                  (std::istreambuf_iterator<char>()));
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



int main(int argc, char **argv) {
  int opt;

  if (argc < 2)
    usage(argc, argv);

  VecRootedTree speciestreesv, genetreesv;
  VecNetwork networksv;

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
  int flag_comparedags = 0;
  int flag_dagshapes = 0;
  int flag_uniquedags = 0;
  int flag_printstsubtrees = 0; // of s
  int flag_uniquedagscnts = 0;

  int flag_preserveroot = 0;
  

  int flag_extra_visible_node_stats = 0;
  int flag_extra_visible_leaf_stats = 0;
  int flag_dpalgorithm = 0;  
  int flag_bbalgorithm = 0;
  int flag_bbalgorithmSTATS = 0;

  int flag_odt_naive_gtvsnet = 0;

  int quasiconsensuscnt = 0;

  int reticulationcnt_R = 0;
  int networkclass = NET_TREECHILD; // default network type

  int runnaiveleqrt_t = 13; // default for DC, based on experiments

  int odtlabelled = 0;

  int randnetuniform = 0;

  
  
  vector<char *> sgenetreesv, sspeciestreesv, snetworksv;

  COSTT bbstartscore = 0;
  int bbstartscoredefined = 0;
  string displaytreesampling = "";

  

  int flag_hcusenaive = 0;
  

  int net_treechild = 1;
  int net_relaxed = 0;
  int net_general = 0;

  int flag_timeconsistent = 0;
  int flag_notimeconsistent = 0;
  timeconsistency = NET_ANY;  // NET_ANY no restrictions, 1-timeconsistent, 2-notimeconsistent
  
  int flag_ptreecostext = 0;
  int flag_noodtfiles = 0;
  int flag_bbtreesearch = 0;
  int flag_bbtsvstats = 0;
  int flag_bbtimestats = 0;
  int flag_hcalgorithm = 0;
  int flag_runnaiveleqrt = 0;
  int flag_hcedit_nni = 0;
  int flag_detectclass = 0;
  int flag_hcrunstats = 0;
  int flag_hcrunstatsext = 0;
  int flag_hcrunstatsalways = 0;
  int flag_hcdetailedsummary = 0;
  int flag_match_after = 0;
  int flag_match_before = 0;
  int flag_match_atpos = 0;
  int flag_print_network_clusters = 0;
  int flag_cutwhendtimproved = 0;
  int flag_netretiterator = 0;
  int flag_iterativeretinsertion = 0;
  int flag_pnetworkretnodescnt = 0;
  int testdisplaytreesampling = 0;

  int flag_icalgorithm = 0;

  int hcmaximprovements = 0;
  int hcstopinit = 0;
  int hcstopclimb = 0;
  int flag_odtcost = 0;
  int flag_bestneworks = 0;
  int pranddisplaytrees = 0;
  char *matchparam = NULL;

  char *opt_guideclusters = NULL;
  char *opt_guidetree = NULL;

  string opt_tester = "";
  string opt_script = "";
  
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
      {"pranddisplaytrees", required_argument, NULL, 'J'},
      
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
      {"pnetworkretnodescnt", no_argument, &flag_pnetworkretnodescnt, 1},
      {"pdetailed", no_argument, &flag_extra_print_detailed, 1},
      {"netretiterator", no_argument, &flag_netretiterator, 1},

      {"DP", no_argument, &flag_dpalgorithm, 1},

      {"odtnaivecost", no_argument, &flag_odt_naive_gtvsnet, 1},
      {"cost", required_argument, NULL, 'C'},
      {"odtcost", no_argument, &flag_odtcost, 1},

      {"BB", no_argument, &flag_bbalgorithm, 1},
      {"bbtreesearch", no_argument, &flag_bbtreesearch, 1},
      {"bbtsvstats", no_argument, &flag_bbtsvstats, 1},
      {"bbtimestats", no_argument, &flag_bbtimestats, 1},
      {"bbstartscore", required_argument, NULL, '1'},

      {"iterativeretinsertion", no_argument, &flag_iterativeretinsertion, 1},

      {"HC", optional_argument, &flag_hcalgorithm, 'o'},
      {"runnaiveleqrt", required_argument, &flag_runnaiveleqrt, 't'},
      {"hcusenaive", no_argument, &flag_hcusenaive, 1},
      {"hcnnimove", no_argument, &flag_hcedit_nni, 1},
      {"hcrunstats", no_argument, &flag_hcrunstats, 1},
      {"hcsamplerstats", no_argument, &flag_hcsamplerstats, 1},
      {"hcrunstatsext", no_argument, &flag_hcrunstatsext, 1},
      {"hcrunstatsalways", no_argument, &flag_hcrunstatsalways, 1},
      {"savewhenimproved", no_argument, &flag_savewhenimproved, 1},
      {"hcsamplingmaxnetstonextlevel", required_argument, NULL, '3'},
      {"hcdetailedsummary", no_argument, &flag_hcdetailedsummary, 1},
      {"hcdetailedsummarydat", no_argument, &flag_hcdetailedsummarydat, 1},
      {"autooutfiles", no_argument, &flag_autooutfiles, 1},
      {"cutwhendtimproved", no_argument, &flag_cutwhendtimproved, 1},

      {"IC", optional_argument, &flag_icalgorithm, 'i'},

      {"bestnetworks", no_argument, &flag_bestneworks, 1},
      
      {"hcmaximprovements", required_argument, NULL, 'F'},
      {"hcignorecostgeq", required_argument, NULL, '4'},
      {"hcsavefinalodt", no_argument, &flag_hcsavefinalodt, 1},
      {"hcstopinit", required_argument, NULL, 'Y'},
      {"hcstopclimb", required_argument, NULL, 'Z'},
      {"hcstoptime", required_argument, NULL, 'P'},

      {"run",required_argument, NULL, 'E'},

      // TODO: implement rules
      {"matchafter", no_argument, &flag_match_after, 1},
      {"matchbefore", no_argument, &flag_match_before, 1},
      {"matchatpos", no_argument, &flag_match_atpos, 1},
      {"matchparam", required_argument, NULL, 'l'},


      // debug
      {"test", required_argument, NULL, 'T'},
      {"guideclusters", required_argument, NULL, 'u'},
      {"guidetree", required_argument, NULL, 'e'},
      {"globaldagcache", no_argument, &flag_globaldagcache, 1 },
      {"saveextnewick", no_argument, &flag_saveextnewick, 1 },
      {"retidprefix", required_argument, NULL, 'Q'},
      {"genetreessimilarity", required_argument, NULL, 'p'},

      {0, 0, 0, 0}

  };

  const char *optstring =
      "c:z:b:U:g:s:G:S:N:n:R:t:m:V:X:Y:Z:l:q:r:A:L:D:O:v:C:1:Wu:e:E:P:T:Q:u:p:i:;";

  while ((opt = getopt_long(argc, argv, optstring, longopts, NULL)) != -1) 
  {

    
    switch (opt) {
    case 'V':
      maxdisplaytreecachesize = atoi(optarg);
      break;

    case 'z':
      randseed = (unsigned int)atoi(optarg);      
      break;

    case 'J':
      pranddisplaytrees = (unsigned int)atoi(optarg);      
      break;

    case '1':
      bbstartscore = atof(optarg);
      break;

    case 'p':
      opt_genetreessimilarity = atof(optarg);
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
      insertstr(sgenetreesv, optarg);
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
      insertstr(sspeciestreesv, optarg);
      break;

    case 'n':
      insertstr(snetworksv, optarg);
      break;

    case 'G': {
      readtrees(optarg, sgenetreesv);
      break;
    }

    case 'O': 
    {
      opt_outfiles = optarg;      
      break;
    }

    case 'T': 
    {
      opt_tester = optarg;      
      break;
    }

    case 'E': 
    {
      opt_script = readscript(optarg);

      break;
    }

    case 'd': 
    {
      opt_outdirectory = optarg;      
      break;
    }

    case 'S': {
      readtrees(optarg, sspeciestreesv);
      break;
    }

    case 'N': {
      readtrees(optarg, snetworksv);
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
        cerr << "Number expected in -q" << endl;
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
    case 'A': 
      {
        int spcnt = 0;
        if (sscanf(optarg, "%d", &spcnt) != 1) {
          cerr << "Number expected in -A" << endl;
          exit(-1);
        }
        setspecies(spcnt);      
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
    } // switch 

    // cout << ((char)opt) << " " << optarg << ' ' << flag_print_networks << endl;


  } // while

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
    timeconsistency = NET_TIMECONSISTENT;
  }
  else 
    if (flag_notimeconsistent) 
    {
      timeconsistency = NET_NOTIMECONSISTENT;
    }

  // Parse species trees
  for (size_t i = 0; i < sspeciestreesv.size(); i++) 
  {
    RootedTree *s = new RootedTree(sspeciestreesv[i]);
    if (!s->bijectiveleaflabelling()) {
      cerr << "Bijective leaf labelling expected in a species tree: " << *s
           << endl;
      exit(-1);
    }
    speciestreesv.push_back(s);
  }

  // Parse gene trees
  for (size_t i = 0; i < sgenetreesv.size(); i++) 
  {    
    RootedTree *gtree = new RootedTree(sgenetreesv[i]);
    gtree->setid(i);
    genetreesv.push_back(gtree);
  }

  // Print species names
  if (flag_extra_print_species_dictionary)
    for (size_t i = 0; i < specnames2id.size(); i++)
      cout << i << " " << specnames[i] << endl;

  Clusters *genetreeclusters = NULL; 

  // Prepare clusters
  if (opt_quasiconsensus || opt_randnetworks) 
  {
    genetreeclusters = new Clusters(genetreesv);     
  }

  // Initialize clusters
  globaltreespace = new TreeSpace(genetreesv, maxdisplaytreecachesize);

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

  // Initialize networks  
  for (size_t i = 0; i < snetworksv.size(); i++) 
  {
    Network *n = new Network(snetworksv[i]);

    if (!n->bijectiveleaflabelling()) 
    {
      cerr << "Bijective leaf labelling expected in a network: " << *n << endl;
      exit(-1);
    }

    networksv.push_back(addrandreticulations(reticulationcnt_R, n, networkclass,
                                           timeconsistency, randnetuniform, NULL, NULL));
  }


  // Gen quasi consensus trees and insert into networksv as networks
  // Add reticulations if -R is set
  // Ignored if flag_hcalgorithm is set
  RootedTree *preserverootst = NULL;
  if (opt_quasiconsensus) 
  {
    if (flag_preserveroot) 
    {
      if (speciestreesv.size() > 0)
        preserverootst = speciestreesv[0];
      else 
      {
        cout << "Exactly one species tree should be defined for "
                "quasi-consensus trees with preserve-root option"
             << endl;
        exit(2);
      }
    }

    if (!genetreesv.size() && !guideclusters && !guidetree) 
    {
        cout << "Gene trees (-g/-G), a guide tree (--guidetree) or guide clusters (--guideclusters) are required to infer quasi consensus" << endl;
        return -1;
    }

    if (!flag_hcalgorithm && !flag_bestneworks && !opt_script.length()) // quasi cons. generated separately
    {
      
      for (int i = 0; i < quasiconsensuscnt; i++)
        networksv.push_back(randquasiconsnetwork(reticulationcnt_R, networkclass,
                                              timeconsistency, genetreeclusters,
                                              preserverootst, 
                                              guideclusters,
                                              guidetree));
    } 
  }

  // Stats on random networks generated on the fly
  // Ignores other options and exits when completed
  if (flag_uniquedags || flag_uniquedagscnts) 
  {
    DagSet dagset(flag_uniquedagscnts, flag_dagshapes);
    // get next network
    Network *n;
    long int i = -1;

    Dag *src;
    NetGenerator netgenerator(&networksv, randomnetworkscnt, quasiconsensuscnt, genetreeclusters,
                          preserverootst, reticulationcnt_R, networkclass,
                          timeconsistency, randnetuniform, guideclusters, 
                           guidetree);
    while ((n = netgenerator.next()) != NULL)
       dagset.add(n, &src);

    cout << dagset;
    cerr << "unique=" << dagset.size() << " all=" << networksv.size() << endl;
    exit(0); 
    // ignore other options
  }

  // Print species trees
  if (flag_print_species_trees) 
  {
    for (auto & stpos : speciestreesv) 
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

  // Print subtrees of species trees
  if (flag_printstsubtrees) 
  {
    for (auto &stpos : speciestreesv) 
    {
      stpos->printsubtrees(cout);
    }
  }

  // Print gene trees
  if (flag_print_gene_trees) 
  {
    for (auto gtpos : genetreesv) 
    {
      cout << *gtpos << endl;
    }
  }

  // Gene tree vs. species tree cost
  if (flag_ptreecost || flag_ptreecostext) 
  {
    for (auto &stpos : speciestreesv) 
    {
      for (auto &gtpos : genetreesv) 
      {
        if (flag_ptreecostext) 
        {
          cout << *stpos << " " << *gtpos << " ";
        }

        cout << costfun->computegt(*gtpos, *stpos) << endl;
      }
    }
  }

  // Gen random trees and store in networksv
  // Add reticulations if -R is set
  if (opt_randnetworks && !flag_hcalgorithm) 
  {

    for (int i = 0; i < randomnetworkscnt; i++)
    {
      networksv.push_back(randnetwork(reticulationcnt_R, networkclass,
                                   timeconsistency, randnetuniform));
    }
  }

  // Print networks (opt. clusters)
  if (flag_print_networks) 
  {
    for (auto &ntpos: networksv) 
    {
      //(*ntpos)->printdeb(cout,2) << endl;
      cout << *ntpos << endl;
      if (flag_print_network_clusters)
      {
        ntpos->printclusters(cout);
      }
    }
  }

  // Print retnodes counts in networks 
  if (flag_pnetworkretnodescnt) 
  {
    for (auto &ntpos: networksv) 
    {      
      cout << (*ntpos).rtcount() << endl;
      
    }
  }
  
  // More printing
  if (flag_extra_print_detailed) 
  {
    if (speciestreesv.size())
    {
      cout << "Species trees:" << endl;
    }

    for (auto &stpos: speciestreesv)
    {
      stpos->printdeb(cout, 2) << endl;
    }

    if (genetreesv.size())
    {
      cout << "Gene trees:" << endl;
    }

    for (auto &gtpos: genetreesv)
    {
      gtpos->printdeb(cout, 2) << endl;
    }

    if (networksv.size())
    {
      cout << "Networks:" << endl;
    }

    for (auto &ntpos: networksv)
    {
      ntpos->printdeb(cout, 2) << endl;
    }
  }

  // Print stats on visible nodes in networks
  if (flag_extra_visible_node_stats) 
  {
    for (auto &ntpos: networksv)
    {
      ntpos->visibilenodestats(1, cout) << endl;
    }
  }

  // Print stats on visible leaves in networks
  if (flag_extra_visible_leaf_stats) 
  {
    for (auto &ntpos: networksv)
      ntpos->visibilenodestats(2, cout) << endl;
  }

  // For every network print all network after inserting one reticulation
  if (flag_netretiterator) 
  {
      for (auto &ntpos: networksv) 
      {        
        NetworkRetIterator netretit((*ntpos), networkclass, timeconsistency, guideclusters, guidetree, "");

        Network *n;
        while ((n = netretit.next())!=NULL)
        {
          cout << *n << endl;
        }
      }
  }

  // Generate dot file with networks and trees
  if (flag_dot) 
  {    
    ostream &s = cout; // todo files
    int dagnum = 0;
    s << "digraph SN {" << endl;

    for (auto &stpos: speciestreesv)
    {
      stpos->printdot(s, dagnum++) << endl;
    }

    for (auto &gtpos: genetreesv)
    {
      gtpos->printdot(s, dagnum++) << endl;
    }

    for (auto &ntpos: networksv)
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
    for (auto &ntpos: networksv) 
    {
      if (flag_odt_naive_gtvsnet == 2) 
      {
        cout << *ntpos << " ";
      }
      cout << ntpos->odtcostnaive(genetreesv, *costfun, odtstats, 0)
           << endl;
    }
  }

  // Compute ODT cost 
  if (flag_odtcost) 
  {
    DagSet *visiteddags = new DagSet();
    ClimbStatsBase *stats =
        new ClimbStatsBase(networkclass, timeconsistency, visiteddags, randseed);
    for (auto &net: networksv) 
    {
      double cost = net->odtcost(genetreesv, *costfun, flag_hcusenaive, runnaiveleqrt_t, stats->getodtstats());
      cout << cost << " " << *net << endl;
    }
    exit(0);
  }

  // Run DP algorithm to compute approx DCE
  if (flag_dpalgorithm) 
  {
    for (auto &ntpos: networksv) 
    {
      for (auto &gtpos: genetreesv) 
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

    for (auto &ntpos: networksv)
    {
      for (auto &gtpos: genetreesv) 
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

  int printstats = 0;
  if (flag_hcrunstats) printstats = 1;
  else if (flag_hcrunstatsext) printstats = 2;
  else if (flag_hcrunstatsalways) printstats = 3;

   

  // ODT heuristic using HC, BB and DP
  if (flag_hcalgorithm || flag_bestneworks) 
  {

      if (flag_bestneworks) 
      {
        opt_hcstoptime = -1; // no climb employed
      }

      ClimbStatsGlobal *globalstats =
        new ClimbStatsGlobal(networkclass, timeconsistency, new DagSet, randseed, flag_savewhenimproved);

      EditOp *editop;
      if (flag_hcedit_nni)
        editop = new NNI();
      else
        editop = new TailMove(networkclass, guideclusters, guidetree);

      if (opt_outfiles.length()) 
      {                
        globalstats->setoutfiles(opt_outdirectory, opt_outfiles, odtlabelled);
      }



      if (verbosealg >= 4) 
      {
        cout << "HC"
           << " start:"
           << " hcusenaive=" << flag_hcusenaive
           << " runnaiveleqrt=" << runnaiveleqrt_t
           << " tailmove=" << !flag_hcedit_nni << endl;
       }

      NetGenerator netgenerator(
                      &networksv, 
                      randomnetworkscnt, 
                      quasiconsensuscnt, 
                      genetreeclusters,
                      preserverootst, 
                      reticulationcnt_R, 
                      networkclass,
                      timeconsistency, 
                      randnetuniform, 
                      guideclusters,
                      guidetree);


     vector<ClimbStatsGlobal*> globalstatsarr;     
     
     if (displaytreesampling.length())
     {
        vector<float> v;
        istringstream iss(displaytreesampling);   
        copy(std::istream_iterator<float>(iss), istream_iterator<float>(), back_inserter(v));
        
        for (auto samplingvalue: v)
        {
          globalstatsarr.push_back(new ClimbStatsGlobalSampler(networkclass, 
            timeconsistency, 
            new DagSet(), randseed, samplingvalue));    
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
        genetreesv,       
        &netgenerator,
        editop,
        costfun,
        printstats,       
        hcstopinit,
        hcstopclimb,
        flag_hcusenaive,
        runnaiveleqrt_t,    
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

  if (flag_print_display_trees || (flag_print_display_trees_with_ids && !pranddisplaytrees)) 
  {
    for (auto &ntpos: networksv) 
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

  if (pranddisplaytrees>0) 
  {
    for (auto &ntpos: networksv) 
    {
      DISPLAYTREEID mxid = ntpos->displaytreemaxid();
      RootedTree *t = NULL;
      for (int i=0; i<pranddisplaytrees; i++)
      {      

        DISPLAYTREEID tid = rand() % mxid;      

        t = ntpos->gendisplaytree(tid, t);

        if (flag_print_display_trees_with_ids)
            cout << tid << " ";
        t->printrepr() << endl;
      }
    }
  }

  if (flag_comparedags) 
  {
    int cnt = 0, cntall = 0;
    for (int i = 0; i < networksv.size(); i++) 
    {      
      Network *n1 = networksv[i];
      cout << endl;
      for (int j = i + 1; j < networksv.size(); j++) 
      {
        Network *n2 = networksv[j];
        bool e1 = n1->eqdags(n2, !flag_dagshapes);
        cout << *(n1) << "\t" << *(n2) << "\t" << e1 << endl;
      }
    }
  }

  if (flag_detectclass) 
  {    
    for (auto &ntpos: networksv) 
    {      
      cout << ntpos->istimeconsistent() <<  " " << ntpos->istreechild() << " " << ntpos->isrelaxed()  << " " << *ntpos << endl;
    }
  }


  if (opt_script.length())
  {
       
      // AlgTokenizer s(opt_script);
      // int tokentype;
      

      // while ((tokentype = s.next())!=0)
      // {   
      //     string token = s.get();
      //     cout << tokentype << ": " << token << endl;

      //     if (tokentype==TLABEL && token=="randnetwork")
      //     {
      //         command = token;
      //     }
      // }

    ExprList *e = parsescript(opt_script);

    Env *env = new Env(genetreesv, 
        networksv,
        costfun,                                    
        guideclusters,
        guidetree,
        preserverootst
        );

    env->set("reticulations",reticulationcnt_R);
    env->set("networkclass",networkclass);
    env->set("timeconsistency", timeconsistency);
    env->set("randnetuniform", randnetuniform);
    env->set("randomnetworks", randomnetworkscnt);
    env->set("quasiconsensusnetworks", quasiconsensuscnt);
    env->set("cutwhendtimproved", flag_cutwhendtimproved);
    env->set("runnaiveleqrt", runnaiveleqrt_t);
    env->set("usenaive", flag_hcusenaive);
    env->set("printstats", printstats);
    env->set("randseed", randseed);
    env->set("outdirectory", opt_outdirectory);
    env->set("outfiles",         opt_outfiles);
    env->set("odtlabelled", odtlabelled);
    env->set("verbosealg",verbosealg);
    env->set("verbosehccost",verbosehccost);
    env->set("genetreessimilarity",opt_genetreessimilarity);

    env->set("stopclimb",hcstopclimb);
    env->set("stopinit",hcstopinit);
    env->set("hceditnni",flag_hcedit_nni);
    env->set("maximprovements",hcmaximprovements);
    env->set("hcdetailedsummary",flag_hcdetailedsummary);    
    env->set("verbosealg", verbosealg);
    env->set("savewhenimproved", flag_savewhenimproved);
        
    e->eval(*env);

    env->finalizestats();
        

  }

  // ----------------- DEBUGs -----------------------

  if (opt_tester=="editnni") 
  {
    testeditnni(networksv);
  }

  if (opt_tester=="edittailmove") 
  {
    testedittailmove(networksv);
  }

  if (opt_tester=="comparedags") 
  {
    testcomparedags(reticulationcnt_R, timeconsistency, randnetuniform, networkclass);

  }

  if (opt_tester=="contract") 
  {
    testcontract(argc, argv, genetreesv, networksv, costfun);
  }

  if (opt_tester=="treerepr") 
  {
    testtreerepr(networksv);    
  }

  if (opt_tester=="iiopt")
  {
    
    DagSet *visiteddags = new DagSet();
    ClimbStatsGlobal *globalstats =
       new ClimbStatsGlobal(networkclass, timeconsistency, visiteddags, randseed, flag_savewhenimproved);
        
    if (opt_outfiles.length()) 
    {
        globalstats->setoutfiles(opt_outdirectory, opt_outfiles, odtlabelled);
    }

    for (int i = 0; i < networksv.size(); i++) 
    {
        // cout << *networksv[i] << endl;
        iterativeretinsertionoptimizer(    
        genetreesv,
        networksv[i],
        costfun,   
        printstats,           
        flag_hcusenaive,
        runnaiveleqrt_t,        
        globalstats, // could be sampler
        flag_cutwhendtimproved,
        networkclass, 
        timeconsistency, 
        guideclusters,
        guidetree          
      );
    }

    vector<ClimbStatsGlobal*> globalstatsarr;

    // merge all data and save odt/dat file(s); optional     
    globalstats->savedatmerged(verbosealg >= 4, globalstatsarr, true);     

    // best dags to file
    globalstats->savebestdags(verbosealg >= 4, true);

    // print summary
    globalstats->print();
    globalstats->printnetworkinfo();     
    cout << endl;
  }

  if (opt_tester=="rnets")
  {
    Network *n; 

    // QuasiConsTreeGenerator qct(10, preserverootst, guideclusters, guidetree, new Clusters(genetreesv));
    // while ((n=qct.next())!=0)
    //   cout << *n << endl;

    QuasiConsTreeGenerator qct2(10, preserverootst, guideclusters, guidetree, new Clusters(genetreesv), opt_genetreessimilarity);

    NetRetGenerator rn(&qct2, reticulationcnt_R, NET_TREECHILD, timeconsistency, randnetuniform, -1, 1, guideclusters, guidetree);
    
    while ((n=rn.next())!=0)
      cout << *n << endl;


  }

    //  supnetheuristic(   
    //     genetreesv,       
    //     &netgenerator,
    //     editop,
    //     costfun,
    //     printstats,       
    //     hcstopinit,
    //     hcstopclimb,
    //     flag_hcusenaive,
    //     runnaiveleqrt_t,    
    //     hcmaximprovements,
    //     globalstatsarr,
    //     flag_cutwhendtimproved
    // );

  


  delete globaltreespace;

  for (auto &stpos: sgenetreesv)  
    free(stpos);

  for (auto &sst: sspeciestreesv)  
    free(sst);  

  for (auto &stpos: speciestreesv)  
    free(stpos);

  delete genetreeclusters; // Clean Clusters

  for (auto &gtpos: genetreesv)
    delete gtpos;

  for (auto &n: networksv)
    delete n;

  cleanspecies();

  return 0;
}
