#ifndef _RANDNETS_
#define _RANDNETS_

#include "network.h"

Network *addrandreticulations(int reticulationcnt, Network *n, int networktype, int timeconsistency, bool uniform);


Network* netiterator(long int &i, VecNetwork &netvec, int &randomnetworkscnt, int &quasiconsensuscnt, 
  TreeClusters *gtc,
  RootedTree *preserverootst,
  int reticulationcnt, int networktype, int timeconsistency, bool randnetuniform);

Network *randquasiconsnetwork(int reticulationcnt, int networktype, int timeconsistency, TreeClusters *gtc, RootedTree *preserverootst);


Network *randnetwork(int reticulationcnt, int networktype, int timeconsistency, bool uniform);

#endif