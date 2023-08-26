#ifndef _RANDNETS_
#define _RANDNETS_

#include "network.h"

Network *addrandreticulations(int reticulationcnt, Network *n, int networktype, int timeconsistency, bool uniform, Clusters *guideclusters, Clusters *guidetree);


Network* netiterator(long int &i, VecNetwork &netvec, int &randomnetworkscnt, int &quasiconsensuscnt, 
  Clusters *gtc,
  RootedTree *preserverootst,
  int reticulationcnt, int networktype, int timeconsistency, bool randnetuniform,
  Clusters *guideclusters,
  Clusters *guidetree);

Network *randquasiconsnetwork(int reticulationcnt, int networktype, int timeconsistency, Clusters *genetreeclusters, RootedTree *preserverootst, 
  Clusters *guideclusters,
  Clusters *guidetree);


Network *randnetwork(int reticulationcnt, int networktype, int timeconsistency, bool uniform);

#endif