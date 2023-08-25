#ifndef _RANDNETS_
#define _RANDNETS_

#include "network.h"

Network *addrandreticulations(int reticulationcnt, Network *n, int networktype, int timeconsistency, bool uniform, Clusters *guidetreeclusters);


Network* netiterator(long int &i, VecNetwork &netvec, int &randomnetworkscnt, int &quasiconsensuscnt, 
  Clusters *gtc,
  RootedTree *preserverootst,
  int reticulationcnt, int networktype, int timeconsistency, bool randnetuniform,
  Clusters *guidetree);

Network *randquasiconsnetwork(int reticulationcnt, int networktype, int timeconsistency, Clusters *genetreeclusters, RootedTree *preserverootst, Clusters *guidetreeclusters);


Network *randnetwork(int reticulationcnt, int networktype, int timeconsistency, bool uniform);

#endif