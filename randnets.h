#ifndef _RANDNETS_
#define _RANDNETS_

#include "network.h"

Network *addrandreticulations(int reticulationcnt, Network *n, int networktype, int timeconsistency, bool uniform, Clusters *guideclusters, Clusters *guidetree);

class NetIterator
{
  long i;
  int randomnetworkscnt; 
  int quasiconsensuscnt;
  VecNetwork &netvec;
  int reticulationcnt;  
  int timeconsistency;
  int randnetuniform;
  RootedTree *preserverootst;
  Clusters *gtc;
  Clusters *guideclusters;
  Clusters *guidetree;
  int networkclass;

public:

  NetIterator(
    VecNetwork &_netvec, 
    int _randomnetworkscnt, 
    int _quasiconsensuscnt, 
    Clusters *_gtc,
    RootedTree *_preserverootst,
    int _reticulationcnt, 
    int _networkclass, 
    int _timeconsistency, 
    bool _randnetuniform,
    Clusters *_guideclusters,
    Clusters *_guidetree)
  : netvec(_netvec),
    randomnetworkscnt(_randomnetworkscnt),
    quasiconsensuscnt(_quasiconsensuscnt),
    reticulationcnt(_reticulationcnt),    
    timeconsistency(_timeconsistency),
    randnetuniform(_randnetuniform),
    preserverootst(_preserverootst),
    guideclusters(_guideclusters),
    guidetree(_guidetree),
    networkclass(_networkclass),
    gtc(_gtc) 

  {
    i = -1;    
  }

  Network *next();

  long pos() { return i; }

};

Network *randquasiconsnetwork(int reticulationcnt, int networktype, int timeconsistency, Clusters *genetreeclusters, RootedTree *preserverootst, 
  Clusters *guideclusters,
  Clusters *guidetree);


Network *randnetwork(int reticulationcnt, int networktype, int timeconsistency, bool uniform);




#endif