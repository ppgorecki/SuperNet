#ifndef _TESTERS_H_
#define _TESTERS_H_

#include "network.h"


void testcomparedags(int reticulationcnt_R, int timeconsistency, int randnetuniform, int networkclass);
void testcontract(int argc, char **argv, VecRootedTree &gtvec,  VecNetwork &netvec, CostFun *costfun);
void testeditnni(VecNetwork &netvec);
void testedittailmove(VecNetwork &netvec);
void testtreerepr(VecNetwork &netvec);
void testnetworkretiterator(VecNetwork &netvec, int networkclass, int timeconsistency, Clusters *guideclusters, Clusters *guidetree);


#endif