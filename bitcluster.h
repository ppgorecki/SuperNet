
#ifndef _BITCLUSTER_
#define _BITCLUSTER_

#include "tools.h"
#include "clusters.h"


#if MAXSPECIES==64
  typedef uint64_t bitcluster;

#elif MAXSPECIES==128
  typedef unsigned __int128 bitcluster; 
  
#else
   // todo  
#endif

extern bitcluster bcsingleton[MAXSPECIES];

#if MAXSPECIES<=128

	#define INTERSECTION(a,b) ((a)&(b))
	
	// a is a subset of b
	#define ISSUBSET(a,b) (INTERSECTION((a),(b))==(a))

	//
	#define EQSETS(a,b) ((a)==(b))
	
	// set difference
	#define DIFFERENCE(a,b) ((a) & (~(b)))

	// set union
	#define UNION(a,b) ((b)|(a))

	// empty set
	#define EMPTYSET() (bitcluster)0

	// add to set
	#define SETADD(a,n)  ((a)|bcsingleton[n])

	// add to set
	#define SETINSERT(a,n)  a|=bcsingleton[n]
	
	// has an elt
	#define SETHAS(a,n)  ((a)& bcsingleton[n])

	// remove from set  
	#define SETDELETE(a,n) ((a)&~bcsingleton[n])	
	
#else
	// todo
#endif 

// print cluster in 0-1 repr
std::ostream &ppbitcluster(std::ostream &os, bitcluster c);

// print elements in cluster
std::ostream &ppbitclusterspecies(std::ostream &os, bitcluster c);

#endif

void clustergraphtester();	
