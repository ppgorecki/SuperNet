
#ifndef _BITCLUSTER_
#define _BITCLUSTER_

#include "tools.h"
#include "clusters.h"

void initbitclusters();

#if MAXSPECIES==64
  typedef uint64_t bitcluster;

#elif MAXSPECIES==128
  typedef unsigned __int128 bitcluster; 

#elif MAXSPECIES==256
  typedef struct bitcluster { unsigned __int128 x, y; } bitcluster;
  
#else
   #error Macro MAXSPECIES must be 64, 128 or 256. Other values are not implemented yet.
#endif

extern bitcluster bcsingleton[MAXSPECIES];

#if MAXSPECIES<=128

	#define NONEMPTYINTERSECTION(a,b) ((a)&(b))	
	
	// a is a subset of b
	#define ISSUBSET(a,b) (((a)&(b))==(a))

	#define EQSETS(a,b) ((a)==(b))
	
	// set difference
	#define DIFFERENCE(a,b) ((a) & (~(b)))

	// set union
	#define UNION(a,b) ((a)|(b))

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
	
#elif MAXSPECIES==256

  #define andx(a,b) (((a).x)&((b).x))
  #define andy(a,b) (((a).y)&((b).y))
	
	#define NONEMPTYINTERSECTION(a,b) (andx(a,b) || andy(a,b))
	
	// a is a subset of b
	#define ISSUBSET(a,b) ( (andx(a,b)==((a).x)) && (andy(a,b)==((a).y)) )

	// set union
	#define UNION(a,b) {((a).x)|((b).x),((a).y)|((b).y)}

	// empty set
	#define EMPTYSET() (bitcluster){0,0}

	// add to set
	#define SETINSERT(a,n)  { a.x|=bcsingleton[n].x; a.y|=bcsingleton[n].y; } 
	
	// has an elt
	#define SETHAS(a,n)  (((a).x & bcsingleton[n].x) | ((a).y & bcsingleton[n].y ) )

	// remove from set  
	// #define SETDELETE(a,n) ((a)&~bcsingleton[n])	
	// TODO	
#endif 

#define EMPTYINTERSECTION(a,b) (!NONEMPTYINTERSECTION(a,b))
// print cluster in 0-1 repr
std::ostream &ppbitcluster(std::ostream &os, bitcluster c);

// print elements in cluster
std::ostream &ppbitclusterspecies(std::ostream &os, bitcluster c);

#endif

void clustergraphtester();	
