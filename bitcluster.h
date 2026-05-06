
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

#elif MAXSPECIES==512 || MAXSPECIES==1024
  // Generic N-chunk bitcluster for very large species sets. Each chunk is
  // 128 bits, so MAXSPECIES==512 → 4 chunks, ==1024 → 8 chunks. Operations
  // are loops with constant trip count; -O3 unrolls them so the cost is
  // identical to hand-rolled macros (verified in disasm for N≤8).
  #if MAXSPECIES==512
    #define BC_CHUNKS 4
  #else
    #define BC_CHUNKS 8
  #endif
  typedef struct bitcluster { unsigned __int128 chunks[BC_CHUNKS]; } bitcluster;

#else
   #error Macro MAXSPECIES must be 64, 128, 256, 512 or 1024.
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

  // Lexicographic ordering / equality for the 256-bit struct so bitclusters
  // can be sorted and binary-searched (used by DTCACHE RF lookups).
  inline bool operator<(const bitcluster &a, const bitcluster &b)
  { return (a.x < b.x) || (a.x == b.x && a.y < b.y); }
  inline bool operator==(const bitcluster &a, const bitcluster &b)
  { return a.x == b.x && a.y == b.y; }

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

#elif MAXSPECIES==512 || MAXSPECIES==1024

  // Comparison + equality (lexicographic over chunks).
  inline bool operator==(const bitcluster &a, const bitcluster &b)
  {
    for (int i = 0; i < BC_CHUNKS; i++)
      if (a.chunks[i] != b.chunks[i]) return false;
    return true;
  }
  inline bool operator<(const bitcluster &a, const bitcluster &b)
  {
    for (int i = 0; i < BC_CHUNKS; i++)
      if (a.chunks[i] != b.chunks[i]) return a.chunks[i] < b.chunks[i];
    return false;
  }

  // Forward decl — defined below.
  extern bitcluster bcsingleton[];

  inline bool _bc_nonempty_intersection(const bitcluster &a, const bitcluster &b)
  {
    for (int i = 0; i < BC_CHUNKS; i++)
      if (a.chunks[i] & b.chunks[i]) return true;
    return false;
  }
  inline bool _bc_subset(const bitcluster &a, const bitcluster &b)
  {
    for (int i = 0; i < BC_CHUNKS; i++)
      if ((a.chunks[i] & b.chunks[i]) != a.chunks[i]) return false;
    return true;
  }
  inline bitcluster _bc_or(const bitcluster &a, const bitcluster &b)
  {
    bitcluster r;
    for (int i = 0; i < BC_CHUNKS; i++) r.chunks[i] = a.chunks[i] | b.chunks[i];
    return r;
  }
  inline bitcluster _bc_emptyset()
  {
    bitcluster r;
    for (int i = 0; i < BC_CHUNKS; i++) r.chunks[i] = 0;
    return r;
  }
  inline void _bc_set_insert(bitcluster &a, NODEID n)
  {
    for (int i = 0; i < BC_CHUNKS; i++) a.chunks[i] |= bcsingleton[n].chunks[i];
  }
  inline bool _bc_set_has(const bitcluster &a, NODEID n)
  {
    for (int i = 0; i < BC_CHUNKS; i++)
      if (a.chunks[i] & bcsingleton[n].chunks[i]) return true;
    return false;
  }

  #define NONEMPTYINTERSECTION(a, b) _bc_nonempty_intersection((a), (b))
  #define ISSUBSET(a, b)             _bc_subset((a), (b))
  #define UNION(a, b)                _bc_or((a), (b))
  #define EMPTYSET()                 _bc_emptyset()
  #define SETINSERT(a, n)            _bc_set_insert((a), (n))
  #define SETHAS(a, n)               _bc_set_has((a), (n))
#endif

#define EMPTYINTERSECTION(a,b) (!NONEMPTYINTERSECTION(a,b))
// print cluster in 0-1 repr
std::ostream &ppbitcluster(std::ostream &os, bitcluster c);

// print elements in cluster
std::ostream &ppbitclusterspecies(std::ostream &os, bitcluster c);

#endif

void clustergraphtester();	
