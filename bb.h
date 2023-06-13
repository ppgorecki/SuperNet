

#ifndef _BB_H_
#define _BB_H_

#if !defined USE_PRIORITY_QUEUE_MINRT && !defined USE_PRIORITY_QUEUE_MINLB && !defined USE_QUEUE_BFS
// default priority queue for BB algorithm; the best performing
#define USE_PRIORITY_QUEUE_MINLB
#endif

#include "tools.h"

#define ALG_DP 1
#define ALG_NAIVE 2


#endif
