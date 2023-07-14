

#ifndef _TOPSORT_
#define _TOPSORT_

#include <bits/stdc++.h>
using namespace std;

#include "dag.h"

void addadjedge(vector<vector<NODEID>> &adj, NODEID v, int NODEID w);

bool isacyclic(vector<vector<NODEID>> &adj, NODEID V);


#endif