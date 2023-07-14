

#include <bits/stdc++.h>
using namespace std;
#include "dag.h"
#include "topsort.h"


static bool dfs(NODEID v, bool visited[], bool used[], vector<vector<NODEID>> &adj)
{
    // Mark the current node as visited.
    visited[v] = used[v] = true;
 
    // Recur for all the vertices
    // adjacent to this vertex
    vector<NODEID>::iterator i;

    for (i = adj[v].begin(); i != adj[v].end(); ++i)
    {
        // cout << "test" << v << "->" << i << endl;
        if (used[*i])
        {
            return false; // cycle
        }

        if (!visited[*i])
        {
            if (!dfs(*i, visited, used, adj))
                return false;
        }
    }
 
    used[v] = false;
    return true;  
}


bool isacyclic(vector<vector<NODEID>> &adj, NODEID V)
{
    bool visited[V];
    bool used[V];

    for (NODEID i = 0; i < V; i++)
    {
        visited[i] = used[i] = false;        
    }
     
    for (NODEID i = 0; i < V; i++)
    {
        if (visited[i] == false)
        {
            if (!dfs(i, visited, used, adj)) 
            {
                return false;
            }
        }
    }
   return true;
}

void addadjedge(vector<vector<NODEID>> &adj, NODEID v, NODEID w)
{
    adj[v].push_back(w);    
}

/*
  vector<vector<NODEID>> adj(4);
  addadjedge(adj,1,0);
  addadjedge(adj,2,1);
  addadjedge(adj,3,2);
  addadjedge(adj,3,1);
  cout << "isacyclic " << isacyclic(adj,4) << endl;
*/