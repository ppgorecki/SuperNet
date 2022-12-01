

#include <iostream>
#include <fstream>
#include <queue>


#include "bb.h"
#include "rtree.h"
#include "tools.h"
#include "network.h"
#include "contrnet.h"
#include "dag.h"
#include "dp.h"


long BBTreeStats::init(int rtnumber, COSTT start_cost)
{
	algnaivecnt = 0;
	algdpcnt = 0;
	algdptime = 0;
	algnaivetime = 0;
	minrtnumber = 100000;
	long id = v.size();
	v.push_back({rtnumber,BB_INIT,-1,start_cost,-1});
	return id;	
}

long BBTreeStats::start(int rtnumber, int algtype, long parent)
{
	long id = v.size();
	v.push_back({rtnumber,algtype,parent,-1,-1,gettime(),false});	 // MEMLEAK 3
	if (minrtnumber>rtnumber)
		minrtnumber=rtnumber;
	return id;
}
		
void BBTreeStats::stop(long id, COSTT cost)
{	
	v[id].cost = cost;	
	v[id].stime = gettime() - v[id].stime;
	if (v[id].type & ALG_NAIVE) 
	{
		algnaivecnt++;
		algnaivetime+=v[id].stime;
	}
	else if (v[id].type & ALG_DP) 
	{
		algdpcnt++;
		algdptime+=v[id].stime;
	}
}

void BBTreeStats::costcut(long id, COSTT wrt_cost)
{
	v[id].wrt_cost = wrt_cost;
	v[id].type |= BB_CUT;	
}

void BBTreeStats::bestupdated(long id, COSTT wrt_cost)
{
	v[id].wrt_cost = wrt_cost;
	v[id].type |= BB_BEST;	
}

void BBTreeStats::exactsolution(long id)
{
	v[id].type |= BB_EXACT;
}

void BBTreeStats::parentcut(long parent, COSTT costcut)
{
	long id = v.size();
	v.push_back({v[parent].rtnumber-1,BB_PARENTCUT,parent,costcut,-1});	 // MEMLEAK 4	
}

bool BBTreeStats::visitedchild(long parent)
{
	if (v[parent].visitedchild) return true;
	v[parent].visitedchild = true;
	return false;
}


void BBTreeStats::savedot()
{
	std::ofstream dotf;
	dotf.open ("bb.dot", std::ofstream::out );   
	dotf << "digraph BB {" << endl;
	for (long i = 0; i < v.size(); i++ ) 
	{
		dotf << "v" << (int)i << " ";
		// dotf << " [label=\"" << "[" << i << "]\nr=" << v[i].rtnumber; 
		// dotf << " [label=\"" << v[i].rtnumber << "R"; 
		dotf << " [label=\""; 

		if (v[i].type == BB_INIT)
		{
			dotf << "initcost=" << v[i].cost;
			dotf << "\n" << v[i].rtnumber << "R";	
			dotf << "\"]" << endl; 
		}
		else if (v[i].type & BB_PARENTCUT)
		{
			dotf << "p" << v[i].cost;	
			dotf << "\", shape=circle, style=filled, fillcolor=\"#ffe0ff\"";
			dotf << "]" << endl; 
		}
		else			
		{
			if (v[i].type & ALG_NAIVE) dotf << "\nnv=" << v[i].cost;
			if (v[i].type & ALG_DP) dotf << "" << v[i].cost;
			if (v[i].type & BB_EXACT) 
				dotf << "\nexact " << v[i].rtnumber << "R";
			// dotf << "\ntime="  << v[i].stime;

			if (v[i].type & BB_CUT) dotf << "\nc" << v[i].wrt_cost;
			if (v[i].type & BB_BEST) dotf << "\nnewbest!";

			dotf << "\"" << endl; 

			if (v[i].type & BB_BEST)
  				dotf << ", style=filled, fillcolor=\"#10ff10\"";  			
  			if (v[i].type & BB_CUT)
  				dotf << ", style=filled,fillcolor=\"#ffa0a0\"";
  			//if (v[i].type & BB_EXACT)
  			//	dotf << ", shape=box";
  			//else 
  			dotf << ", shape=box";   			
			dotf << "]" << endl; 

		}
		
		if (v[i].parent>=0)
		{
			dotf << "v" << (int)v[i].parent << " -> v" << i << endl;
			dotf << " [label=\"" << i << "\"]\n";
		}
	}
 	dotf << "}" << endl;
    dotf.close();
    // cout << "bb.dot saved" << endl;

}

void BBTreeStats::savetsv()
{
	// for (long i = 0; i < v.size(); i++ ) 
	// {

	// }
}