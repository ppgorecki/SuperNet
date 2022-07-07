

#include <queue>

#include "bb.h"
#include "rtree.h"
#include "tools.h"
#include "network.h"
#include "contrnet.h"
#include "dag.h"
#include "dp.h"

// DC via BB
// #define _DEBUG_DPBB_
COSTT Network::mindc(RootedTree &genetree, int runnaiveleqrt, AdaptiveBB *adaptivebb)
{      
	
	RootedTree *t = gendisplaytree(0,NULL);
	if (!t)
	{
	   cerr << "No display tree generated" << endl;
	   exit(-2);
	}

	COSTT best_cost = genetree.cost(*t,COSTDEEPCOAL)+genetree.lf*2-2; // convert to non-classic 

#ifdef _DEBUG_DPBB_
	cout << " displaytree-startcost:" << best_cost << endl;
#endif	
	delete t;

    long branching_count = 0;
    int max_depth = 0;
    int best_depth = 0;
    long dp_called = 0;
    long naive_called = 0;

    if (!adaptivebb) adaptivebb = new AdaptiveBB();

    long bbnodeid = adaptivebb->init(rtcount(), best_cost);

    typedef struct {
		ContractedNetwork *src;
		SPID rtid;
		bool left;
		long bbparnodeid;
        COSTT parcost;
    } QData;

    queue<QData> q;
   
    q.push({ .src=NULL, .rtid=this->rt, .left=false, .bbparnodeid = bbnodeid, .parcost=0}); // init

    while (!q.empty())
    {
    	QData qdata = q.front();    	
   		q.pop();        

   		ContractedNetwork *srcc = qdata.src;
   		SPID rtid = qdata.rtid;
   		bool left = qdata.left;
   		long bbparnodeid = qdata.bbparnodeid;

        // New optimization
        if (qdata.parcost>=best_cost)
        {
          // parent lower bound is worst than the current best cost 
          // this branch cannot improve the cost
          if (bbparnodeid>=0)
            adaptivebb->parentcut(bbparnodeid,best_cost);
#ifdef _DEBUG_DPBB_
        cout << "Parent-Cut=" << qdata.parcost << " vs " << best_cost << endl;
#endif      
          continue;

        }

   		//cout << "SRC:" << srcc << " " << left << endl;

   		RETUSAGE retusage;
   		ContractedNetwork *c;

   		if (!srcc)
   		{   			   		
   			c = new ContractedNetwork(*this, true);
   		}  		
   		else
   		{
   			// build a new contracted network
   			c = new ContractedNetwork(*srcc, true);
        	RETUSAGE retusage;
        	emptyretusage(retusage);
        	if (left)
        		addleftretusage(retusage,rtid);
        	else addrightretusage(retusage,rtid);
        	c->contract(retusage);        
   		}	

   		// run DP

		COSTT cost;
		bool naivecomputed = false;
		int rtnum = c->rtcount();

   		if (runnaiveleqrt>rtnum)   		
   		{
   			      
        
   			bbnodeid = adaptivebb->start(rtnum, ALG_NAIVE, bbparnodeid);

   			// compute naive (exact)
   			cost =  c->odtcostnaive(&genetree, COSTDEEPCOAL) + genetree.lf*2 - 2;  
   			adaptivebb->stop(bbnodeid, cost);		
   			naive_called++;
   			naivecomputed = true;

#ifdef _DEBUG_DPBB_
        cout << "Naive-cost=" << cost << endl;
#endif
   		}
   		else 
   		{ 
   			bbnodeid = adaptivebb->start(rtnum, ALG_DP, bbparnodeid);

   			// compute via DP (lower bound)
   			cost = c->approxmindcusage(genetree, retusage);
   			adaptivebb->stop(bbnodeid, cost);
   			dp_called++;
#ifdef _DEBUG_DPBB_
        cout << "DP-cost=" << cost << endl;
#endif
   		}
   		
   		// If a cost for a potentially incomplete tree is higher than the one we achieved,
        // we can stop at that point
        if (cost >= best_cost)
        {
#ifdef _DEBUG_DPBB_        	
        	cout << " Cut branch" << endl;
#endif

        	adaptivebb->costcut(bbnodeid, best_cost);

        	if (left) delete srcc;
        	delete c;
       		continue;
       	}

       	// Solution has all reticulation parents set, we can use the cost
#ifdef _DEBUG_DPBB_        	       	
       	cout << " Retusage " << retusage << endl;
#endif       	

       	// if computed exactly - store and terminate branch
     	if (naivecomputed || !conflicted(retusage)) 
     	{
            
#ifdef _DEBUG_DPBB_                     
            cout << " No-conflict: " << cost << endl;
#endif          

            if (cost < best_cost)
            {
#ifdef _DEBUG_DPBB_                     
            cout << " Best-cost update" << endl;
#endif          
                best_cost = cost;            
                adaptivebb->bestupdated(bbnodeid, best_cost); 
            }
            adaptivebb->exactsolution(bbnodeid); 

            if (left)  
            	delete srcc;    


            delete c;   
            continue;
        }    


        // Conflicted solution from DP

        // Find conflicted reticulation
        rtid = c->getconflictedreticulation(retusage);

        // Branching
        adaptivebb->branch(bbnodeid, rtid); 


#ifdef _DEBUG_DPBB_ 
        cout << "CR:" << rtid << " " << spid2retlabel[rtid+rtstartid] << endl;
#endif        
 		
 		// Insert contracted networks 
 		q.push({.src=c, .rtid=rtid, .left=false, .bbparnodeid = bbnodeid, .parcost=cost}); 
 		q.push({.src=c, .rtid=rtid, .left=true, .bbparnodeid = bbnodeid, .parcost=cost}); 

        if (left) delete srcc;

    }
	
    return best_cost - genetree.lf*2 + 2; 

}