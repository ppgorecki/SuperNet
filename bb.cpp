#include <queue>

#include "bb.h"
#include "rtree.h"
#include "tools.h"
#include "costs.h"
#include "network.h"
#include "contrnet.h"
#include "dag.h"
#include "dp.h"

extern int verbosealg;

// DCE via BB with naive switching
// #define _DEBUG_DPBB_

COSTT Network::mindce(
    RootedTree &genetree, 
    int runnaiveleqrt_t, 
    CostFun &costfun,
    ODTStats &odtstats,
    BBTreeStats *bbtreestats,
    COSTT bbstartscore,
    bool bbstartscoredefined)
{      
	
    COSTT best_cost = 0;

    COSTT dc2dce = 0;

    if (costfun.costtype()==COSTDEEPCOAL) dc2dce = 2*genetree.lf-2; // adjust DC to DCE


    if (bbstartscoredefined)
    {
#ifdef _DEBUG_DPBB_
        cout <<"#INITSCORE " << bbstartscore << endl;
#endif        
        best_cost = bbstartscore;
    }
    else 
    {

        // initialize cost value
    	RootedTree *t = gendisplaytree(0,NULL);
    	if (!t)
    	{
    	   cerr << "No display tree generated" << endl;
    	   exit(-2);
    	}

    	best_cost = costfun.computegt(genetree, *t)+dc2dce; 
        // convert to     non-classic 
        delete t;
    }

#ifdef _DEBUG_DPBB_
	cout << " displaytree-startcost:" << best_cost << endl;
#endif	

    if (verbosealg>=6) 
    { 
        cout << " BB gid=" << genetree.getid() << " startcost=" << best_cost << endl;
    }

    long branching_count = 0;
    int max_depth = 0;
    int best_depth = 0;

    bool bbstatsallocated = false;

    if (!bbtreestats) { 
        bbtreestats = new BBTreeStats(); // MEMLEAK 1
        bbstatsallocated = true;
    }

    long bbnodeid = bbtreestats->init(rtcount(), best_cost);

    typedef struct {
        long nodeid;
		ContractedNetwork *src;
		NODEID rtid;
		bool left;
		long bbparnodeid;
        COSTT parcost;
        int rtnumber;        
    } QData;

#ifdef USE_PRIORITY_QUEUE_MINRT    

    // Min RT network - first
    struct CompQData 
    {
        bool operator()(QData const& p1, QData const& p2)
        {        
            if ((p1.bbparnodeid == p2.bbparnodeid) && !p2.left) 
                return true;    
            return p1.rtnumber > p2.rtnumber;
        };
    };

    priority_queue<QData, vector<QData>, CompQData> q;


#elif defined USE_PRIORITY_QUEUE_MINLB

    // Min lower bound RT network - first
    struct CompQData 
    {
        bool operator()(QData const& p1, QData const& p2)
        {        
            if ((p1.bbparnodeid == p2.bbparnodeid) && !p2.left) 
                return true;    
            return p1.parcost > p2.parcost;
        };
    };

    priority_queue<QData, vector<QData>, CompQData> q;

#elif defined USE_QUEUE_BFS
    // Classic BFS 
    queue<QData> q;
#else
    // Yes, error
#endif

    q.push({ .nodeid=bbnodeid, .src=NULL, .rtid=this->rt, .left=false, 
        .bbparnodeid = bbnodeid, .parcost=0, 
        .rtnumber=rtcount()}); // init

    while (!q.empty())
    {
        //cout <<"#score " << best_cost << endl;
#if defined USE_PRIORITY_QUEUE_MINRT || defined USE_PRIORITY_QUEUE_MINLB

        QData qdata = q.top();                

#if DEBUG_BB_QUEUES
        // print priority queue data
        priority_queue<QData, vector<QData>, CompQData> _q = q;    
        cout << "#[ ";
        for(int i = 0; i < q.size(); ++i) {
            cout << _q.top().parcost  << " ";
            _q.pop();
        }
        cout << "] " << endl;

#endif 

#else    
    	QData qdata = q.front();    	
#endif

   		q.pop();  

        // printf("#bestscore=%ld startscore=%ld %d\n",best_cost, bbstartscore,bbstartscoredefined);

   		ContractedNetwork *srcc = qdata.src;
   		NODEID rtid = qdata.rtid;
   		bool left = qdata.left;
   		long bbparnodeid = qdata.bbparnodeid;

        // cout << qdata.nodeid << " lft=" << left << " par=" << bbparnodeid << " rt=" << qdata.rtnumber << endl;      

        if (verbosealg>=6)
        {
            cout << " BB-main-loop parcost=" << qdata.parcost << " best=" << best_cost << endl; 
        }

        // New optimization
        if (qdata.parcost>=best_cost)
        {
          // parent lower bound is worst than the current best cost 
          // this branch cannot improve the cost
          if (bbparnodeid>=0)
          {

            bbtreestats->parentcut(bbparnodeid,best_cost);
          }

          if (verbosealg>=6)
          {
                cout << " --parent-cut=" << qdata.parcost << " vs " << best_cost << endl;
          }

          continue;
        }

   		RETUSAGE retusage;
   		ContractedNetwork *c;

   		if (!srcc)
   		{   			   		
   			c = new ContractedNetwork(*this, true);
   		}  		
   		else
   		{
   			// build a new contracted network
   			c = new ContractedNetwork(*srcc, true); // MEMLEAK 2
        	RETUSAGE retusage;
        	emptyretusage(retusage);
        	if (left)
        		addleftretusage(retusage,rtid);
        	else addrightretusage(retusage,rtid);
        	c->contract(retusage);        
   		}	

// #define _VISNETS_


   		// run DP


		COSTT cost;
		bool naivecomputed = false;
		int rtnum = c->rtcount();

   		if (runnaiveleqrt_t>rtnum)   		
   		{
   			                  
            

   			bbnodeid = bbtreestats->start(rtnum, ALG_NAIVE, bbparnodeid);

   			// compute naive (exact)

   			cost =  c->odtcostnaive(&genetree, costfun, odtstats) + dc2dce;  

            if (verbosealg>=6) 
            { 
                cout << " --naive=" << cost << " rt=" << rtnum <<  endl;
            }

   			bbtreestats->stop(bbnodeid, cost);		   			
   			naivecomputed = true;

   		}
   		else 
   		{ 
            
   			bbnodeid = bbtreestats->start(rtnum, ALG_DP, bbparnodeid);

   			// compute via DP (lower bound)
   			cost = c->approxmindceusage(genetree, retusage, costfun);
   			bbtreestats->stop(bbnodeid, cost);

            if (verbosealg>=6) 
            { 
                cout << " --dp=" << cost << endl;
            }
   			

#ifdef _VISNETS_
        //cout << "BBDP[" << genetree.gettid() << "] " << *c << endl;
        cout << "BBDP[" << genetree.>getid() << "] ";
        c->printrepr(cout) << endl;
#endif        

   		}
   		
   		// If a cost for a potentially incomplete tree is higher than the one we achieved,
        // we can stop at that point
        if (cost >= best_cost)
        {

            if (verbosealg>=6) 
            { 
                cout << " --cut branch " 
                    << cost << ">=" << best_cost
                    << endl;
            }
            

        	bbtreestats->costcut(bbnodeid, best_cost);

        	if (bbtreestats->visitedchild(bbparnodeid)) delete srcc;            
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
        
            if (cost < best_cost)
            {
                if (verbosealg>=6)
                    cout << " --no-conflict & best-cost update" << cost << endl;

                best_cost = cost;            
                bbtreestats->bestupdated(bbnodeid, best_cost); 

            }            
            bbtreestats->exactsolution(bbnodeid); 
            
            if (bbtreestats->visitedchild(bbparnodeid)) delete srcc;            

            delete c;   
            continue;
        }    


        // Conflicted solution from DP

        // Find conflicted reticulation
        rtid = c->getconflictedreticulation(retusage);

        
#ifdef _DEBUG_DPBB_ 
        cout << "CR:" << rtid << " " << spid2retlabel[rtid+rtstartid] << endl;
#endif        
 		
 		// Insert contracted networks 
 		q.push({.nodeid=bbnodeid, .src=c, .rtid=rtid, .left=false, .bbparnodeid = bbnodeid, .parcost=cost, .rtnumber=rtnum-1 }); 
 		q.push({.nodeid=bbnodeid, .src=c, .rtid=rtid, .left=true, .bbparnodeid = bbnodeid, .parcost=cost, .rtnumber=rtnum-1}); 

        if (verbosealg>=6)
            cout << " --conflicted" << cost << endl;

        // if (left) delete srcc;
        if (bbtreestats->visitedchild(bbparnodeid)) delete srcc;            
    }


    odtstats.addbbstats(bbtreestats->stats);

    if (bbstatsallocated)
        delete bbtreestats;
	
    return best_cost; 

}