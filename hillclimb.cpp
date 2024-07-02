#include "hillclimb.h"
#include "network.h"
#include "costs.h"
#include "netgen.h"


extern int verbosehccost;
extern int verbosealg;
extern int flag_globaldagcache;
extern int flag_hcsamplingmaxnetstonextlevel;
extern float opt_hcstoptime;	



double HillClimb::climb(
		EditOp &op, 
		Network *net, 
		CostFun &costfun, 
		ClimbStats &nhcstats, 
		bool usenaive, 
		int runnaiveleqrt_t, 		
		int hcmaximprovements,
		int hcstopclimb,
		float displaytreesampling,
		bool cutwhendtimproved,
		bool stopatcostdefined,
		float stopatcost
		)
{

	double lasthcimprovementtime = gettime();
	int timeconsistency = nhcstats.gettimeconsistency();

	if (timeconsistency!=NET_ANY)
	{
			bool tc = net->istimeconsistent();
			if (timeconsistency==NET_TIMECONSISTENT && !tc)
			{
				 cerr << "Initial network in HC is not time consistent" << endl;
				 exit(-1);
			}

			if (timeconsistency==NET_NOTIMECONSISTENT && tc)	
			{
				 cerr << "Initial network in HC is time consistent" << endl;
				 exit(-1);
			}
	}

	if (!displaytreesampling && cutwhendtimproved)
	{
		cutwhendtimproved = false;
	}

	// compute the first odt cost
	double optcost = net->odtcost(genetrees, costfun, usenaive, runnaiveleqrt_t, nhcstats.getodtstats(), displaytreesampling);

	nhcstats.addnewbest(*net, optcost); // save the first network

	int improvementscnt = 0;   // number of improvements
	int noimprovementstep = 0; // number of steps without improvement

	double curcost = optcost; 	
	std::ofstream odtf;

	bool first = true;
	bool timeout = false;


	if (!stopatcostdefined || (stopatcost<optcost))			
	while (1)
	{
		op.init(net);

		if (verbosealg>=5 || verbosehccost>=1) 
		{				
				// The first - initial network
				if (first)
				{
	      	cout << "   i: ";				
	      	cout << *net << " cost=" << optcost << endl;
	      }
	  }
	  else if (verbosehccost>=2)	
	  {
	  		// Continued HC from equal cost network
	      cout << "   e: ";
	      cout << *net << " cost=" << optcost << endl;
	  }	  

	  first = false;

	  if (opt_hcstoptime<0) break; //
	    
		while (op.next())
		{			
			

			if (opt_hcstoptime>0 && (gettime()-lasthcimprovementtime > opt_hcstoptime))
			{
				timeout = true;				
				break;
			}

			if (!net->checktimeconsistency(timeconsistency))
			{  			
				continue; // skip network
			}
			
			if (flag_globaldagcache) 
				if (nhcstats.alreadyvisited(*net))
				{								
						continue;
				}	

			if (nhcstats.inbestdags(*net))
			{
					// cout << "In best!" << endl;
					continue;
			}

			double curcost = net->odtcost(genetrees, costfun, usenaive, runnaiveleqrt_t, nhcstats.getodtstats(), displaytreesampling, cutwhendtimproved, optcost);
				 
			nhcstats.step();




			// Worse cost
			if ((verbosehccost==3) && (curcost>optcost))
			{
				cout << "   <: " << *net << " cost=" << curcost << endl;				
			}
			
			// Equal cost network
			if (curcost==optcost)
			{
				if (verbosehccost>=2)					
				{
					cout << "   =: " << *net << " cost=" << curcost << endl;				
				}
				
				nhcstats.addeq(*net);	    
			}
				
			// Yeah, new better network
			if (curcost<optcost)
			{
					optcost = curcost; 						
					if (verbosehccost>=1)
					{
						cout << "   >: " << *net << " cost=" << optcost << endl;	
					}

					// new optimal; forget old   
					nhcstats.addnewbest(*net, optcost);					

					// Stopping condition for maximprovements
					if (hcmaximprovements && (improvementscnt==hcmaximprovements))
						break;

					improvementscnt++;
				  noimprovementstep = 0; // reset counter

				  // reset 
				  lasthcimprovementtime = gettime();

					// search in a new neighbourhood
					op.reset();	

					if (stopatcostdefined)		
					{
						if (stopatcost>=optcost)
							break;
					}
			}
			else
			{
					noimprovementstep++;
			}
			

			if (hcstopclimb && (noimprovementstep>=hcstopclimb))
			{
			 	 break;
			}


		} // op.next()


		if (timeout)
		{				
				break;
		}
		// check if new equal cost networks are present

		if (nhcstats.haseqdags())
		{			
				Dag *d = nhcstats.popeqdag();				
				net = new Network(*d,false);				
		}
		else break;



	} // while (1)

	if (verbosealg>=4) 
	{

		  cout << "HC run";
		  if (displaytreesampling)
		  		cout << " with DTsampling=" << displaytreesampling;	

			if (!timeout)
      	cout << " completed: " << *net;
     	else
     		cout << " stopped due to timeout,";
     	cout << " cost=" << optcost << endl;
  }


	return optcost;

}

extern int flag_hcsamplerstats;


void iterativeretinsertionoptimizer(		
		vector<RootedTree*> &gtvec,				
		Network *startnet,		
		CostFun *costfun,		
		int printstats,						
		bool usenaive, 
		int runnaiveleqrt_t, 				
		ClimbStatsGlobal *globalstats, // could be sampler
		bool cutwhendtimproved,
		int networkclass, 
		int timeconsistency, 
		Clusters *guideclusters, 
		Clusters *guidetree, 
		string _retid
		)

{
 	  DagSet visiteddags;		
		ClimbStats climbstats(visiteddags, globalstats);  
		float displaytreesampling = globalstats->getsampling();

	  double optcost = startnet->odtcost(gtvec, *costfun, usenaive, runnaiveleqrt_t, climbstats.getodtstats(), displaytreesampling);

	  climbstats.addnewbest(*startnet, optcost); // save the first network

		NetworkRetIterator netretit(*startnet, networkclass, timeconsistency, guideclusters, guidetree, "");
   	Network *net;
   	   	
    while ((net = netretit.next())!=NULL)
    {
      	double curcost = net->odtcost(gtvec, *costfun, usenaive, runnaiveleqrt_t, climbstats.getodtstats(), displaytreesampling);

				// Equal cost network
				if (curcost==optcost)
				{
					if (verbosehccost>=2)					
					{
						cout << "   =: " << *net << " cost=" << curcost << endl;				
					}
					
					climbstats.addeq(*net);	    
				}
				
			// Yeah, new better network
			if (curcost<optcost)
			{
					optcost = curcost; 						
					if (verbosehccost>=1)
					{
						cout << "   >: " << *net << " cost=" << optcost << endl;	
					}

					// new optimal; forget old   
					climbstats.addnewbest(*net, optcost);					

					
					// improvementscnt++;
				  // noimprovementstep = 0; // reset counter

				  // reset 
				  // lasthcimprovementtime = gettime();
					
			}


      	// cout << cost << " " << *n << endl;
    }

		climbstats.finalize();

		globalstats->merge(climbstats, printstats, false);
}

void supnetheuristic(		
		vector<RootedTree*> &gtvec,				
		NetGen *netgenerator,
		EditOp *op,
		CostFun *costfun,		
		int printstats,				
		int hcstopinit,
		int hcstopclimb,
		bool usenaive, 
		int runnaiveleqrt_t, 		
		int hcmaximprovements,		
		vector<ClimbStatsGlobal*> globalstatsarr,
		bool cutwhendtimproved,
		int multipleoptima,
		bool stopatcostdefined,
		float stopatcost,
		DagScoredQueue *scoreddags
		)
{
      
		HillClimb hc(gtvec);    
    DagSet visiteddags;

    vector<pair<Network *, int>> samplerarr;

    long int hccnt = -1;
    int lastimprovement = 0;

    int cnt = 0;
    while (1) 
    {
      // stopping criterion
      if (hcstopinit && (hccnt - lastimprovement) > hcstopinit)
        break; // stop

      // get next initial network
      Network *n = netgenerator->next();
      
      if (!n)
      {
        if (!cnt)
        {
          cerr << "No network defined. Use -q, -r, --guidetree or --guideclusters to generate some starting networks." << endl;
          exit(-1);
        }
        break;
      }

      cnt++;
      hccnt++;

      // preprocess using samplers

      unsigned int maxdt = n->displaytreemaxid();

      samplerarr.push_back(make_pair(n, 0));
      
			while (!samplerarr.empty())
			{

					pair<Network*, int> p = samplerarr.back();					
					samplerarr.pop_back();					


					DagSet visiteddags;				
					ClimbStatsGlobal *globalstats = globalstatsarr[p.second];

					if (flag_hcsamplerstats)
					{
						cout << "Level " <<  p.second << " ";
						globalstats->info(cout) <<  endl;
					}

					ClimbStats nhcstats(visiteddags, globalstats);      
					nhcstats.start();	

					double cost = hc.climb(*op, p.first, *costfun, nhcstats, 
														 usenaive,
                             runnaiveleqrt_t, 
                             hcmaximprovements, 
                             hcstopclimb,
                             globalstats->getsampling(),
                             cutwhendtimproved,
                             stopatcostdefined,
														 stopatcost);



					if (scoreddags && !globalstats->issampler())
					{
							
						  // insert best from 						  
						  int cnt=multipleoptima;
						  for (auto a : *nhcstats.getbestdags())
						  {
						  	 cnt--;
								 scoreddags->push(a, cost);
								 if (!cnt) break;
							}
					}

					if (globalstats->issampler())
					{
						// all best dags to the next level sampler

						DagSet *dagset = nhcstats.getbestdags();
						int nextsampler = p.second + 1;

						int ncnt = 0;
						//for (std::vector<Dag*>::iterator a = dagset->begin(); a!=dagset->end(); a++ )			
						for (auto a : *dagset)
						{													
								if ((flag_hcsamplingmaxnetstonextlevel>0) && ncnt == flag_hcsamplingmaxnetstonextlevel)
									break;
								ncnt++;

								Network *net = new Network(*a, false);			
								//cout << *net << endl;
								samplerarr.push_back(make_pair(net, nextsampler));
						}

						if (flag_hcsamplerstats)
						{
							cout << "" << ncnt << " (out of " << dagset->size() << ") network(s) moved to next level: ";
							if (globalstatsarr.size() == nextsampler+1 )
								cout << " exact climb" << endl;
							else 
							 	cout << " sampler " << nextsampler <<  endl;
						}
					}

					nhcstats.finalize();

					if (!globalstats->issampler())
					{
						 hccnt++; // Double hccnt++ TODO - check
					}

					if (globalstats->merge(nhcstats, printstats, false)) 
      		{
      				if (!globalstats->issampler())
        				lastimprovement = hccnt;
      		}
			}      

			delete n;
      
      
      
    } 

    
}




