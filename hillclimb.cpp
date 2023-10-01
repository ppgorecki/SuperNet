#include "hillclimb.h"
#include "network.h"
#include "costs.h"
#include "randnets.h"

void EditOp::init(Network *net) 
{ 
	source = net; 
	reset();
}

extern int verbosehccost;
extern int verbosealg;
extern int flag_globaldagcache;
extern int flag_hcsamplingmaxnetstonextlevel;
extern float opt_hcstoptime;

//#define TAILMOVE_DEBUG

void TailMove::reset()
{
	u = source->lf-1;
	v = 0;
	if (!vreachable)
		vreachable = new bool[source->size()];	

	s = source->size();	
	t = 0;
	moved = false;
}


// #defind TreeNode_ChildSiblingPTR(child,parent) 

// Get parent address
// child - par must be an edge
#define ParentAddr(child,par) ((par)==source->parent[(child)])?(source->parent+(child)):(source->retparent+(child))

// Get child address
// child - par must be an edge
// Assumption retchild == leftchild
#define ChildAddr(child,par) ((child)==source->leftchild[(par)])?(source->leftchild+(par)):(source->rightchild+(par))


// Get sibling address of a child 
#define TreeNodeSiblingChildAddr(child,par) (source->leftchild[(par)]==(child))?(source->rightchild+(par)):(source->leftchild+(par))	

bool TailMove::next() 
{ 
	// TODO: uniform random sampling; store all candidates and uniformly draw

#ifdef TAILMOVE_DEBUG	
	int tmcnt = 0 ;
#endif 

	if (moved)
	{
		//perform reverse move
		move(u,v,p,q);
		moved = false;		
	}

	while (1)
	{	

		// Find four candidate nodes 
		do  
		{	
			
			if (s >= source->size())
			{
				// gen new u-t edge
				if ((u<source->lf) || (v == source->rightchild[u]))
				{
					// new u 
					if (++u == source->rtstartid) return false; // end of neighbourhood			
					v = source->leftchild[u];
				}
				else	
				{
					// try next edge with the same parent u
					v = source->rightchild[u];
				}

				// gen v-reachable 		
				source->getreachablefrom(v, vreachable);

				// initialize s and search for new s			
				s = source->lf;			
				
			}
			else
			{
				if (!vreachable[s] && s<source->rtstartid && source->leftchild[s]==t) 
				{ 
					t=source->rightchild[s]; 
					break; // old s (tree node), new t
				}
				else			
					s++; // search for new s-t		
			}

			// Search for s among non-v-reachable nodes
			// (if s is v-reachable then kids of s are also v-reachable)		

			for (; s<source->size(); s++)
				if (!vreachable[s]) 
				{
					if (s<source->rtstartid)
						t = source->leftchild[s];
					else
						t = source->retchild[s];
					break; 	
				}


			
		} while (s==source->size());  // do while - search for candidate nodes

	  // Candidates:
	  // 	u-v fixed 
	  // 	s-t fixed

		p = source->parent[u];
		q = *(TreeNodeSiblingChildAddr(v,u));

#ifdef TAILMOVE_DEBUG	

		cout << ++tmcnt << ". Candidate: p=" << p << " q=" << q << " || " <<
			" u=" << u << " v=" << v << " | s=" << s << " t=" << t << "  ";
#endif
	
		// Check conditions

		if (v>=source->rtstartid && (source->parent[v]==q || source->retparent[v]==q) && s==q) 
		{ 
#ifdef TAILMOVE_DEBUG					
			cout << "2b2" << endl;
#endif			
			continue;
		}
		
		if (v==t) 
		{ 
#ifdef TAILMOVE_DEBUG								
			cout << "1a" << endl;
#endif			
			continue;
		}


		if (s==u && t==q) 
		{ 
#ifdef TAILMOVE_DEBUG								
			cout << "2a2" << endl;
#endif						
			continue;
		}

		if (s==p && t==u) 
		{ 
#ifdef TAILMOVE_DEBUG								
			cout << "2a1" << endl;
#endif						
			continue;
		}

		if (q>=source->rtstartid && (source->parent[q]==p || source->retparent[q]==p)) 
		{ 
#ifdef TAILMOVE_DEBUG								
			cout << "1b" << endl;
#endif						
			continue;
		}


		if (v>=source->rtstartid && (source->parent[v]==p || source->retparent[v]==p) && t==p) 
		{ 
#ifdef TAILMOVE_DEBUG								
			cout << "2b1" << endl;
#endif			

			continue;
		}

#ifdef TAILMOVE_DEBUG	
		cout << endl;
#endif

		if (networkclass==NET_TREECHILD || networkclass==NET_CLASS1RELAXED)
		{
			NODEID rtstartid = source->rtstartid;

			// tree child conditions 3a-3c
			if (v>=rtstartid && t>=rtstartid)
			{
#ifdef TAILMOVE_DEBUG								
				cout << "3a" << endl;
#endif			
				continue;
			}

			if ((networkclass==NET_TREECHILD) && (p>=rtstartid && q>=rtstartid))
			{
#ifdef TAILMOVE_DEBUG								
				cout << "3b" << endl;
#endif			
				continue;
			}

			if (q>=source->rtstartid) 
			{
				if (p >= source->lf)  // p - is a non-leaf tree node
				{				
					NODEID pr = MAXNODEID;
					NODEID rc = source->rightchild[p]; 
					NODEID lc = source->leftchild[p];

					if ((lc >= rtstartid) && (rc >= rtstartid))
					{
						// network is not tree-child

						cerr << "Tree-child network expected" << endl;
						exit(-1);

						continue;
					}

					if (lc >= rtstartid) pr = lc;
					if (rc >= rtstartid) pr = rc;

					if (pr != MAXNODEID) // ret. child is found
					{
						if ((s!=p) || (t!=pr))
						{

#ifdef TAILMOVE_DEBUG								
							cout << "3c" << endl;
#endif	
							continue;
						}
					}
				}				
			} // 3c

		} // limittotreechild or relaxed

		// move!
		move(u,v,s,t);

		// check guideclusters condition		
		if (guideclusters)
		{
				
				if (!source->hasclusters(guideclusters))
				{											
#ifdef TAILMOVE_DEBUG			
						cout << tmcnt << " guideclusters_rejected:" << *source << ":::";
						source->printclusters(cout);
						cout << endl;
#endif						
						//move(u,v,s,t); // revert						
						move(u,v,p,q);

#ifdef TAILMOVE_DEBUG			
						cout << tmcnt <<  "Reverted" << endl;
#endif
						continue;
				}
		}

		// check guidetree condition		
		if (guidetree)
		{
				
				if (!source->hastreeclusters(guidetree))
				{											
#ifdef TAILMOVE_DEBUG			
						cout << tmcnt << " guidetree_rejected:" << *source << ":::";
						// source->printclusters(cout);
						cout << tmcnt << endl;
#endif						
						//move(u,v,s,t); // revert	
						move(u,v,p,q);					
#ifdef TAILMOVE_DEBUG			
						cout << "Reverted" << endl;
#endif						
						continue;
				}
		}
		
#ifdef TAILMOVE_DEBUG			
		cout << tmcnt << "ACCEPT" << *source << endl;
#endif		
		
		moved = true;

#ifdef TAILMOVE_DEBUG		
		cout << endl << tmcnt << "VER:" << source->verifychildparent() << endl;;
		source->printdeb(cout,2) << endl;
#endif		

		// Done
		return true;
	}
}

void TailMove::move(NODEID un, NODEID vn, NODEID sn, NODEID tn)
{

	NODEID pn = source->parent[un];
	// NODEID qn = source->sibling(vn); - wrong!
	NODEID qn = *(TreeNodeSiblingChildAddr(vn,un));

	NODEID *qpar = ParentAddr(qn,un);
	*qpar = pn; 

#ifdef TAILMOVE_DEBUG
	cout << " In move: moved=" << moved << " TAILMOVE: p=" << pn << " q=" << qn << " || " <<
			" u=" << un << " v=" << vn << " | s=" << sn << " t=" << tn << "  " << endl;
#endif

	// move!
	// disconnect u; connect p -> q
	if (pn==MAXNODEID) // root	
		source->root = qn;			
	else
	{
		NODEID *pchild = ChildAddr(un,pn);		
		*pchild = qn;
	}
		
	source->parent[u] = sn;
	if (sn==MAXNODEID)
	{
		//above root move
		source->root = un; // new root
		source->parent[tn] = un;
	}
	else
	{
		NODEID *schild = ChildAddr(tn,sn);
		NODEID *tpar = ParentAddr(tn,sn);

		// attach s->u->t
		*schild = un;
		*tpar = un;
	}

	// reattach u->q => u->t
	NODEID *qchild = ChildAddr(qn,un);
	*qchild = tn;
		
}

void NNI::reset() 
{ 
	curnode = source->lf-1; // search from the next
	cycle = 0;
}

bool NNI::next() 
{ 
	if (!cycle)
	{	
		curnode++;
		NODEID p;
		for (; curnode<source->rtstartid; curnode++)
			if (source->parent[curnode]<source->rtstartid) 
			{
				// Set adresses of children (subtrees to be rotated)
				ac = source->leftchild+curnode;
				bc = source->rightchild+curnode;
				p = source->parent[curnode];

				cc = TreeNodeSiblingChildAddr(curnode,p);
				// if (source->leftchild[p]==curnode) 
				// 	cc=source->rightchild+p;
				// else cc=source->leftchild+p;		

				if (*ac==*bc || *ac==*cc || *bc==*cc)				
					// avoid possible double edges
					continue;
				

				break; // found
			}

		if (curnode==source->rtstartid) 
			return false; // no more candidates

		// Set adresses of parents (from the subtrees)
		ap = ParentAddr(*ac,curnode);
		bp = ParentAddr(*bc,curnode);
		cp = ParentAddr(*cc,p);

		cycle=3;
	}

	// cout << " NNI node="<< (int) curnode 
	// 			<< " cycle=" << cycle << 
	// 		" " << (*ac) << " " << (*bc) << " " << (*cc) << " | " <<
	// 		" " << (*ap) << " " << (*bp) << " " << (*cp) << endl;

	// rotate children
	NODEID d = *ac;
	*ac = *bc;
	*bc = *cc;
	*cc = d;

	// rotate parent addresses
	NODEID *p = ap;
	ap = bp;
	bp = cp;
	cp = p; 

	// assign parent links
	*ap = curnode;
	*bp = curnode;
	*cp = source->parent[curnode];
	
	cycle--;

	if (!cycle)	
		// Rotations of the current edge are completed 
		// and the same network is reached 
		// Search for the next 
		return next();

	return true; 
}


double HillClimb::climb(
		EditOp &op, 
		Network *net, 
		CostFun &costfun, 
		NetworkHCStats &nhcstats, 
		bool usenaive, 
		int runnaiveleqrt_t, 		
		int hcmaximprovements,
		int hcstopclimb,
		float displaytreesampling,
		bool cutwhendtimproved
		)
{

	double lasthcimprovementtime = gettime();
	int timeconsistent = nhcstats.gettimeconsistency();
	

	if (timeconsistent)
	{
			bool tc = net->istimeconsistent();
			if (timeconsistent==TIMECONSISTENT && !tc)
			{
				 cerr << "Initial network in HC is not time consistent" << endl;
				 exit(-1);
			}

			if (timeconsistent==NOTIMECONSISTENT && tc)	
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

			if (timeconsistent)			
			{  
				bool tc = net->istimeconsistent();
				if (timeconsistent==TIMECONSISTENT && !tc) continue; // skip network
				if (timeconsistent==NOTIMECONSISTENT && tc) continue; // skip network
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
			}
			else
			{
					noimprovementstep++;
			}

			if (hcstopclimb && (noimprovementstep==hcstopclimb))
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
			if (!timeout)
      	cout << "HC run completed: " << *net << " cost=" << optcost << endl;
     	else
     		cout << "HC run stopped due to timeout, cost=" << optcost << endl;
  }


	return optcost;

}

extern int flag_hcsamplerstats;

void supnetheuristic(		
		vector<RootedTree*> &gtvec,				
		NetIterator *netiterator,
		EditOp *op,
		CostFun *costfun,		
		int printstats,				
		int hcstopinit,
		int hcstopclimb,
		bool usenaive, 
		int runnaiveleqrt_t, 		
		int hcmaximprovements,		
		vector<NetworkHCStatsGlobal*> globalstatsarr,
		bool cutwhendtimproved
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
      Network *n = netiterator->next();
      
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
					NetworkHCStatsGlobal *globalstats = globalstatsarr[p.second];

					if (flag_hcsamplerstats)
					{
						cout << "Level " <<  p.second << " ";
						globalstats->info(cout) <<  endl;
					}

					NetworkHCStats nhcstats(visiteddags, globalstats);      
					nhcstats.start();					

					double cost = hc.climb(*op, p.first, *costfun, nhcstats, 
														 usenaive,
                             runnaiveleqrt_t, 
                             hcmaximprovements, 
                             hcstopclimb,
                             globalstats->getsampling(),
                             cutwhendtimproved);

					if (globalstats->issampler())
					{
						// all best dags to the next level sampler

						DagSet *dagset = nhcstats.getbestdags();
						int nextsampler = p.second + 1;

						int ncnt = 0;
						for (std::vector<Dag*>::iterator a = dagset->begin(); a!=dagset->end(); a++ )			
						{													
								if ((flag_hcsamplingmaxnetstonextlevel>0) && ncnt == flag_hcsamplingmaxnetstonextlevel)
									break;
								ncnt++;

								Network *net = new Network(**a, false);			
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
						 hccnt++;
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