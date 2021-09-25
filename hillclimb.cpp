#include "hillclimb.h"
#include "network.h"

#define HC_DEBUG

void EditOp::init(Network *net) 
{ 
	source = net; 
	reset();
}


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


	if (moved)
	{
		//peform reverse move
		move(u,v,p,q);
		moved = false;
	}

	while (1)
	{	
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


			
		} while (s==source->size());  // do while

	// Candidates:
	// 	u-v fixed 
	// 	s-t fixed

		p = source->parent[u];
		q = *(TreeNodeSiblingChildAddr(v,u));

#ifdef TAILMOVE_DEBUG	
		cout << " p=" << p << " q=" << q << " || " <<
			" u=" << u << " v=" << v << " | s=" << s << " t=" << t << "  ";
#endif
		// Now check conditions

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

		if (limittotreechild)
		{
			SPID rtstartid = source->rtstartid;

			// tree child conditions 3a-3c
			if (v>=rtstartid && t>=rtstartid)
			{
#ifdef TAILMOVE_DEBUG								
				cout << "3a" << endl;
#endif			
				continue;
			}

			if (p>=rtstartid && q>=rtstartid)
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
					SPID pr = MAXSP;
					SPID rc = source->rightchild[p];
					SPID lc = source->leftchild[p];

					if ((lc >= rtstartid) && (rc >= rtstartid))
					{
						// network is not tree-child

						cerr << "Tree-child network expected" << endl;
						exit(-1);

						continue;
					}

					if (lc >= rtstartid) pr = lc;
					if (rc >= rtstartid) pr = rc;

					if (pr != MAXSP) // ret. child is found
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

		} // limittotreechild

		// Finally, move!
		move(u,v,s,t);
		moved = true;

#ifdef TAILMOVE_DEBUG		
		cout << endl << "VER:" << source->verifychildparent() << endl;;
		source->printdeb(cout,2) << endl;
#endif		

		// Done
		return true;
	}
}

void TailMove::move(SPID un, SPID vn, SPID sn, SPID tn)
{

	SPID pn = source->parent[un];
	// SPID qn = source->sibling(vn); - wrong!
	SPID qn = *(TreeNodeSiblingChildAddr(vn,un));

	SPID *qpar = ParentAddr(qn,un);
	*qpar = pn; 

#ifdef TAILMOVE_DEBUG
	cout << dagcnt << " moved=" << moved << " TAILMOVE: p=" << pn << " q=" << qn << " || " <<
			" u=" << un << " v=" << vn << " | s=" << sn << " t=" << tn << "  ";
#endif

	// move!
	// disconnect u; connect p -> q
	if (pn==MAXSP) // root	
		source->root = qn;			
	else
	{
		SPID *pchild = ChildAddr(un,pn);		
		*pchild = qn;
	}
		
	source->parent[u] = sn;
	if (sn==MAXSP)
	{
		//above root move
		source->root = un; // new root
		source->parent[tn] = un;
	}
	else
	{
		SPID *schild = ChildAddr(tn,sn);
		SPID *tpar = ParentAddr(tn,sn);

		// attach s->u->t
		*schild = un;
		*tpar = un;
	}

	// reattach u->q => u->t
	SPID *qchild = ChildAddr(qn,un);
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
		SPID p;
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
	SPID d = *ac;
	*ac = *bc;
	*bc = *cc;
	*cc = d;

	// rotate parent addresses
	SPID *p = ap;
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


double HillClimb::climb(EditOp &op, Network *net, int costfunc)
{

	double starttime = gettime();

	// init edit operation
	op.init(net);

	// compute the first odt cost
	double optcost = net->odtnaivecost(genetrees, costfunc);
	
	double curcost; 

	long steps = 0;
	long improvements = 0;
	long equalnets = 0;
	std::ofstream odtf;

	bool saveodtfile = odtfile.length();

	if (saveodtfile) 
	{
		odtf.open (odtfile, std::ofstream::out);
		odtf << *net << endl; // save the first network
	}
	if ((verbose==3) && (curcost>optcost))
		cout << " = " << *net << " cost=" << curcost << endl;				
    

	while (op.next())
	{		
		
		double curcost = net->odtnaivecost(genetrees, costfunc);
		steps++;

		if ((verbose==3) && (curcost>optcost))
			cout << " < " << *net << " cost=" << curcost << endl;				
		
		if (curcost==optcost)
		{
			equalnets++;
			if (verbose>=2)					
				cout << " = " << *net << " cost=" << curcost << endl;				
			if (saveodtfile) 
				odtf << *net << endl;
    
		}
			
		// Yeah, new better network
		if (curcost<optcost)
		{
			optcost = curcost; 
			
			equalnets = 1;
			improvements++;

			if (verbose>=1)
				cout << " > " << *net << " cost=" << optcost << endl;	

			if (saveodtfile) 
			{
				odtf.close();
				odtf.open (odtfile, std::ofstream::out);
				odtf << *net << endl;		
			}

			// search in a new neighbourhood
			op.reset();
		}

	}
	string odtfiledat; 

	if (saveodtfile) 
	{
		odtf.close();

		//write dat file
		odtfiledat = odtfile.substr(0,odtfile.find_last_of('.'))+".dat";

		odtf.open ( odtfiledat, std::ofstream::out);
		odtf << optcost << endl; // cost
		odtf << (gettime() - starttime) << endl; // time
		odtf << equalnets << endl; // networks
		odtf << improvements << endl; // improvements
		odtf << steps << endl; // steps
		odtf.close();

	}
	

	if (printstats)
	{
		cout << "Steps (networks evaluated): " << steps << endl;
		cout << "Improvements: " << improvements << endl;
		cout << "Equal-cost networks: " << equalnets << endl;
	}

	if (verbose>=1 && saveodtfile)
	{
		cout << "Optimal networks saved: " << odtfile << endl;	
		cout << "Stats data save to: " << odtfiledat << endl;
	}


	return optcost;

}