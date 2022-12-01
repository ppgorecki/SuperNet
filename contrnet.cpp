
#include "rtree.h"
#include "contrnet.h"
#include "network.h"

#define CON_VISITED 4
#define CON_DELPARENT 1
#define CON_DELRETPARENT 2

void ContractedNetwork::contractabove(SPID v)
{
#ifdef _CONTRDEBUG_
	cout << "contractabove:" << v << endl;
#endif	
	if (v==root) return;
	SPID p = mapup[parent[v]];
		
	contractedge(v,p);  // search

	if (v>=rtstartid)
	{
		addleftretusage(retdeleted, v-rtstartid);
		addrightretusage(retdeleted, v-rtstartid);
		p = mapup[retparent[v]];
		contractedge(v,p);  
	}
}

// void ContractedNetwork::propagate_maps(SPID v)
// {
// 	SPID par = mapup[v]; // tr/rt
//   	SPID bot = mapdn[v]; 

// #ifdef _CONTRDEBUG_
//   	cout << " propagate_maps " << v << " " << par << " " << bot << endl;
// #endif  	

//   	if (bot<rtstartid)
// 	{
// 		//bot is a tree node; one parent
// 		if (mapup[parent[bot]]==v) 
// 			mapup[parent[bot]]=par;		
// 	}
// 	else
// 	{
// 		//bot is reticulation
// 		if (mapup[parent[bot]]==v) 
// 			mapup[parent[bot]]=par;		
// 		else if (mapup[retparent[bot]]==v) 
// 			mapup[retparent[bot]]=par;		
// 	}
// }


void ContractedNetwork::propagate_maps(SPID v)
{
  	SPID par = mapup[v]; // tr/rt
  	SPID bot = mapdn[v]; 

#ifdef _CONTRDEBUG_
  	cout << " propagate_maps " << v << " " << par << " " << bot << endl;
#endif  	

	// propagede pred
  	if (par<nn) 
  	{
	  	if (par<rtstartid)
		{
			if (mapdn[leftchild[par]]==v) 
				mapdn[leftchild[par]]=bot;
			else if (mapdn[rightchild[par]]==v) 
					mapdn[rightchild[par]]=bot;
		}
		else
		{				
			if (mapdn[retchild[par]]==v) 
				mapdn[retchild[par]]=bot;				
		}
	}


	// propagede succ
  	if (bot<rtstartid)
	{
		//bot is a tree node; one parent
		if (mapup[parent[bot]]==v) 
			mapup[parent[bot]]=par;		
	}
	else
	{
		//bot is reticulation
		if (mapup[parent[bot]]==v) 
			mapup[parent[bot]]=par;		
		else if (mapup[retparent[bot]]==v) 
			mapup[retparent[bot]]=par;		
	}


}


void ContractedNetwork::contractedge(SPID v, SPID p)
{

#ifdef _CONTRDEBUG_
	cout << "contractedge:" << v << " " << p << endl;
#endif

	if (p>=nn) return; 
	
	if (p>=rtstartid) 
	{ 
		contractabove(p);
		mapdn[p]=mapdn[retchild[p]];
		return;
	}

	SPID lc = mapdn[leftchild[p]];
	SPID rc = mapdn[rightchild[p]];

#ifdef _CONTRDEBUG_
	cout << "           " << p << ":" << lc << " " << rc << endl;
#endif

	if (lc==v) 
	{ 

#ifdef _CONTRDEBUG_		
		cout << "ce-lc" << endl;
#endif		
		if (mapdn[p]!=p) contractabove(p);
	  	else { 
	  		mapdn[p]=rc;	   
	  		if (parent[p]<nn)
	  			mapup[p]=mapup[parent[p]];	   
	  		else mapup[p]=MAXSP;
	  		propagate_maps(p);	  		

	  	}
	   	return;
	}

	if (rc==v) 
	{

#ifdef _CONTRDEBUG_				
	  cout << "ce-rc" << endl;
#endif

   	  if (mapdn[p]!=p) contractabove(p);
   	  else { 
   	  	mapdn[p]=lc;	   
   	  	if (parent[p]<nn)
   	  		mapup[p]=mapup[parent[p]];	   
   	  	else mapup[p]=MAXSP;

   	  	propagate_maps(p);   	  	
   	  }
   	
	}
}

void ContractedNetwork::contract(RETUSAGE &retcontract)
{	
	SPID rtid = rtstartid;
	for (int i=0; i<rt; i++, rtid++)
	{	
		if (leftret(retcontract, i) && !leftret(retdeleted, i))  // (rtstartid+i, parent) -> remove
		{
			addleftretusage(retdeleted, i);
			SPID p = mapup[parent[rtid]];		
#ifdef _CONTRDEBUG_						
			cout << "\nRTL:" << " " << i << " r=" << rtid << " p=" << p << endl;
#endif			
			contractedge(rtid, p);										
			mapdn[rtid]= mapdn[retchild[rtid]];
			mapup[rtid]= mapup[retparent[rtid]];

			propagate_maps(rtid);
			
#ifdef _CONTRDEBUG_
			checkmaps();
			cout << newickrepr() << endl;
#endif			
		}

		if (rightret(retcontract, i) && !rightret(retdeleted, i))  
		{
			addrightretusage(retdeleted, i);			
			SPID p = mapup[retparent[rtid]];	
#ifdef _CONTRDEBUG_										
			cout << "\nRTR:" << i << " r=" << rtid << " p=" << p << endl;
#endif			
			contractedge(rtid, p);													
			mapdn[rtid]= mapdn[retchild[rtid]];
			mapup[rtid]= mapup[parent[rtid]];

			propagate_maps(rtid);
			
#ifdef _CONTRDEBUG_
			checkmaps();
			cout << newickrepr() << endl;
#endif			
		}
	}
	initdid();

	// correct root 
	while (mapdn[root]!=root) root=mapdn[root];

}


void ContractedNetwork::_checkmaps(SPID v, SPID p)
{
	SPID par = mapup[v]; 
  	SPID bot = mapdn[v]; 

  	if (mapup[v]!=v || mapdn[v]!=v)
  	{
  		cerr << "Non-contracted node - map error" << endl;
  		exit(-1);
  	}

  	if (v>=rtstartid)
  	{	
  		SPID cretchild = mapdn[retchild[v]];
  		SPID cparent = mapdn[parent[v]];
  		SPID cretparent = mapdn[retparent[v]];

  		if (v==cparent) _checkmaps(cretchild,v);
  		else
  		{
			if (v!=cretparent) 			
			{
				cerr << "Incorrect par-map in " << v << endl;
  				exit(-1);			
			}
  		}

  		return;
  	}

  	SPID bp = parent[v];
  	if (bp==MAXSP)
  	{
  		if (p!=MAXSP)
  		{ 
  			cerr << "Incorrect root-map (1) in " << v << endl;
			exit(-1);			
		}		
		// v is root		
  	}

  	if (p==MAXSP && bp!=MAXSP)
  	{
  		// v is root of contr. network, but not the root of the network
  		if (mapup[bp]!=MAXSP)
  		{
  			cerr << "Incorrect root-map (2) in " << v << endl;
			exit(-1);			
  		}  		
  	}

  	
	if (bp!=MAXSP && v!=mapup[bp])
	{
		cerr << "Incorrect leaf-parent map in " << v << endl;
		exit(-1);			
	}  		

	if (v<lf) return; // leaf case - done
  	
  	
  	// internal tree node
  	// check children  	
  	_checkmaps(mapdn[leftchild[v]],v);
  	_checkmaps(mapdn[rightchild[v]],v);

}


void ContractedNetwork::gendotcontracted(ostream &s)
{
 	int dagnum=2;

 	SPID cnodes[nn];
 	SPID cnodescnt=0;

    for (SPID i = 0; i < nn; i++ ) 
    {
	    	if (mapdn[i]==i && mapup[i]==i)
	    	{
	    		cnodes[cnodescnt++]=i;
				s << "v" << (int)i << "x" << dagnum << " [";
	    	
				if (i>=rtstartid)
	  			  			
	  				 s << "shape=box,color=red,";
	  		  			
		  		// if (mapdn[i]!=i || mapup[i]!=i)
		  		// 	s << "style=filled,fillcolor=\"#f0f0f0\",";

		  		s << "label=\"" << (int)i;  		
		  		if (i<lf)    
		  		  	s << " " << species(lab[i]) << " ";
		  		s << "\"]" << endl;  
  			}
  	}

  	for (SPID i = 0; i < nn; i++ ) 
    {
	    if (mapdn[i]==i && mapup[i]==i)
	    {
	    	

  			SPID _iparent = MAXSP;
  			while (getparentiter(i,_iparent))     
  			{
  				SPID iparent = _iparent;
  				while (iparent!=MAXSP && mapup[iparent]!=iparent) iparent = mapup[iparent];
  				if (iparent==MAXSP) continue;

		  		s << "v" << (int)iparent << "x" << dagnum << " -> v" << (int)i << "x" << dagnum;
		        s << " [ penwidth=1";
		          // if (leftchild[iparent]==i)  // leftchild - thicker edge
		          //   s << ",arrowhead=vee"; 
		          // else s << "";

				if (i>=rtstartid)
				    						    
				   s << ",label=\"" << spid2retlabel[i] << "\"";				    		   	
				    	     
		  		s << "]" << endl;
		  	}
		  }

	}
}

void ContractedNetwork::gendot(ostream &s)
{    
    int dagnum=1;
    for (SPID i = 0; i < nn; i++ ) 
    {     	
		s << "v" << (int)i << "x" << dagnum << " [";
		if (i>=rtstartid)
  		{
  			if (leftret(retdeleted, i-rtstartid) && rightret(retdeleted, i-rtstartid))  			
  				s << "shape=box,color=\"#f0f0f0\",";
  			else s << "shape=box,color=red,";
  		}

  		if (mapdn[i]!=i || mapup[i]!=i)
  			s << "style=filled,fillcolor=\"#f0f0f0\",";

  		s << "label=\"" << (int)i;
  		if (mapup[i]!=i) s << " u" << mapup[i] << " "; 
  		if (mapdn[i]!=i) s << " d" << mapdn[i] << " "; 
  
  		if (i<lf)    
  		  	s << " " << species(lab[i]) << " ";
  		s << "\"]" << endl;  

  		SPID iparent = MAXSP;
  		while (getparentiter(i,iparent))     
  		{
  		  s << "v" << (int)iparent << "x" << dagnum << " -> v" << (int)i << "x" << dagnum;
          s << " [ penwidth=1";
          if (leftchild[iparent]==i)  // leftchild - thicker edge
            s << ",arrowhead=vee"; 
          else s << "";

		  if (i>=rtstartid)
		    	{	
		    		if (iparent==parent[i] && leftret(retdeleted, i-rtstartid)) s << ",style=dashed";  // from a leaf #A  		    		
		    		s << ",label=\"" << spid2retlabel[i];
		    		if (iparent==retparent[i]) s<<" rp"; 
		    			s << "\"";  
		    		if (iparent==retparent[i] && rightret(retdeleted, i-rtstartid)) 
		    			s << ",style=dotted";  // from a leaf #A  		    				    		
		    	}     
		    else
		    {				
		    	if (mapdn[i]!=i || mapup[i]!=i)
		    	{
		    		//s << ",color=\"#808080\"";
		    	}
		    }   
  		    s << "]" << endl;
  		}

	}

	s << " info [ shape=plaintext, label=\"";
	print(s);
	s << "\"] " << endl;
	s << "info2 [ shape=plaintext, label=\"" << retdeleted << " DTmax=" << displaytreemaxid();
	s << "\"] " << endl;
	s << "info2  ->  info" << endl;
	s << "info" << " -> " << "v" << (int)root << "x" << dagnum  << endl;	
}


SPID ContractedNetwork::rtcount()
{
	SPID cnt=rt;
	for (SPID i=0; i<rt; i++)		
		 if (leftret(retdeleted,i) || rightret(retdeleted,i)) cnt--;
	return cnt;
}


string ContractedNetwork::_newickrepr(SPID v, SPID p)
{
	while (mapdn[v]!=v) v = mapdn[v];
	if (v<lf) return specnames[lab[v]];
	if (v<rtstartid)
		return "("+_newickrepr(leftchild[v],v)+","+_newickrepr(rightchild[v],v) +")";

	// retid
	SPID rp = retparent[v];
	while (mapup[rp]!=rp) rp = mapup[rp];
	if (rp==p)
		return "("+_newickrepr(retchild[v],v)+")"+spid2retlabel[v];
	return spid2retlabel[v];
	
		
}

string ContractedNetwork::newickrepr()
{	
	return _newickrepr(root,MAXSP);
}

void ContractedNetwork::_init()
{	
	emptyretusage(retdeleted);
	mapdn = new SPID[nn];  
	mapup = new SPID[nn];  	
	for (int i = 0; i<nn; i++)
	{
		mapdn[i]=i;
		mapup[i]=i;
	}
	initdid();
	localbitmask=NULL;
}


ContractedNetwork::ContractedNetwork(ContractedNetwork &net, int shallowcopy) : Network(net, shallowcopy) 
{ 
	mapdn = new SPID[nn];  
	mapup = new SPID[nn];  	
	memcpy ( mapdn, net.mapdn, sizeof(SPID)*nn );
	memcpy ( mapup, net.mapup, sizeof(SPID)*nn );
	retdeleted = net.retdeleted;
	localbitmask=NULL;
}	

ContractedNetwork::~ContractedNetwork()
{
	delete[] mapup;
	delete[] mapdn;
	if (localbitmask)
		delete[] localbitmask;
}


bool ContractedNetwork::getnodeiter(SPID &i)
{
  if (i==MAXSP) 
  { 
    i=0; 
    return true;  // leaf
  }
  i++;
  while (i<nn && mapdn[i]!=i) i++;  
  return i<nn;  
}

SPID ContractedNetwork::getconflictedreticulation(RETUSAGE &retusage)
{
	SPID r = rtstartid;
	for (SPID i=0; i<rt; i++,r++)
	{
		if (leftret(retusage,i) && rightret(retusage,i))
			if (mapdn[r]==r) return i;
	}
	cerr << "No conflicted reticulations located." << endl;
	return MAXSP;
}

extern DISPLAYTREEID bitmask[]; 

RootedTree* ContractedNetwork::gendisplaytree(DISPLAYTREEID id, RootedTree *t)
{
	// prepare local bitmask
	if (!id)
	{			
		if (!localbitmask) 
			localbitmask = new DISPLAYTREEID[rt];
		int bit = 0;
		for (SPID i=0; i<rt; i++)
		{
			if (!(leftret(retdeleted,i) || rightret(retdeleted,i)))
			{
				localbitmask[i]=bitmask[bit];											
				bit++;
			}
		}
	}

	return Network::gendisplaytree(id, t);
}

bool ContractedNetwork::_skiprtedge(SPID i, SPID iparent, DISPLAYTREEID id)
{
	return ((iparent != getparent(i)) == bool(id & localbitmask[i-rtstartid]));
}
