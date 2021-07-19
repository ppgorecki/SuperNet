

#include "tools.h"
#include "dag.h"

void Dag::init(int _lf, int _rt)
{
    if (2*_lf+2*_rt-1 >= MAXSP) 
    {
      cout << "Too many nodes in a dag (lf=" << _lf << " rt=" << _rt <<"). Compile with SPLARGE or SPMED macro." << endl;
      exit(-1);
    }    

    lf=_lf;
    rt=_rt;
    nn=2*lf+2*rt-1; 

    retchild = leftchild = new SPID[nn];
    rightchild = new SPID[nn - lf - rt]; // proper indexing   
    rightchild -= lf;
    rtstartid = nn - rt;

    parent = new SPID[nn];
    lab = new SPID[nn];

    if (rt)
    {
      retparent = new SPID[rt];    
      spid2retlabel  = new string[rt];
      retparent -= rtstartid; // shift to obtain easy adressing   
      spid2retlabel -= rtstartid;
    }
    else     
    {
      retparent = NULL; 
      spid2retlabel = NULL;
    }
}

void Dag::parse(char *s)
{
    int _rt = strcount(s,'#')/2;
    int _lf = strcount(s,',')+1-_rt;
    //cerr << _lf << endl;
    Dag::init(_lf, _rt);
    int p = 0;    

    SPID freeleaf = 0;
    SPID freeint = lf;
    SPID freeret = rtstartid; // last rt nodes are reticulations
    SPID *_p = parent; 
    
    map<string, SPID> retlabel2spid;
    root = _parse(s, p, 0, freeleaf, freeint, freeret, _p, retlabel2spid);   
    parent[root] = MAXSP;     
    checkparsing(s+p); 
}

SPID Dag::_parse(char *s, int &p, int num, 
    SPID &freeleaf, SPID &freeint, SPID &freeret,         
    SPID* &parentset, map<string, SPID> &retlabel2spid)
{
  char *cur = getTok(s, p, num);
  
  if (cur[0] == '(')
  {
    SPID *parenta = parent;
    SPID *parentb = parent;

    int a = _parse(s, p, num, freeleaf, freeint, freeret, parenta,retlabel2spid);
    char *token = getTok(s, p, num);
    if (token[0]==')')
    {
      char *cur = getTok(s, p, num);        
      string t = getstringn(cur, s + p - cur);
      SPID retnode; 
      if (retlabel2spid.count(t)) retnode = retlabel2spid[t];
      else retlabel2spid[t] = retnode = freeret++;             
      spid2retlabel[retnode] = t;
      retchild[retnode] = a;
      parenta[a] = retnode;                
      return retnode;
    }
    else
    {        
      int b = _parse(s, p, num, freeleaf, freeint, freeret, parentb,retlabel2spid);
      getTok(s, p, num);
      if (freeint >= nn)
      {
        cerr << "Is it a binary network? Too many nodes in rooted network" << endl;
        exit(-1);
      }
      rightchild[freeint] = b;
      leftchild[freeint] = a;              
      parenta[a] = freeint;
      parentb[b] = freeint;        
      return freeint++;
    }
  }

  // leaf processing 
  if (freeleaf >= lf)
  {
    cerr << "Too many leaves in a binary network. " << endl; 
    exit(-1);
  }
  
  if (cur[0]=='#') 
  {
    // network leaf      
    parentset = retparent; // use reticulation parent 
    string t = getstringn(cur, s + p - cur);
    if (retlabel2spid.count(t)) return retlabel2spid[t];    
    retlabel2spid[t] = freeret;
    return freeret++;
  }
  
  // normal leaf
  lab[freeleaf] = getspecies(cur, s + p - cur);        
  return freeleaf++;
 
}


Dag::Dag(char *s, double dagweight): weight(dagweight)
{
	parse(s);
	setexactspecies();
}


Dag::Dag(int _lf, SPID *labels, double dagweight) : weight(dagweight)
{
  init(_lf,0);

  // set labels
  for (SPID i=0; i<lf; i++) lab[i] = labels[i]; 
  
  if (lf>1)
    {

        SPID left = 0; //left child
        SPID in = lf;  //parent
        for (SPID right=1; right<lf; right++) // iterate over right children
        {
            leftchild[in] = left;
            rightchild[in] = right;
            parent[left] = parent[right] = in;             
            left = in++;            
        }         
    }
    root = nn-1;
    parent[root] = MAXSP;
}

// Find leaf label
SPID Dag::findlab(SPID slab, int stoponerr)
{
    if (exactspecies && slab < lf) return slab;
      
    for (SPID i=0; i<lf; i++)
    	if (slab == lab[i]) 
    		return i;

    if (stoponerr)
    {
      cerr << "Species lab " << specnames[slab] << " ($" << (int)slab << ") not found in: " << *this << endl;

      exit(-1);
    }
    return MAXSP;
}

// Returns true if the leaf labelling is bijective
bool Dag::bijectiveleaflabelling()
{            
  SPID occ[specnames.size()];
  for (size_t i=0; i<specnames.size(); i++) occ[i]=0;
  for (SPID i=0; i<lf; i++) 
    if (occ[lab[i]]++ == 1) return false;
  return true;
}

ostream& Dag::printdot(ostream&s)
{
    for (SPID i = 0; i < size(); i++ ) 
    {     
		s << "v" << (int)i;
		if (i<lf)    
		  s << " [label=\"" << species(lab[i]) << " " << (int)i << "\"]" << endl;  
		else if (i>=rtstartid)
		{
			s << " [shape=box,color=red,label=\"" << (int)i << "\"]" << endl;  
		}
		else
			s << " [label=\"" << (int)i << "\"]" << endl;  	
		s << endl;

		SPID iparent = MAXSP;
		while (getparent(i,iparent))          
		{
		    s << "v" << (int)iparent << "->  v" << (int)i;      
		    if (i>=rtstartid)
		    	{	
		    		s << " [label=\"" << spid2retlabel[i] << "\"," << endl;  
		    		if (iparent==retparent[i])
		    			s << " color=red]";  // from a leaf #A
		    		else
		    			s << " color=green]"; // from internal #A  		    	
		    	}
		    s << endl;
		}

	}
	return s;
}



ostream& Dag::printsubtree(ostream&s, SPID i, SPID iparent) 
{ 
	// leaf
	if (i < lf) return s << species(lab[i]); 

	// internals via getchild	

	if (i>=rtstartid)	
		// reticulation 
		if (iparent == retparent[i]) 
			// convert to leaf and return if coming from retparent[i]
			return   s << spid2retlabel[i];

	s << "(";
	bool first = true;
	SPID ichild = MAXSP;
	while (getchild(i,ichild))          
	{        
	    if (!first) s << ",";
	    else first = false;
	    printsubtree(s, ichild, i);
	}
	s << ")";

	if (i>=rtstartid) 
		// reticulation; insert reticulation label
		s << spid2retlabel[i];
	return s; 
}


// Returns the parents; to get all parents use:
// SPID p=MAXSP;
// while (getparent(i,p)) { .. p is a parent ... }
bool Dag::getparent(SPID i, SPID &rparent)
{
	if (rparent==MAXSP) 
	{
		// first run  		
		rparent = parent[i];
		return rparent!=MAXSP;  	
	}

	if (rparent==parent[i])
	{
		// second run; try reticulation
		if (i>=rtstartid) 
		{
			rparent=retparent[i];
			return true;
		}
	}
	return false;
}

// Returns the parents; to get all parents use:
// SPID ichild=MAXSP;
// while (getchild(i,ichild)) { .. ichild is a child of i ... }
bool Dag::getchild(SPID i, SPID &ichild)
{
	if (i<lf) return false;

	if (ichild==MAXSP) 
	{
		// first run  		
		ichild = leftchild[i]; // also retchild[i] 
		return true;
	}

	if (ichild==leftchild[i] && i<rtstartid)
	{
		// second run; 
		ichild=rightchild[i];
		return true;		
	}

	return false;
}


void Dag::setexactspecies()
{
    exactspecies = 1;
    for (SPID i=0; i<lf; i++) 
    	if (lab[i] != i) 
    		{ 
    			exactspecies = 0; 
    			break; 
    		}
}

ostream& Dag::printdebstats(ostream&s) 
{
    s << " Nodes(nn)=" << (int)size() 
      << " Leaves(lf)=" << (int)lf ;
    if (rt)      
        s << " Reticulations(rt)=" << (int)rt;
    
    s << endl;

    s << " Leaves:0.."<< (int)(lf-1)
      << " TreeNodes:" << (int)lf << ".." << (int)rtstartid-1;
    if (rt) s << " RetNodes:" << (int)rtstartid << ".." << nn-1;

    s << endl;

    if (rt) s << " rtstartid=" << (int)rtstartid;
    s << " root=" << (int)root;           
    s << " exactspecies=" << exactspecies;
    return s;
}

ostream& Dag::printdebarrays(ostream&s) 
{
    s << " parent= "; 
    for (SPID i=0; i<nn; i++) s << " " << (int)i << ":" << parent[i]; 
    s << endl;

    s << " leftchild= "; 
    for (SPID i=lf; i<rtstartid; i++) s << " " << (int)i << ":" << leftchild[i]; 
    s << endl;

    s << " rightchild= "; 
    for (SPID i=lf; i<rtstartid; i++) s << " " << (int)i << ":" << leftchild[i]; 
    s << endl;

    s << " lab= "; 
    for (SPID i=0; i<lf; i++) s << " " << int(i) << ":" << lab[i]; 
    s << endl;   

    if (rt)
    {
      s << " retchild= "; 
      for (SPID i=rtstartid; i<nn; i++) s << " " << (int)i << ":" << retchild[i]; 
      s << endl;

      s << " retparent= "; 
      for (SPID i=rtstartid; i<nn; i++) s << " " << (int)i << ":" << retparent[i]; 
      s << endl;   

      s << " spid2retlabel= "; 
      for (SPID i=rtstartid; i<nn; i++) s << " " << int(i) << ":" << spid2retlabel[i]; 
      s << endl;   


    }
      
    return s;
}

ostream& Dag::printdeb(ostream&s, int gse, string tn) 
{
	if (gse&2)
	{	  
    printdebstats(s);

    s << endl;
	  for (SPID i = 0; i < size(); i++ ) 
	  {	  		
	  	   	s << setiosflags(ios::left) << setw(2) << (int)i << " ";	  	   	
	  	   		  	   	
    		if (parent[i] == MAXSP) 
    			s << "Root";
    		else { 
    			if (i<lf) 
    				s << "Leaf";
    			else if (i>=rtstartid)
    				s << "Reti";
    			else s << "Tree";
    			// s  << " par=" << setw(3) << (int)parent[i];
    		} 	

    		SPID ic = MAXSP;    		
    		while (getchild(i,ic)) s << setw(0) << " c=" << (int)ic;    			
   		
    		SPID ip = MAXSP;
    		while (getparent(i,ip)) s << setw(0) << " p=" << (int)ip;   	    	

    		if (i>=rtstartid) 
    		{
    			//s << setw(3) << (int)retparent[i];
    			s << " retlabel=" << spid2retlabel[i];
    		}
    		
    		if (i < lf)
      			s << " $" << (int)lab[i] << " " << species(lab[i]);
    			
    		// else if (i<rtstartid)
    		// {
    		// 	s << " ch=" << setw(3) << (int)leftchild[i] << setw(3) << (int)rightchild[i];    		
    		// } else
    		// 	s << " ch=" << setw(3) << (int)retchild[i] << setw(3) << "";
    		
    		s << " subtree=";
    		printsubtree(s,i); 
    		s << endl;
	  }	    

    printdebarrays(s);
	}
  
    // if (gse&1)
    // {
    //   s << "&s (";
    //   char cnt='A';
    //   NODESLOOP if (par[i]==MAXSP)
    //   {    
    //        if (cnt!='A') s << ",";          
    //        printsubtree(s,i) << " edgelab=\"" << cnt << "\"" ;
    //        cnt++;
    //   }
    //   s << ")" << " treename=\"" << _gsetrees++ << " " << tn << "\"" ;
    // }
    return s;
  }