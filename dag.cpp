
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

    retchild = leftchild = new SPID[nn-lf];
    retchild -= lf;
    leftchild -= lf;
    rightchild = new SPID[nn - lf - rt]; 
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
      if (cur[0]!='#')
      {
        cerr << "Parse error: # expected in " << cur << endl;
        exit(-1);
      }
      
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
      token = getTok(s, p, num);

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
  
  if (cur[0]=='#') 
  {
    // network leaf      
    parentset = retparent; // use reticulation parent 
    string t = getstringn(cur, s + p - cur);
    if (retlabel2spid.count(t)) return retlabel2spid[t];    
    retlabel2spid[t] = freeret;
    return freeret++;
  }


  // leaf processing 
  if (freeleaf >= lf)
  {
    cerr << "Too many leaves in a binary network. " << endl; 
    exit(-1);
  }
  
  // normal leaf
  lab[freeleaf] = getspecies(cur, s + p - cur);        
  return freeleaf++;
 
}


Dag::Dag(const char *s, double dagweight): weight(dagweight)
{
  char *t = strdup(s);
	parse(t);
  free(t);
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

  setexactspecies();


}

// Find leaf label
SPID Dag::findlab(SPID slab, int stoponerr)
{
    if (exactspecies && (slab < lf)) return slab;
      
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

int Dag::verifychildparent()
{

   if (parent[root]!=MAXSP) { 
      cerr << (int)root << ": wrong parent (" << (int)parent[root] << ")of the root ";
      return 1; // wrong parent of the root
   }


   if (lf==1 && nn==1) return 0; // single noded tree; OK

   for (SPID i = 0; i < nn; i++ ) 
   {
      SPID p = parent[i];

      if (p==MAXSP && i!=root) { 
        cerr << (int)i << ": non-root parent is MAXSP";
        return 2; // incorect parent 
      }
      if (p==MAXSP && i>=rtstartid) { 
        cerr << (int)i << ": reticulation parent is MAXSP";
        return 3; 
      }

      if (parent[i]==i) { 
        cerr << (int)i << " and its parent are looped";
        return 11; // self - loop
      }
      if (i<lf) continue;
      if (i<rtstartid)
      {
        if (leftchild[i]==i) { 
          cerr << (int)i << " and its leftchild are looped";
          return 12; // self - loop
        }
        if (rightchild[i]==i) { 
            cerr << (int)i << " and its rightchild are looped";
            return 13; // self - loop
        }
      } else {        
        if (retchild[i]==i) { 
          cerr << (int)i << " and its retchild are looped";
          return 14; // self - loop      
        }
      }       
   }

   for (SPID i = 0; i < nn; i++ ) 
   {
      SPID p = parent[i];
      if (p==MAXSP) continue;      
      if (p<lf) { 
        cerr << (int)i << ": parent " << (int)p << "is a leaf?";
        return 23; // Parent is a leaf?
      }

      if (p>=rtstartid) // parent is a reticulation
      {
        if (retchild[p]!=i) 
        {
          cerr << (int)i << ": parent->retchild disconnected " << (int)p;
          return 24; // wrong retchild
        }
      }
      else
      {
        if (i==leftchild[p] && i==rightchild[p]) { 
            cerr << (int)i << ": parent[i] has equal left and right children" << (int)p;
            return 25; // double edge? 
        }
        if (i==leftchild[p] || i==rightchild[p]) continue; // OK        
        cerr << (int)i << ": no connection from parent[i]" << (int)p;
        return 26; // no child in the parent
      }
   }

   // check ret. node (ugly repetition)
   for (SPID i = rtstartid; i < nn; i++) 
   {
      SPID rp = retparent[i];
      if (rp == MAXSP) { 
        cerr << (int)i << ": retparent[i] is MAXSP";
        return 4; // retparent is undefined
      }
      if (rp == parent[i]) { 
        cerr << (int)i << ": retparent[i] equals parent[i]";
        return 5; // double edge?
      }

      if (rp>=rtstartid) // parent is a reticulation
      {
        if (retchild[rp]!=i) 
        {
            cerr << (int)i << ": wrong retchild of retparent[i]";
            return 34; // wrong retchild
        }
      }
      else
      {
        if (i==leftchild[rp] && i==rightchild[rp]) { 
          cerr << (int)i << ": retparent[i] has equal left and right children" << (int)rp;
          return 35; // double edge? 
        }
        if (i==leftchild[rp] || i==rightchild[rp]) continue; // OK 

        cerr << (int)i << ": no connection from retparent[i]" << (int)rp;       
        return 36; // no child in the parent
      }

   }

   return 0; // OK

}

ostream& Dag::printdot(ostream&s, int dagnum)
{
    for (SPID i = 0; i < size(); i++ ) 
    {     
		
      s << "v" << (int)i << "x" << dagnum;
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
  		    s << "v" << (int)iparent << "x" << dagnum << " -> v" << (int)i << "x" << dagnum;

          s << " [ penwidth=";
          if (leftchild[iparent]==i)  // leftchild - thicker edge
            s << "3"; 
          else s << "1";

  		    if (i>=rtstartid)
  		    	{	
  		    		s << ",label=\"" << spid2retlabel[i] << "\"";  
  		    		if (iparent==retparent[i])
  		    			s << ",color=red";  // from a leaf #A  		    		
  		    	}        
  		    s << "]" << endl;
  		}

	}
	return s;
}



ostream& Dag::printsubtree(ostream&s, SPID i, SPID iparent, int level) 
{ 
	// leaf
  if (level>nn)
  {
    // stop infinite prints
    s << "...loop?" << endl;
    return s;
  }
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
  int childcnt=0;
	while (getchild(i,ichild))          
	{        
	    if (!first) s << ",";
	    else first = false;
	    printsubtree(s, ichild, i, level+1);
      childcnt++;
      if (childcnt>2)
      {
        s << ",...?"; 
        break;
      }
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
      if (retparent[i]!=rparent) 
      {
         rparent=retparent[i];      
			   return true;
      }
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

	if (i<rtstartid && ichild==leftchild[i] && leftchild[i]!=rightchild[i])
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
    for (SPID i=lf; i<rtstartid; i++) s << " " << (int)i << ":" << rightchild[i]; 
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
    		while (getchild(i,ic)) 
          s << setw(0) << " c=" << (int)ic;    			                  
   		
    		SPID ip = MAXSP;        
    		while (getparent(i,ip)) 
          s << setw(0) << " p=" << (int)ip;   	    	
          
       

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

Dag::Dag(Dag *d, SPID v, SPID p, SPID w, SPID q, string retid, double dagweight)
{
    init(d->lf,d->rt+1); 

    for (SPID i=0; i<d->lf; i++) lab[i] = d->lab[i]; 

#define CID(id)  ((id)>=d->rtstartid)?((id)+1):(id)

    // leaves and tree nodes - parent
    for (SPID i=0; i<d->rtstartid; i++)            
      parent[i] = CID(d->parent[i]);      
    
    // ret. nodes - parent
    for (SPID i=d->rtstartid; i<d->nn; i++)        
      parent[i+1]=CID(d->parent[i]);     

    // children of tree nodes
    for (SPID i=lf; i<d->rtstartid; i++)  
    {          
      leftchild[i] = CID(d->leftchild[i]);
      rightchild[i] = CID(d->rightchild[i]);
    }

    // children of retnodes
    for (SPID i=d->rtstartid; i<d->nn; i++)            
      retchild[i+1] = CID(d->retchild[i]);

    // children of retnodes
    for (SPID i=d->rtstartid; i<d->nn; i++)            
      retparent[i+1] = CID(d->retparent[i]);

    // ret labels        
    for (SPID i=d->rtstartid; i<d->nn; i++)            
      spid2retlabel[i+1] = d->spid2retlabel[i];

    root = d->root;

    parent[root] = MAXSP;

    if (!retid.length())    
        retid = "#"+std::to_string(rt);
    

    // insert ret. edge 
    // new retid = nn-1 
    // new treenode id = rtstartid-1

    SPID nr = nn - 1;
    SPID nt = rtstartid - 1;

    spid2retlabel[nr] = retid;


    p=CID(p);
    v=CID(v);
    q=CID(q);
    w=CID(w);

    parent[nt] = p;

    parent[nr] = nt;
    retparent[nr] = q;

    leftchild[nt] = nr; 
    rightchild[nt] = v;
    retchild[nr] = w;

    if (leftchild[p]==v) leftchild[p]=nt;
    else rightchild[p]=nt;

    if (parent[v]==p) parent[v]=nt;
    else retparent[v]=nt;

    if (leftchild[q]==w) leftchild[q]=nr;
    else rightchild[q]=nr;

    if (parent[w]==q) parent[w]=nr;
    else retparent[w]=nr;

    verifychildparent();

    setexactspecies();

}

// copy constructor
Dag::Dag(const Dag &d)
{
  init(d.lf,d.rt);
  memcpy ( lab, d.lab, sizeof(SPID)*lf );
  memcpy ( parent, d.parent, sizeof(SPID)*nn );
  memcpy ( leftchild+lf, d.leftchild+lf, sizeof(SPID)*(nn-lf) );
  memcpy ( rightchild+lf, d.rightchild+lf, sizeof(SPID)*(nn-lf-rt) );
  if (rt)
    memcpy ( retparent+rtstartid, d.retparent+rtstartid, sizeof(SPID)*(nn-rtstartid) );
  

  for (int i=rtstartid;i<nn;i++)
    spid2retlabel[i]=d.spid2retlabel[i];

  root=d.root;  
  weight=d.weight;     
  exactspecies=d.exactspecies; 

  // if (verifychildparent())
  //   { exit(-1); } 
  
}

Dag::~Dag() 
  {
    
    leftchild+=lf;
    delete[] leftchild;
    rightchild+=lf;
    delete[] rightchild;
    delete[] parent;
    delete[] lab;
    if (rt)
    {
      retparent+=rtstartid;
      spid2retlabel+=rtstartid;
      delete[] retparent;
      delete[] spid2retlabel;
    }


  }
