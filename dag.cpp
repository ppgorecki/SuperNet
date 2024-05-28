
#include <algorithm>
#include <iterator>
#include<bits/stdc++.h>
using namespace std;
#include "tools.h"
#include "dag.h"
#include "topsort.h"

#include <bits/stdc++.h>
using namespace std;

extern string flag_retidprefix;

void Dag::init(int _lf, int _rt)
{
    if (2*_lf+2*_rt-1 >= MAXNODEID) 
    {
      cout << "Too many nodes in a dag (lf=" << _lf << " rt=" << _rt <<"). Compile with SPMED macro." << endl;
      exit(-1);
    }    

    lf=_lf;
    rt=_rt;
    nn=2*lf+2*rt-1; 

#define GUARD 1  

    leftchild = new NODEID[nn-lf+GUARD];  
    leftchild -= lf;
    retchild = leftchild;

    rightchild = new NODEID[nn - lf - rt + GUARD]; 
    rightchild -= lf;
    rtstartid = nn - rt;

    parent = new NODEID[nn+GUARD];
    lab = new NODEID[nn+GUARD];

    if (rt)
    {
      retparent = new NODEID[rt+GUARD];    
      spid2retlabel  = new string[rt+GUARD];
      retparent -= rtstartid; // shift to obtain easy adressing   
      spid2retlabel -= rtstartid;
    }
    else     
    {
      retparent = NULL; 
      spid2retlabel = NULL;
    }
    shallow = false;
}

void Dag::parse(char *s)
{
    int _rt = strcount(s,'#')/2;
    int _lf = strcount(s,',')+1-_rt;
    //cerr << _lf << endl;
    Dag::init(_lf, _rt);
    int p = 0;    

    NODEID freeleaf = 0;
    NODEID freeint = lf;
    NODEID freeret = rtstartid; // last rt nodes are reticulations
    NODEID *_p = parent; 
    
    map<string, NODEID> retlabel2spid;
    root = _parse(s, p, 0, freeleaf, freeint, freeret, _p, retlabel2spid);   
    parent[root] = MAXNODEID;     
    checkparsing(s+p); 
}


NODEID Dag::_parse(char *s, int &p, int num, 
    NODEID &freeleaf, NODEID &freeint, NODEID &freeret,         
    NODEID* &parentset, map<string, NODEID> &retlabel2spid)
{
  char *cur = getTok(s, p, num);

  if (cur[0] == '(')
  {
    NODEID *parenta = parent;
    NODEID *parentb = parent;

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
      NODEID retnode; 
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
    if (retlabel2spid.count(t)) 
      return retlabel2spid[t];    
    retlabel2spid[t] = freeret;
    return freeret++;
  }


  // leaf processing 
  if (freeleaf >= lf)
  {
    
    cerr << "Too many leaves in a binary network. " << endl; 
    exit(-1);
  }
  
  // leaf  
  lab[freeleaf] = getspecies(cur, s + p - cur);    
  char *cur2 = seeTok(s, p, num);
  if (cur2[0]!=',' && cur2[0]!=')')
  {     
    cerr << ", or ) expected after leaf definition " << cur << endl;
    exit(-1);
  }
  
  return freeleaf++;
 
}


Dag::Dag(const char *s, double dagweight): weight(dagweight)
{
  char *t = strdup(s);
	parse(t);
  free(t);
	setexactspecies();
}


inline bool hrcompare(const std::pair<NODEID, NODEID> &p1, const std::pair<NODEID, NODEID> &p2)
{
    return p1.first < p2.first;
}


void Dag::_dagrtreplace(NODEID s, NODEID d)
{
  cout << "REPL: " << s << "->" << d << endl; 

  NODEID p = parent[s];
  NODEID q = retparent[s];
  NODEID c = retchild[s];

  if (p>=rtstartid) retchild[p]=d;
  else if (leftchild[p]==s) leftchild[p]=d;
       else rightchild[p]=d;

  if (q>=rtstartid) retchild[q]=d;
  else if (leftchild[q]==s) leftchild[q]=d;
       else rightchild[q]=d;

  if (c>=rtstartid) 
    { 
      if (parent[c]==s) parent[c]=d;
      else retparent[c]=d; 
    }
  else parent[c]=d;

  retchild[d] = retchild[s];
  parent[d] = parent[s];
  retparent[d] = retparent[s];

}

void Dag::sortrtnodes()
{
  if (rt<=1) return;

  int height[nn];
  memset(height,0,sizeof(int)*nn);
  for (int i=0;i<lf;i++) height[i]=1;
  for (NODEID i=rtstartid;i<nn;i++) _height(i, height);
  vector<pair<NODEID, NODEID>> ren;
  int j=0;
  for (NODEID i=rtstartid; i<nn; i++)
    ren.push_back(make_pair(height[i],i));

  std::sort( std::begin( ren ), std::end( ren ), hrcompare );
  cout << rtstartid << " " << nn << endl;
  for (int i=0; i<ren.size(); i++)
  {        
        cout << i+rtstartid << " " << ren[i].first << " "
             << ren[i].second << endl;        
  }

  // set first to the destination rtid
  for (int i=0; i<ren.size(); i++)  
    ren[ren[i].second-rtstartid].first=i+rtstartid;
  
  cout << "====" << endl;
  for (int i=0; i<ren.size(); i++)
  {
      NODEID newid = ren[i].first; 
      cout << i+rtstartid << " " <<  ren[i].first  << " "  << "newid=" << newid << " " 
             << ren[i].second << endl;        
  }

  // gen rtid permutation cycles
  NODEID rtcycles[rt*2];
  rtcycles[0]=0; // guard
  NODEID rtc = 1;
  for (int i=0; i<ren.size(); i++)
  {
       if (ren[i].second==0) continue;
       NODEID start = i+rtstartid;
       NODEID next = ren[i].first;
       if (start==next) continue; // singleton cycle
       rtcycles[rtc++] = start;       
       NODEID prev = start;
       while (start!=next)
       {
          rtcycles[rtc++]=next;
          ren[next-rtstartid].second=0;
          prev = next;
          next = ren[next-rtstartid].first;            
       }
       rtcycles[rtc++]=0; // guard
  }


  // NODEID *parent;  // Parent array
  // NODEID *leftchild, *rightchild; // Left and right child arrays
  // NODEID *retchild, *retparent;   

  cout << "Cycles: ";
  for (int i=0; i<rtc; i++)
    cout << rtcycles[i] << " "; 
  cout << endl;



  


  NODEID i = rtc-1;
  while (i>0)
  {
    if (rtcycles[i]==0)
    { 
      i--;
      NODEID last = rtcycles[i];      
      _dagrtreplace(last,nn);
      i--;
      while (rtcycles[i])
      {                        
        _dagrtreplace(rtcycles[i],rtcycles[i+1]);        
        i--;
      }
      _dagrtreplace(nn,rtcycles[i+1]);        
    }
  }
    
  //     NODEID p = parent[prev];
  //     NODEID q = retparent[prev];
  //         NODEID c = retchild[prev];

  //         if (rp>=rtstartid) retchild[p]=next;
  //         else if (leftchild[p]==prev) leftchild[p]=next;
  //              else rightchild[p]=next;

  //         if (leftchild)

  // }
  // // cout << endl;



  string _spid2retlabel[rt];
  for (NODEID oldid=rtstartid; oldid<nn; oldid++)
  {
    NODEID newid = ren[oldid-rtstartid].first; 
    _spid2retlabel[newid-rtstartid] = spid2retlabel[oldid];
    // replace oldid by newid
  }
  for (NODEID rtid=rtstartid; rtid<nn; rtid++)
  {
    spid2retlabel[rtid] = _spid2retlabel[rtid-rtstartid];    
  }


  verifychildparent();



  // NODEID *parent;  // Parent array
  // NODEID *leftchild, *rightchild; // Left and right child arrays
  // NODEID *retchild, *retparent;   
  // // In implementation retchild == leftchild 
  // // Reticulation child and parent (only for ret. nodes)

  // string* spid2retlabel;        // Dict. of reticulation labels (spid -> string)


  // convert?



}

inline int Dag::_height(NODEID v, int heightarr[])
{
  if (v<lf) return 1;
  if (!heightarr[v]) 
  {
    if (v>=rtstartid)   
       heightarr[v]=_height(retchild[v],heightarr)+1;
    else  
    {
      int a = _height(leftchild[v],heightarr);
      int b = _height(rightchild[v],heightarr);
      heightarr[v]=1+max(a,b);
    }
  }
  return heightarr[v];
    
}

Dag::Dag(int _lf, NODEID *labels, double dagweight) : weight(dagweight)
{
  init(_lf,0);

  // set labels
  for (NODEID i=0; i<lf; i++) lab[i] = labels[i]; 
  
  if (lf>1)
    {

        NODEID left = 0; //left child
        NODEID in = lf;  //parent
        for (NODEID right=1; right<lf; right++) // iterate over right children
        {
            leftchild[in] = left;
            rightchild[in] = right;
            parent[left] = parent[right] = in;             
            left = in++;            
        }         
    }
    root = nn-1;
    parent[root] = MAXNODEID;

  setexactspecies();


}

// Find a leaf by a label
NODEID Dag::findlab(NODEID slab, int stoponerr)
{
    if (exactspecies && (slab < lf)) return slab;
      
    for (NODEID i=0; i<lf; i++)
    	if (slab == lab[i]) 
    		return i;

    if (stoponerr)
    {
      
      cerr << "Species lab id:" << slab;
      if (slab>=0 && slab<specnames.size()) 
        cerr << " " << specnames[slab] << " ";
      cerr << " not found in "
           << *this << endl;  

      exit(-1);
    }
    return MAXNODEID;
}

// Returns true if the leaf labelling is bijective
// TODO: optimize to attribute
bool Dag::bijectiveleaflabelling()
{            
  NODEID occ[specnames.size()];
  for (size_t i=0; i<specnames.size(); i++) occ[i]=0;
  for (NODEID i=0; i<lf; i++)    
    if (occ[lab[i]]++ == 1) return false;
  return true;
}

bool Dag::belongtoclass(int netclass)
{
  NODEID i=MAXNODEID;
  while (getnodeiter(i)) 
  {   
    int retcnt = 0;
    int chcnt = 0;
    NODEID ic=MAXNODEID;
    while (getchild(i,ic)) 
    { 
      if (ic>=rtstartid)
        retcnt++;
      chcnt++;
    }
    if ((netclass==NET_TREECHILD) && (chcnt>0) && chcnt==retcnt) return false;
    if ((netclass==NET_CLASS1RELAXED) && (retcnt>1)) return false;
  }  
  return true;
}

int Dag::verifychildparent()
{

   if (parent[root]!=MAXNODEID) { 
      cerr << (int)root << ": wrong parent (" << (int)parent[root] << ")of the root ";
      return 1; // wrong parent of the root
   }


   if (lf==1 && nn==1) return 0; // single noded tree; OK

   for (NODEID i = 0; i < nn; i++ ) 
   {
      NODEID p = parent[i];

      if (p==MAXNODEID && i!=root) { 
        cerr << (int)i << ": non-root parent is MAXNODEID";
        return 2; // incorect parent 
      }
      if (p==MAXNODEID && i>=rtstartid) { 
        cerr << (int)i << ": reticulation parent is MAXNODEID";
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

   for (NODEID i = 0; i < nn; i++ ) 
   {
      NODEID p = parent[i];
      if (p==MAXNODEID) continue;      
      if (p<lf) { 
        cerr << (int)i << ": parent " << (int)p << " is a leaf?";
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
   for (NODEID i = rtstartid; i < nn; i++) 
   {
      NODEID rp = retparent[i];
      if (rp == MAXNODEID) { 
        cerr << (int)i << ": retparent[i] is MAXNODEID";
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
    for (NODEID i = 0; i < size(); i++ ) 
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

  		NODEID iparent = MAXNODEID;
  		while (getparentiter(i,iparent))          
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



ostream& Dag::printsubtree(ostream&s, NODEID i, NODEID iparent, int level) 
{ 
	// leaf
  if (level>nn)
  {
    // stop infinite prints
    s << "...loop?" << endl;
    return s;
  }
  // cout << "HERE" << i << endl;
	if (i < lf) return s << species(lab[i]); 

	// internals via getchild	

	if (i>=rtstartid)	
		// reticulation 
		if (iparent == retparent[i]) 
			// convert to leaf and return if coming from retparent[i]
			return   s << spid2retlabel[i];

	s << "(";
	bool first = true;
	NODEID ichild = MAXNODEID;
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

bool Dag::getnodeiter(NODEID &i)
{
  if (i==MAXNODEID) { 
    i=0; 
    return true; 
  }
  return ++i<nn;  
}

// Returns the parents; to get all parents use:
// NODEID p=MAXNODEID;
// while (getparent(i,p)) { .. p is a parent ... }
bool Dag::getparentiter(NODEID i, NODEID &rparent)
{
	if (rparent==MAXNODEID) 
	{
		// first run  		
		rparent = parent[i];
		return rparent!=MAXNODEID;  	
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
// NODEID ichild=MAXNODEID;
// while (getchild(i,ichild)) { .. ichild is a child of i ... }

bool Dag::getchild(NODEID i, NODEID &ichild)
{
	if (i<lf) return false;

	if (ichild==MAXNODEID) 
	{
		// first run  		
		ichild = getleftchild(i); // also retchild[i] 
		return true;
	}

	if (i<rtstartid && ichild==getleftchild(i) && getleftchild(i)!=getrightchild(i))
	{
		// second run; 
		ichild=getrightchild(i);
		return true;		
	}

	return false;
}


void Dag::setexactspecies()
{
    exactspecies = 1;
    for (NODEID i=0; i<lf; i++) 
    	if (lab[i] != i) 
    		{ 
    			exactspecies = 0; 
    			break; 
    		}
}

ostream& Dag::printdebstats(ostream&s) 
{    
    s << " Nodes(nn)=" << (int)size() 
      << " Leaves(lf)=" << (int)lf;
    if (rt)      
        s << " Reticulations(rt)=" << (int)rt << " rtcount()=" << rtcount(); 
    
    s << endl;

    s << " Leaves:0.."<< (int)(lf-1)
      << " TreeNodes:" << (int)lf << ".." << (int)rtstartid-1;
    if (rt) s << " RetNodes:" << (int)rtstartid << ".." << nn-1;

    s << endl;

    if (rt) s << " rtstartid=" << (int)rtstartid;
    s << " root=" << (int)getroot();           
    s << " exactspecies=" << exactspecies;
    return s;
}

ostream& Dag::printdebarrays(ostream&s) 
{
    s << " parent= "; 
    for (NODEID i=0; i<nn; i++) s << " " << (int)i << ":" << parent[i]; 
    s << endl;

    s << " leftchild= "; 
    for (NODEID i=lf; i<rtstartid; i++) s << " " << (int)i << ":" << leftchild[i]; 
    s << endl;
 
    s << " rightchild= "; 
    for (NODEID i=lf; i<rtstartid; i++) s << " " << (int)i << ":" << rightchild[i]; 
    s << endl;

    s << " lab= "; 
    for (NODEID i=0; i<lf; i++) s << " " << int(i) << ":" << lab[i]; 
    s << endl;   

    if (rt)
    {
      s << " retchild= "; 
      for (NODEID i=rtstartid; i<nn; i++) s << " " << (int)i << ":" << retchild[i]; 
      s << endl;

      s << " retparent= "; 
      for (NODEID i=rtstartid; i<nn; i++) s << " " << (int)i << ":" << retparent[i]; 
      s << endl;   

      s << " spid2retlabel= "; 
      for (NODEID i=rtstartid; i<nn; i++) s << " " << int(i) << ":" << spid2retlabel[i]; 
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
	  for (NODEID i = 0; i < size(); i++ ) 
	  {	  		
	  	   	s << setiosflags(ios::left) << setw(2) << (int)i << " ";	  	   	
	  	   		  	   	
    		if (parent[i] == MAXNODEID) 
    			s << "Root";
    		else 
        { 
    			if (i<lf) 
    				s << "Leaf";
    			else if (i>=rtstartid)
    				s << "Reti";
    			else s << "Tree";
    			// s  << " par=" << setw(3) << (int)parent[i];
    		} 	

    		NODEID ic = MAXNODEID;            
    		while (getchild(i,ic)) 
          s << setw(0) << " c=" << (int)ic;    			                  
   		
    		NODEID ip = MAXNODEID;        
    		while (getparentiter(i,ip)) 
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
    //   NODESLOOP if (par[i]==MAXNODEID)
    //   {    
    //        if (cnt!='A') s << ",";          
    //        printsubtree(s,i) << " edgelab=\"" << cnt << "\"" ;
    //        cnt++;
    //   }
    //   s << ")" << " treename=\"" << _gsetrees++ << " " << tn << "\"" ;
    // }
    return s;
  }

Dag::Dag(Dag *d, NODEID v, NODEID p, NODEID w, NODEID q, string retid, double dagweight)
{
    init(d->lf,d->rt+1); 

    for (NODEID i=0; i<d->lf; i++) lab[i] = d->lab[i]; 

#define CID(id)  ((id)>=d->rtstartid)?((id)+1):(id)

    // leaves and tree nodes - parent
    for (NODEID i=0; i<d->rtstartid; i++)            
      parent[i] = CID(d->parent[i]);      
    
    // ret. nodes - parent
    for (NODEID i=d->rtstartid; i<d->nn; i++)        
      parent[i+1]=CID(d->parent[i]);     

    // children of tree nodes
    for (NODEID i=lf; i<d->rtstartid; i++)  
    {          
      leftchild[i] = CID(d->leftchild[i]);
      rightchild[i] = CID(d->rightchild[i]);
    }

    // children of retnodes
    for (NODEID i=d->rtstartid; i<d->nn; i++)            
      retchild[i+1] = CID(d->retchild[i]);

    // children of retnodes
    for (NODEID i=d->rtstartid; i<d->nn; i++)            
      retparent[i+1] = CID(d->retparent[i]);

    // ret labels        
    for (NODEID i=d->rtstartid; i<d->nn; i++)            
      spid2retlabel[i+1] = d->spid2retlabel[i];

    root = d->root;

    parent[root] = MAXNODEID;

    if (!retid.length())    
      retid = "#"+flag_retidprefix + std::to_string(rt);
      

    // insert ret. edge 
    // new retid = nn-1 
    // new treenode id = rtstartid-1

    NODEID nr = nn - 1;
    NODEID nt = rtstartid - 1;

    spid2retlabel[nr] = retid;

    bool rootcase = v==root;

    v=CID(v);
    q=CID(q);
    w=CID(w);

    if (rootcase) // start node above the root
    {
      root = nt;
      parent[root] = MAXNODEID;
      parent[v] = root; 
    }
    else
    {
      p=CID(p);
      parent[nt] = p;

      if (leftchild[p]==v) leftchild[p]=nt;
      else rightchild[p]=nt;

      if (parent[v]==p) parent[v]=nt;
      else retparent[v]=nt;
    }
    
    parent[nr] = nt;
    retparent[nr] = q;

    leftchild[nt] = nr; 
    rightchild[nt] = v;
    retchild[nr] = w;    
    
    if (leftchild[q]==w) leftchild[q]=nr;
    else rightchild[q]=nr;

    if (parent[w]==q) parent[w]=nr;
    else retparent[w]=nr;

    verifychildparent();
    setexactspecies();

}

// Copy constructor (includes shallow copy)
Dag::Dag(const Dag &d, bool shallowcopy)
{
  if (shallowcopy)
  {
    shallow = true;
    lab = d.lab;
    parent = d.parent;
    leftchild = d.leftchild;
    rightchild = d.rightchild;
    retparent = d.retparent;
    retchild = d.retchild;
    spid2retlabel = d.spid2retlabel;
    lf = d.lf;
    rt = d.rt;
    nn = d.nn;
    usagecount = d.usagecount; // TODO: check
    rtstartid = d.rtstartid;
  }
  else
  {
    shallow = false;
    init(d.lf,d.rt);
    memcpy ( lab, d.lab, sizeof(NODEID)*lf );
    memcpy ( parent, d.parent, sizeof(NODEID)*nn );
    memcpy ( leftchild+lf, d.leftchild+lf, sizeof(NODEID)*(nn-lf) );
    memcpy ( rightchild+lf, d.rightchild+lf, sizeof(NODEID)*(nn-lf-rt) );  
    if (rt)
      memcpy ( retparent+rtstartid, d.retparent+rtstartid, sizeof(NODEID)*(nn-rtstartid) );
    for (int i=rtstartid; i<nn; i++)
      spid2retlabel[i]=d.spid2retlabel[i];
  }

  root=d.root;  
  weight=d.weight;     
  exactspecies=d.exactspecies;   

}

Dag::~Dag() 
{
    if (shallow) return; 
    
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



// #define _DEBUG_COMPR_

SPID* Dag::compressedrepr(SPID *r)
{
    if (!r)
        r = new SPID(compressedreprsize());

#ifdef _DEBUG_COMPR_
      printdeb(cout,2);
      cout << "RT=" << rtcount() << endl;
#endif    

    int rpos = 0;

    r[rpos++] = lf;
    r[rpos++] = rtcount();

    NODEID stack[nn];
    int last=0;

    stack[last++] = getroot();

    int speciescnt = specnames.size();

    if (speciescnt+rtcount()>=MAXSPID)
    {
      cerr << "MAXSPID is to small to represent networks. " << speciescnt + rtcount() << endl;
    }

#define PPSTACK cout << "Stack:"; for(int i=0;i<last;i++) cout << stack[i] << " "; cout << endl;     

    // compute mm's for every node
    SPID mm[nn];
    for (NODEID i=0; i<lf; i++) mm[i]=lab[i];
    memset(mm+lf,MAXSPID,sizeof(SPID)*(nn-lf));

    // take ret-mm's 
    pair<NODEID,SPID> retmms[rtcount()]; // storing <reticulationnode, mm from the subtree>

    int retn = 0;
    while (last>0)
    {
        NODEID v = stack[--last];
        if (v>=rtstartid)
        {  
          // reticulation
          NODEID c = getretchild(v);
          
          if (mm[c]==MAXSPID)
          {
              stack[last++] = v;
              stack[last++] = c;
          }
          else          
          {
              mm[v]=speciescnt+mm[c]; // reticulation forest id                            
          }
          
          continue;
        }

        // internal
        NODEID lc = getleftchild(v);
        NODEID rc = getrightchild(v);

        if ((mm[lc]!=MAXSPID) && (mm[rc]!=MAXSPID))
        {
            mm[v] = min(mm[lc],mm[rc]);
            if ((v<=rtstartid) && (mm[v]>=speciescnt))            
              cerr << "Warning: reticulation forest requires tree-child network (use -ot)" << endl;                            
            continue;
        }
        stack[last++] = v;

        if (mm[rc]==MAXSPID) stack[last++] = rc;
        if (mm[lc]==MAXSPID) stack[last++] = lc;            
    }


#ifdef _DEBUG_COMPR_
    // print mm's    
    for (int i=0;i<nn;i++)
    {
      cout << "mm[" << i << "]=$" << (int)mm[i] << " ";
      if (mm[i]<specnames.size()) 
          cout << specnames[mm[i]]; 
      cout << endl;
    }    
#endif

    // aggregate mm's for reticulations
    for(int i=0;i<rtcount();i++)   
    {
      retmms[i].first = i;
      retmms[i].second = mm[getretchild(i)];
    }

    // sort ret-pairs
    sort(retmms, retmms+rtcount(), [](const pair<NODEID,SPID> &pair1, 
                             const pair<NODEID,SPID> &pair2){                              
                              if(pair2.second < pair1.second) return 1;
                              return 0;
                        });
    
    // insert into the stack ret-roots (child of ret. node)
    for(int i=0;i<rtcount();i++)         
       stack[last++]=getretchild(retmms[i].first+rtstartid);       

    // insert root    
    stack[last++] = getroot();

    // gen repr
    while (last>0)
    {
      NODEID v = stack[--last];        
#ifdef _DEBUG_COMPR_
      cout << "Stack:" << int(v) << endl;
#endif      
      if (v==REPRCLOSE)
      {
        r[rpos++]=REPRCLOSE; 
        continue;
      }

      if (v>=rtstartid)
      {         
        // reticulation
        r[rpos++]=speciescnt + mm[getretchild(v)];
        continue;
      }
      if (v<lf)
      {
        // leaf      
        r[rpos++]=lab[v];
        continue;
      }

      r[rpos++]=REPROPEN;
      stack[last++] = REPRCLOSE;
      SPID lc = getleftchild(v);;
      SPID rc = getrightchild(v);
      if (mm[lc]>mm[rc])
      {
        stack[last++] = lc;
        stack[last++] = rc;
      }
      else
      {
        stack[last++] = rc;
        stack[last++] = lc; 
      }      
    }


#ifdef _DEBUG_COMPR_
    for (int i=2;i<rpos;i++)
    {
      if (r[i]==REPROPEN) cout << "[";
      else if (r[i]==REPRCLOSE) cout << "]";
      else if (r[i]<specnames.size()) cout << specnames[r[i]];
      else cout << (int)r[i];
    }
    cout << endl;
#endif    

    // cout << " " << rpos << " " << 3*lf+4*rtcount()-2 << endl;
    // ofstream myfile;
    // myfile.open ("g.dot");
    // myfile << "digraph {" << endl;
    // printdot(myfile,1);
    // myfile << "}" << endl;
    // myfile.close();
  
    return r;
}



bool Dag::istimeconsistent()
{
   vector<vector<NODEID>> adj(nn);
   NODEID e[nn][2];
   NODEID nmap[nn]; 
   NODEID stack[nn];
   int last=0;

   for (NODEID i=0; i<nn; i++)
   {
      e[i][0] = i;
      e[i][1] = MAXNODEID;
      nmap[i] = i;
   }

// #define TCDEBUG

    // glue ret parents
    // 
    for (NODEID i=rtstartid; i<nn; i++)
    {
      NODEID iparent=MAXNODEID;
      getparentiter(i,iparent);
      NODEID p1 = iparent;
      getparentiter(i,iparent);
      NODEID p2 = iparent;
      // cout << " p1=" << p1;            
      // cout << " p2=" << p2;            
      // cout << " mp1=" << e[p1][0] << "." << e[p1][1];            
      // cout << " mp1=" << e[p2][0] << "." << e[p2][1];            
            
      if (p1!=p2)
      {
        if (e[p1][0]==p1) e[p1][0]=p2;
        else e[p1][1]=p2;

        if (e[p2][0]==p2) e[p2][0]=p1;
        else e[p2][1]=p1;    
      }      
      // cout << " amp1=" << e[p1][0] << "." << e[p1][1];            
      // cout << " amp2=" << e[p2][0] << "." << e[p2][1];                  
    }
//#define TCDEBUG

#ifdef TCDEBUG      
    for (NODEID i=0; i<nn; i++)
    {            
      cout << i << "~~" << e[i][0] << " " << e[i][1] << endl;    
    }
#endif    

    // merging parent is not sufficient if network is general
    // fix to min values, when general network is present
    
    for (NODEID i=lf; i<nn; i++)
    {         
       last = 0;
       auto curmap = nmap[i];

       if (curmap<=i) 
       {
          // propagate
          stack[last++] = i;
          nmap[i] = MAXNODEID;

          while (last)
          {
             auto cur = stack[--last];            
             //cout << "STACK" << i <<  " cur=" << cur << endl; 
             if (nmap[cur] != curmap)
             {
                nmap[cur] = curmap;
                stack[last++] = e[cur][0];
                if (e[cur][1]<MAXNODEID) 
                  stack[last++] = e[cur][1];                
             }
          }
          
       }       
    }


#ifdef TCDEBUG      
    for (NODEID i=0; i<nn; i++)
    {            
      cout << i << ". " << e[i][0] << "~~" << e[i][1] << "   " << nmap[i] << endl;    
    }
#endif        
               


     NODEID i = MAXNODEID;    
     while (getnodeiter(i)) 
     { 
        NODEID iparent=MAXNODEID;            
        while (getparentiter(i, iparent))
        {          
            addadjedge(adj, nmap[iparent], nmap[i]);
            // cout << nmap[iparent] << " -> " <<  nmap[i] << endl;
        }
     }
     return isacyclic(adj,nn);
}


void Dag::_getclusters(NODEID v, GTCluster **res, Clusters *clusters)
{
  if (res[v]) return;

  NODEID i=MAXNODEID;
  NODEID *spclu = NULL;
  int cnt = 0;
  while (getchild(v,i)) 
  {
      _getclusters(i, res, clusters);
      if (spclu)  
        spclu = joinspclusters(spclu, res[i]->spcluster, NULL, false);
      else spclu = res[i]->spcluster;
      cnt++;      
  }
  GTCluster *gc = clusters->has(spclu);
  
  if (gc)
  {
    // present, clean if joinspclusters was used
    if (cnt>1) deletespcluster(spclu);      
  }
  else
  {
      // cnt>1, (otherwise cluster was present)
      gc = clusters->add(spclu);
  }
  res[v] = gc;
}

// 0 - not subtree, 1 - subtree, 2 - unknown
// set true if a node is the root of a subtree
bool Dag::_marksubtrees(NODEID v, short *subtree)
{
  
  if (subtree[v]!=2) return subtree[v]==1;
  bool sb = v<rtstartid;
  
  NODEID i=MAXNODEID;
  while (getchild(v,i)) 
  {      
      sb = _marksubtrees(i, subtree) && sb;  // do not optimize!
  }

  if (sb) subtree[v] = 1;
  else subtree[v] = 0;

  return sb;
}

void Dag::getclusters(Clusters *clusters, bool subtreesonly)
{
  GTCluster *res[nn];

  for (NODEID i=0; i<lf; i++) 
  {      
      res[i]=spec2gtcluster[lab[i]];
      clusters->add(res[i]->spcluster); //todo: optimize
  }

  for (NODEID i=lf; i<nn; i++) 
      res[i]=NULL;

  if (subtreesonly)
  {
      short subtree[nn];    
      for (NODEID i=0; i<lf; i++) subtree[i] = 1;
      for (NODEID i=lf; i<nn; i++) subtree[i] = 2;      
      _marksubtrees(getroot(), subtree);

      for (NODEID v=0; v<nn; v++) 
        if (subtree[v]==1)
          _getclusters(v, res, clusters);     
  }
  else
  {
      _getclusters(getroot(), res, clusters);  
  }
}


// Returns true if arg clusters are present in network
bool Dag::hasclusters(Clusters *clusters)
{ 
  bool ok=true;
  Clusters netclusters;  
  getclusters(&netclusters);
  return netclusters.hasall(clusters);
}

// Returns true if arg clusters are present in network as clusters of subtrees
bool Dag::hastreeclusters(Clusters *clusters)
{ 
  bool ok=true;
  Clusters netclusters;  
  getclusters(&netclusters, true);
  // cout <<  "TC" << *this << endl << netclusters << endl;
  // cout << "GC" << *clusters << endl;
  return netclusters.hasall(clusters);
}

void Dag::printclusters(ostream &os)
{
  Clusters netclusters;  
  getclusters(&netclusters);
  os << netclusters;
}
