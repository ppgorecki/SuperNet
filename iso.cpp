#include <iostream>     // std::cout
#include <algorithm>    // std::next_permutation, std::sort

#include "tools.h"
#include "dag.h"


// Check if internal node can be inserted into isomap
bool Dag::_addisomap(SPID src, SPID dest, Dag *d, SPID *isomap, SPIDPair* cands, int& first,int& last)
{

    if (src==MAXSP && dest==MAXSP) return true;         
    if (src==MAXSP) return false;
    if (dest==MAXSP) return false;
     
    if (isomap[src]!=MAXSP) 
    {
      // already set
      if (isomap[src]==dest) return true; 
      return false; // diff      
    }

#define rtype(n) (((n)<lf)?1:(((n)<rtstartid)?2:3))

    // check types locally
    int srctype = rtype(src);
    int desttype = rtype(dest);
    if (srctype!=desttype) return false;

    // if (src>=rtstartid && dest<rtstartid) return false; 
    // if (src<rtstartid && dest>=rtstartid) return false; 

    // if (src<lf && dest>=lf) return false; 
    // if (src>=lf && dest<lf) return false; 


    SPIDPair csrc[3];    
    int lst=0;


#define inscand(a,b)  { csrc[lst][0]=a; csrc[lst++][1]=b; } 

    // check children
    if (src>=rtstartid)
    {
        // reticulation
        if (isomap[retchild[src]]!=MAXSP)
        {
            if (isomap[retchild[src]]!=d->retchild[dest]) return false; // wrong retchild            
        }     
        else
        {
          if (rtype(retchild[src])!=rtype(d->retchild[dest]))
            return false;
          // children must be assigned 
          inscand(retchild[src],d->retchild[dest]);          
        }   
    }    
    else if (src>=lf)
      {
        // tree node
        SPID dl = isomap[leftchild[src]];
        SPID dr = isomap[rightchild[src]];

        int tsl = rtype(leftchild[src]);
        int tsr = rtype(rightchild[src]);
        int tdl = rtype(d->leftchild[dest]);
        int tdr = rtype(d->rightchild[dest]);

        bool ll = (tsl==tdl) && (dl==d->leftchild[dest]);
        bool lr = (tsl==tdr) &&  (dl==d->rightchild[dest]);
        bool rl = (tsr==tdl) &&  (dr==d->leftchild[dest]);
        bool rr = (tsr==tdr) &&  (dr==d->rightchild[dest]);

        if (dl!=MAXSP && dr!=MAXSP)
        {
            // two variants
            if  (!((ll && rr) || (lr && rl))) return false; 
        }
        else if (dl!=MAXSP)
          {            
          
            // dr == isomap[rightchild[src]] == MAXSP
            if (ll && (tsr==tdr))
              { inscand(rightchild[src],d->rightchild[dest]); }          
            else if (lr && (tsr==tdl))
              { inscand(rightchild[src],d->leftchild[dest]); }          
            else return false;          
          } 
        else if (dr!=MAXSP)
          {          
            // isomap[leftchild[src]] == MAXSP
            if (rr && (tsl==tdl))
              { inscand(leftchild[src],d->leftchild[dest]); }          
            else if (rl && (tsl==tdr))
              { inscand(leftchild[src],d->rightchild[dest]); }         
            else return false;
          }
        else
        {
            //check only types
            if (!(((tsl==tdl) && (tsr==tdr)) || ((tsr==tdl) && (tsl==tdr))))
              return false;

        }

      }
    else
      {
        // leaf case
        // skip
      }

    // check parents
    if (root==src && d->root==dest) return true; // root case
    if (root==src) return false;
    if (d->root==dest) return false;

    // non-root case
    if (src>=rtstartid)
    {
        // ret node
        SPID dl = isomap[parent[src]];
        SPID dr = isomap[retparent[src]];


        // types leaf/internal/ret
        int tsl = rtype(parent[src]);
        int tsr = rtype(retparent[src]);
        int tdl = rtype(d->parent[dest]);
        int tdr = rtype(d->retparent[dest]);

        bool ll = (tsl==tdl) && (dl==d->parent[dest]);
        bool lr = (tsl==tdr) && (dl==d->retparent[dest]);
        bool rl = (tsr==tdl) && (dr==d->parent[dest]);
        bool rr = (tsr==tdr) && (dr==d->retparent[dest]);        

        if (dl!=MAXSP && dr!=MAXSP)
        {
            if  (!((ll && rr) || (lr && rl))) return false; 
        }
        else if (dl!=MAXSP)
          {
            if (ll && (tsr==tdr))
              { inscand(retparent[src],d->retparent[dest]); }          
            else if (lr && (tsr==tdl))
              { inscand(retparent[src],d->parent[dest]); }          
            else return false;          
          } 
        else if (dr!=MAXSP)
          {
            if (rr && (tsl==tdl))
              { inscand(parent[ src],d->parent[dest]); }          
            else if (rl && (tsl==tdr))
              { inscand(parent[src],d->retparent[dest]); }         
            else return false;        
          }
        else
        {
            //check only types
            if (!(((tsl==tdl) && (tsr==tdr)) || ((tsr==tdl) && (tsl==tdr))))
              return false;        
        }

    }
    else 
    {
      // tree node/leaf
      SPID pr = isomap[parent[src]];
      if (pr!=MAXSP)
        { if (d->parent[dest]!=pr) return false; }
      else
      {
         if (rtype(parent[src])!=rtype(d->parent[dest]))
            return false;
         inscand(parent[src],d->parent[dest]); 
      }
    }    


    if (lst && cands)
    {
      //copy candidates      
      for (int i=0;i<lst;i++)
       {
         cands[last][0] = csrc[i][0];
         cands[last++][1] = csrc[i][1];         
       }
       cout << endl;
    }

    return true;    
}

#define ppisomap() cout << "{"; for (SPID q=0;q<nn;q++) cout << isomap[q] << " "; cout << "}" << endl;

// only tree child and class 1 testing (to be proved)
bool Dag::eqdagstc1(Dag *d)
{
  if (d->lf!=lf || d->rt!=rt || d->nn!=nn) return false;

  SPID isomap[nn];
  SPIDPair cands[2*nn];
  int cfirst=0, clast=0;

  for (SPID i=0; i<nn; i++) isomap[i]=MAXSP;

 
  // leaves
  for (SPID i=0; i<lf; i++) 
    {
      isomap[i]=d->findlab(lab[i]);
      if (isomap[i]==MAXSP) return false;    
      cands[clast][0]=parent[i];
      cands[clast++][1]=d->parent[isomap[i]];
    }

  // root 
  if (!_addisomap(root,d->root,d,isomap,cands,cfirst,clast))   
  		return false;
  
  
  // inserted=0
  while (cfirst<clast)
  { 
    SPID src = cands[cfirst][0];
    SPID dest = cands[cfirst++][1];        

    if (!_addisomap(src,dest,d,isomap,cands,cfirst,clast))   
    {
#ifdef _DEBUG_ISO_      
    	cout << "Cannot set: " << src << "=" << dest << endl;
  		ppisomap();
#endif      
  		return false;
  	}
    isomap[src] = dest; // set new map
                     
  }

#ifdef _DEBUG_ISO_
  ppisomap();
#endif  

  if (validisomap(d,isomap))
    return true;

  // map maybe incomplete

  return false; // unknown
}


bool Dag::eqdagsbypermutations(Dag *d)
{

  if (d->lf!=lf || d->rt!=rt || d->nn!=nn) return false;

  SPID isomap[nn];
  for (SPID i=0; i<nn; i++) isomap[i]=MAXSP;  

  // leaves
  for (SPID i=0; i<lf; i++) 
    {
      isomap[i]=d->findlab(lab[i]);
      if (isomap[i]==MAXSP) return false;          
    }

  isomap[root] = d->root;

  return _eqdagsbypermutations(d,isomap);

}


bool Dag::eqdags(Dag *d, bool maplabels)
{

  if (d->lf!=lf || d->rt!=rt || d->nn!=nn) return false;

  SPID isomap[nn];
  for (SPID i=0; i<nn; i++) isomap[i]=MAXSP;  

  // leaves
  if (maplabels)
    for (SPID i=0; i<lf; i++) 
      {
        isomap[i]=d->findlab(lab[i]);
        if (isomap[i]==MAXSP) return 0;          
      }

  int _;

  if (!_addisomap(root,d->root,d,isomap,NULL,_,_)) 
  {
#ifdef _DEBUG_ISO_  
  		cout << "Cannot set roots" << endl;
  	  ppisomap();
#endif      
  		return false;
  }

  return _eqdags(d,isomap);

}


bool Dag::validisomap(Dag *d, SPID *isomap)
{

  for (SPID i=lf;i<rtstartid;i++)
  {
      //tree node

      SPID di = isomap[i];      
      if (di==MAXSP) return false; // prog error?      
      if (parent[i]==MAXSP && d->parent[di]!=MAXSP) return false;      
      if ((parent[i]!=MAXSP) && isomap[parent[i]]!=d->parent[di]) return false;              

      SPID lt=isomap[leftchild[i]];
      SPID rt=isomap[rightchild[i]];
      SPID ld=d->leftchild[di];
      SPID rd=d->rightchild[di];
      
      if (!((lt==ld && rt==rd) || (rt==ld && lt==rd))) return false;       
  }


  for (SPID i=rtstartid;i<nn;i++)
  {
      //tree node
      SPID di = isomap[i];
      if (di==MAXSP) return false; // prog error?      

      if (isomap[retchild[i]]!=d->retchild[di]) return false;       

      SPID lt=isomap[parent[i]];
      SPID rt=isomap[retparent[i]];
      SPID ld=d->parent[di];
      SPID rd=d->retparent[di];

      if (!((lt==ld && rt==rd) || (rt==ld && lt==rd))) return false;  
  }

  return true;
  

}


// infer maps by permutations (costly)
bool Dag::_eqdagsbypermutations(Dag *d, SPID *isomap)
{

  int pos[nn], dest[nn];
  bool used[nn];
  int cnt=0;
   
  for (int i=0;i<nn;i++) used[i]=false;

  for (int i=0;i<nn;i++)  
    if (isomap[i]==MAXSP) pos[cnt++]=i;
    else used[isomap[i]]=true;
  

  int dcnt=0;
  for (int i=0;i<nn;i++)  
    if (!used[i]) dest[dcnt++]=i;
  
  if (cnt!=dcnt)
  {
    cerr << "Wrong array sizes in _eqdagsbruteforce: " << cnt << " vs. " << dcnt << endl;
    exit(-1);
  }

  std::sort (dest,dest+cnt);
 
  do {
    
      for (int i=0;i<cnt;i++)
        isomap[pos[i]] = dest[i];
    
      if (validisomap(d,isomap)) 
        {
#ifdef _DEBUG_ISO_          
          ppisomap();
#endif          
          return true;
        }    

   } while ( std::next_permutation(dest,dest+cnt) );
    
  return false;
}

// infer maps by combinations (less costly)
bool Dag::_eqdags(Dag *d, SPID *isomap)
{

  int pos[nn], dest[nn];
  bool used[nn];
  int cnt=0;

  for (int i=0;i<nn;i++) used[i]=false;

  for (int i=0;i<nn;i++)  
    if (isomap[i]==MAXSP) pos[cnt++]=i;
    else used[isomap[i]] = true;  

  int dcnt=0;
  for (int i=0;i<nn;i++)  
    if (!used[i]) dest[dcnt++]=i;
  
  if (cnt!=dcnt)
  {
    cerr << "Wrong array sizes in _eqdagsbruteforce: " << cnt << " vs. " << dcnt << endl;
    exit(-1);
  }

  SPID candlist[nn][cnt];
  SPID candlistsizes[nn];  

  int _;
  int fixed=1;
  int left = cnt;

  while (fixed)
  {

#ifdef _DEBUG_ISO_    
  	cout << "=================" << endl;
#endif    
  	fixed=0;
  	SPID fixedarr[cnt];
	  for (int i=0;i<cnt;i++)
	  {
	    SPID cur = pos[i];
	    if (cur==MAXSP) continue;

	    int s=0;

      // build candidates
	    for (int j=0; j<cnt; j++)	         
	      if (dest[j]!=MAXSP && _addisomap(cur,dest[j],d,isomap,NULL,_,_))
	        	candlist[cur][s++]=j;
	    
	    if (!s) // no candidates; non-identical dags 
	    {
#ifdef _DEBUG_ISO_        
	    	cout << "No cands: " << cur << " :::";
	    	ppisomap();
#endif        
	    	return 0; 
	    }

	    if (s==1) // one candidate -> fix
	    {        
	    	isomap[cur] = dest[candlist[cur][0]];

#ifdef _DEBUG_ISO_        
	    	cout <<"FIx!" << cur << "==>" << isomap[cur] << endl;
#endif        
	    	dest[candlist[cur][0]]=MAXSP;
	    	fixed++;
	    	pos[i] = MAXSP; 
        left--;
	    }

	    candlistsizes[cur]=s;


#ifdef _DEBUG_ISO_        
	    cout << cur << ":len=" << s << " :";
	    for (int j=0; j<s; j++)	    	     
	      cout << dest[candlist[cur][j]] << " ";
#endif      
	  	
	  }
  }

  if (left>0)
  {
    // recursive  

    // find cand to 
    int v=0;
    for (;v<cnt;v++)      
      if (pos[v]!=MAXSP) break;

    SPID cur = pos[v];    

    for (int i=0;i<candlistsizes[cur];i++)
    {

      SPID isomap2[nn];    
      memcpy ( isomap2, isomap, sizeof(SPID)*nn );
      isomap2[cur] = dest[candlist[cur][i]]; // fix 

#ifdef _DEBUG_ISO_
      cout << " Set: " << cur << ":=" << isomap2[cur];
#endif      

      if (_eqdags(d, isomap2)) return true;    
    }

    return false;

  }

  return validisomap(d,isomap);
}




