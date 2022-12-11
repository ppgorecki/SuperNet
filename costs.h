#ifndef _COSTFUN_H
#define _COSTFUN_H

#include "tools.h"

#define COSTDEEPCOAL 1   // deep coalescence (classic)
#define COSTDEEPCOALEDGE 2   // deep coalescence (edge based)
#define COSTDUPLICATION 3    // duplication
#define COSTLOSS 4    // loss
#define COSTDUPLICATIONLOSS 5   // dl
#define COSTROBINSONFOULDS 6   // robinson-foulds



class CostFun
{
	public:
		virtual long computegt(RootedTree &genetree, RootedTree &speciestree)
		{
			speciestree.initlca();
    		speciestree.initdepth();    

    		SPID *lcamap = genetree.getlcamapping(speciestree);

    		double v = compute(genetree, speciestree, lcamap);

		    if (lcamap)
		        delete[] lcamap;

		    return v;
		}

		virtual long compute(RootedTree &genetree, RootedTree &speciestree,  SPID *lcamap)=0;
		virtual COSTT lowerbound(RootedTree &genetree, RootedTree &speciestree) { return 0; }
		virtual COSTT lowerboundnet(RootedTree &genetree, Network &speciestree) { return 0; }

		virtual int costtype() { return 0; }
}; 

class CFDuplication: public CostFun
{
public:

	virtual long compute(RootedTree &genetree, RootedTree &speciestree, SPID *lcamap)
	{      
  		long s = 0;
  		for (SPID i=genetree.lf; i<genetree.nn; i++)        
    		if (lcamap[genetree.leftchild[i]]==lcamap[i] || lcamap[genetree.rightchild[i]]==lcamap[i]) 
            	s++; 
  		return s;
	}
	virtual int costtype() { return COSTDUPLICATION; }
}; 


class CFDeepCoalescence: public CostFun
{

public: 

	long compute(RootedTree &genetree, RootedTree &speciestree, SPID *lcamap)
	{      
	  long s = 0;            
	  for (SPID i=0; i<genetree.nn; i++) 	  
	        if (i!=genetree.root) 
	            s+=speciestree.depth[lcamap[i]]-speciestree.depth[lcamap[genetree.parent[i]]]-1; 
	  return s;
	}

	virtual COSTT lowerbound(RootedTree &genetree, RootedTree &speciestree) 
	{ 
		if (genetree.bijectiveleaflabelling()) return 0;
		return -genetree.nn+1; // if the gene tree contains the same labels; TODO: better estimate
	}

	virtual COSTT lowerboundnet(RootedTree &genetree, Network &net) 
	{ 

		if (genetree.bijectiveleaflabelling()) return 0;
		return -genetree.nn+1; // if the gene tree contains the same labels; TODO: better estimate
	}
	virtual int costtype() { return COSTDEEPCOAL; }
};


class CFDeepCoalescenceEdge: public CostFun
{

public: 

	long compute(RootedTree &genetree, RootedTree &speciestree, SPID *lcamap)
	{      
	  long s = 0;            
	  for (SPID i=0; i<genetree.nn; i++) 	  
	        if (i!=genetree.root) 
	            s+=speciestree.depth[lcamap[i]]-speciestree.depth[lcamap[genetree.parent[i]]]; 
	  return s;
	}

	virtual COSTT lowerbound(RootedTree &genetree, RootedTree &speciestree) 
	{ 
		if (genetree.bijectiveleaflabelling()) return genetree.nn-1;
		return 0; // if the gene tree contains the same labels; TODO: better estimate
	}

	virtual COSTT lowerboundnet(RootedTree &genetree, Network &net) 
	{ 

		if (genetree.bijectiveleaflabelling()) return genetree.nn-1;
		return 0; // if the gene tree contains the same labels; TODO: better estimate
	}
	virtual int costtype() { return COSTDEEPCOALEDGE; }
};

class CFDuplicationLoss: public CostFun
{

	public: 

	long compute(RootedTree &genetree, RootedTree &speciestree, SPID *lcamap) // code repeated for efficiency
	{      
	  long s = 0;            
	  for (SPID i=0; i<genetree.nn; i++) 	  
	  {
	        if (i!=genetree.root) 
	            s+=speciestree.depth[lcamap[i]]-speciestree.depth[lcamap[genetree.parent[i]]]-1;  // dc comp

	        if (i>=genetree.lf && (lcamap[genetree.leftchild[i]]==lcamap[i] || lcamap[genetree.rightchild[i]]==lcamap[i])) // 3*dupl
            	s+=3; 
      }
	  return s;
	}
	virtual int costtype() { return COSTDUPLICATIONLOSS; }
};


class CFLoss: public CostFun
{

	public: 

	long compute(RootedTree &genetree, RootedTree &speciestree, SPID *lcamap) // code repeated for efficiency
	{      
	  long s = 0;            
	  for (SPID i=0; i<genetree.nn; i++) 	  
	  {
	        if (i!=genetree.root) 
	            s+=speciestree.depth[lcamap[i]]-speciestree.depth[lcamap[genetree.parent[i]]]-1;  // dc comp

	        if (i>=genetree.lf && (lcamap[genetree.leftchild[i]]==lcamap[i] || lcamap[genetree.rightchild[i]]==lcamap[i])) // 3*dupl
            	s+=2; 
      }
	  return s;
	}
	virtual int costtype() { return COSTLOSS; }
};



class CFRobinsonFoulds: public CostFun
{

	public: 

	long compute(RootedTree &genetree, RootedTree &speciestree, SPID *lcamap) // code repeated for efficiency
	{      
		cerr << "RF is not implemented yet" << endl;
		exit(-1);	  	
	}
	virtual int costtype() { return COSTROBINSONFOULDS; }
};




// long RootedTree::costrobinsonfoulds(RootedTree &speciestree, SPID *lcamap)
// {
//     //TODO
//     cerr << "RF is not implemented yet" << endl;
//     exit(-1);    
// }







#endif
