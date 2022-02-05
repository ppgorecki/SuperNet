
#include "rtree.h"
#include "tools.h"
#include "network.h"
#include "dag.h"

#define MAXRTNODES 32

typedef long COSTT;
COSTT INFTY = 1000000;

#if  MAXRTNODES < 33
	
	typedef uint64_t RETUSAGE;

	#define emptyretusage(r) r=0
	#define unionretusage(r,a,b) r=((a)|(b))

	#define addleftretusage(r, rtid)  r|=(1<<((rtid)*2)) // via parent
	#define addrightretusage(r, rtid)  r|=(1<<((rtid)*2-1)) // via retparent


#elif MAXRTNODES < 65

	typedef struct { uint64_t lft, rgh; } RETUSAGE;

	#define emptyretusage(r) { r.lft=0; r.rgh=0; }
	#define unionretusage(r,a,b) { r.rgh=a.rgh|b.rgh; r.lft=a.lft|b.lft; }

	#define addleftretusage(r, rtid)  r.lft|=1<<(rtid) // via parent
	#define addrightretusage(r, rtid)  r.rgh|=1<<(rtid) // via retparent

 	ostream& operator<<(ostream& os, const RETUSAGE& r) { return os << r.lft <<"|" << r.rgh; }
#else

 	// TODO!

#endif

#define DELTA 1
#define DELTAUP0 2
#define DELTAUP1 4

#define min3(x,y,z) x < y ? (x < z ? x : z) : (y < z ? y : z)


class DP
{
	RETUSAGE *deltaretusage;
	RETUSAGE *deltaup0retusage;
	RETUSAGE *deltaup1retusage;
	COSTT *delta;
	COSTT *deltaup0;
	COSTT *deltaup1;	
	int8_t *computed;
	RootedTree &genetree;
	Network &network;
	size_t sz;

#define naddr(n) ((n)*(genetree).nn)

public:
	DP(RootedTree &g, Network &n) : genetree(g), network(n), sz(genetree.nn*network.nn) { alloc(); clean(); }


	void preprocess()
	{
		memset(deltaretusage,0,sizeof(RETUSAGE)*sz);
		memset(deltaup0retusage,0,sizeof(RETUSAGE)*sz);
		memset(deltaup1retusage,0,sizeof(RETUSAGE)*sz);	

		// assign leaves compute delta, deltaup0, deltaup1 for leaves
		// todo: optimize leaves arrays
    	for (SPID g=0; g < genetree.lf; g++)    
    	{    		
		    SPID gmap = network.findlab(genetree.lab[g]);
		    COSTT c = 0;
	    	for (SPID s=0; s < network.nn; s++)	
		    {	    		
	    		SPID idx = g+naddr(s);		    
		    	emptyretusage ( deltaretusage[idx] );
		    	if (s==gmap) c = 0;		    			    	
		    	else c = INFTY;		    			    		    	
		    	delta[idx] = c;		    	     		
	     		computed[idx]=DELTA;
	   		}	   			   		
	   	}	   	
	}

	void clean()
	{
		memset(computed,0,sizeof(int8_t)*sz);		
	}

	void print(ostream &c)
	{
		for (SPID g=0; g < genetree.nn; g++)    
	    	for (SPID n=0; n < network.nn; n++)	
	    	{
		    	size_t idx = g+naddr(n);				
				if (computed[idx]&DELTA && delta[idx]<INFTY) 
		    	c << " delta[" << g << " " << n << "]="	<< " " << delta[idx] << "@" << deltaretusage[idx]<< endl;

				if (computed[idx]&DELTAUP0 && deltaup0[idx]<INFTY) 
		    	c << " deltaup0[" << g << " " << n << "]="	<< " " << deltaup0[idx] << "@" << deltaup0retusage[idx]<< endl;
		    
				if (computed[idx]&DELTAUP1 && deltaup1[idx]<INFTY) 
		    	c << " deltaup1[" << g << " " << n << "]="	<< " " << deltaup1[idx] << "@" << deltaup1retusage[idx] << endl;
		    }	    	
	}

	void alloc()
	{
		sz = genetree.nn*network.nn;
		deltaretusage = (RETUSAGE*)malloc(sizeof(RETUSAGE)*sz);
		deltaup1retusage = (RETUSAGE*)malloc(sizeof(RETUSAGE)*sz);
		deltaup0retusage = (RETUSAGE*)malloc(sizeof(RETUSAGE)*sz);		
		delta = (COSTT*)malloc(sizeof(COSTT)*sz);
		deltaup1 = (COSTT*)malloc(sizeof(COSTT)*sz);
		deltaup0 = (COSTT*)malloc(sizeof(COSTT)*sz);
		computed = (int8_t*)malloc(sizeof(int8_t)*sz);			
	}


	~DP() {
		free(deltaretusage);
		free(deltaup1retusage);
		free(deltaup0retusage);
		free(delta);
		free(deltaup0);
		free(deltaup1);
		free(computed);
	}


	inline COSTT _delta(SPID g, SPID n) { 		
		size_t nidx = naddr(n);
		size_t idx = g+nidx;

#ifdef _DPDEBUG_
		cout << "DE:" << g << ":" << n << " c=" << (computed[idx]&DELTA) << endl;
#endif		

		if (computed[idx]&DELTA) return delta[idx];

		
		SPID lft = genetree.leftchild[g];
		SPID rgh = genetree.rightchild[g];
		COSTT resL = _deltaup1(lft,n);
		COSTT resR = _deltaup1(rgh,n);;

		computed[idx]|=DELTA;

		if (resL==INFTY || resR==INFTY)		
			delta[idx] = INFTY;		
		else
		{
			unionretusage( deltaretusage[idx], deltaup1retusage[lft+nidx],  deltaup1retusage[rgh+nidx]);
			delta[idx] = resL+resR;			
		}		
		return delta[idx];
			
	}
		
	inline COSTT _deltaup1(SPID g, SPID n) { 
		
		size_t nidx = naddr(n);
		size_t idx = g+nidx;

#ifdef _DPDEBUG_
		cout << "D1:" << g << ":" << n << " c=" << (computed[idx]&DELTAUP1) << endl;
#endif		

		if (computed[idx]&DELTAUP1) 
			return deltaup1[idx];

		COSTT res;
		RETUSAGE retusage;

		
		if (n>=network.rtstartid || n<network.lf)
		{
			// reticulation or a leaf -> go to delta0
			res = _deltaup0(g,n);
			retusage = deltaup0retusage[idx];			
		}
		else 
		{
			// n is a tree node with two children
			
			SPID s0 = network.leftchild[n];
			SPID s1 = network.rightchild[n];

			int s0reticulation = (s0>=network.rtstartid)?1:0;
			int s1reticulation = (s1>=network.rtstartid)?1:0;

			SPID idx0 = naddr(s0)+g;
			SPID idx1 = naddr(s1)+g;

	        COSTT du0 = _deltaup0(g,s0);
	        COSTT du1 = _deltaup0(g,s1);

	        // ignore edge <s,s0/s1> if s0/s1 is a reticulation
	        res = _delta(g,n);
	        COSTT res0 = 1 - s0reticulation + du0;
	        COSTT res1 = 1 - s1reticulation + du1;                   	
	        res = min3(res, res0, res1);                    
					
	        COSTT sc = 0, optval = 0;
	        bool has = false;
	        RETUSAGE retusagec;

	        if (res == res0)
	        {                    
	            sc = s0;
	            retusagec = deltaup0retusage[idx0];      
	            has = true;
	        }
	        else if (res == res1)
	        {
	            sc = s1;
	            retusagec = deltaup0retusage[idx1];
	            has = true;
	        }
	        else
	        {
	 			retusage = deltaretusage[idx];
	        }

	        // res is min
	        // problem: multiple min's
	        if (has)
	        {
	            // min is achieved by a kid
	            if (sc>=network.rtstartid)
	            {
	                // kid is a reticulation node
	                SPID rtid = network.rtstartid-sc;
	                retusage = retusagec;
	                if (n == network.parent[sc])
	                	addleftretusage( retusage, rtid);                                
	                else
	                	addrightretusage( retusage, rtid);
	            }
	            else
	                retusage = retusagec;
        	}
        }

		deltaup1[idx] = res;
		deltaup1retusage[idx] = retusage;
		computed[idx]|=DELTAUP1;

		return res;

	}
	
	inline COSTT _deltaup0(SPID g, SPID n) 
	{ 
		size_t nidx = naddr(n);
		size_t idx = g+nidx;

#ifdef _DPDEBUG_
		cout << "D0:" << g << ":" << n << " c=" << (computed[idx]&DELTAUP1) << endl;
#endif	

		if (computed[idx]&DELTAUP0) 
			return deltaup0[idx];	

		COSTT res;
		RETUSAGE retusage;

		if (n>=network.rtstartid) // reticulation
    	{
    		SPID sc = network.retchild[n];
    		size_t scidx = g+naddr(sc);	
    		res = _deltaup0(g,sc);
    		retusage = deltaup0retusage[scidx]; 
    		
    		if (sc>=network.rtstartid) // special new case nonTC1
    		{
    			SPID rtid = network.rtstartid-sc;
    			if (network.parent[sc] == n)
    				addleftretusage( retusage, rtid );
    			else
    				addrightretusage( retusage, rtid );	
    		}
    		else
    			res++;   
    	} // ned rtnetwork case
    	else
    	{
    		// treenode cases
			res = _delta(g, n); 
			retusage = deltaretusage[idx];	

			// internal node case
			if (n>=network.lf)
			{
		
	    		// s has 2 children s0 and s1
	    		
				SPID s0 = network.leftchild[n];
				SPID s1 = network.rightchild[n];						

	    		int s0reticulation = (s0>=network.rtstartid)?1:0;
	    		int s1reticulation = (s1>=network.rtstartid)?1:0;

	    		SPID idx0 = naddr(s0)+g;
				SPID idx1 = naddr(s1)+g;			

	        	COSTT du0 = _deltaup0(g,s0);
	        	COSTT du1 = _deltaup0(g,s1);

	        	RETUSAGE retusage0 = deltaup0retusage[idx0]; 
				RETUSAGE retusage1 = deltaup0retusage[idx1]; 

	        	COSTT res0 = du0, res1 = du1;

	        	if (!s0reticulation && !s1reticulation)
	        	{
	        		res0++;
	        		res1++;
	        	}
	        	else if (s0reticulation && s1reticulation)
	        	{
	        		cerr << "Tree child network expected" << endl;
	        		exit(-1);
	        	}

	        	COSTT res01 = min(res0,res1);

	        	if (res>res01)
	        	{

	        		bool used = false, notused = false;
	        		res = res01;

	        		if (res==res0)
	        		{
	        			if (s0reticulation) used = s0;
	        			if (s1reticulation) notused = s1;
	        			retusage = retusage0;
	        		}
	        		else if (res==res1)
	        		{
	        			if (s1reticulation) used = s1;
	        			if (s0reticulation) notused = s0;
	        			retusage = retusage1;
	        		}

	        		if (used)
	        		{        			
						SPID rtid = network.rtstartid-used;
						if (network.parent[used] == n)
							addleftretusage( retusage, rtid);
						else
							addrightretusage( retusage, rtid);	    			
	       
	        		}

	        		if (notused)
	        		{
	        			SPID rtid = network.rtstartid-notused;
						if (network.parent[notused] == n)
							addrightretusage( retusage, rtid);
						else
							addleftretusage(retusage, rtid);	    			
	                }
	        	}
	        } // internal node (if)	    	
    	}

		deltaup0[idx] = res;
	    deltaup0retusage[idx] = retusage;
		computed[idx]|=DELTAUP0;
		return res;
	}

};



double Network::retmindc(RootedTree &genetree)
{

#ifdef _DPDEBUG_	
	cout << " GENE TREE  -------------------------------------- " << endl;
	genetree.printdeb(cout,2);
	cout << " NETWORK ----------------------------------------- " << endl;
	printdeb(cout,2);
	cout << " ----------------------------------------- " << endl;

#endif	

	SPID gnn = genetree.nn;
	SPID glf = genetree.lf;
    
    DP dp(genetree,*this);

    dp.preprocess();

   	COSTT res = INFTY;
	for (SPID s=0; s<nn; s++)	
	{			
		COSTT r = dp._delta(genetree.root, s);
		if (res > r) res = r; 
	}

#ifdef _DPDEBUG_
	dp.print(cout);
#endif	

	return res;
       
}