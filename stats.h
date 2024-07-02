#ifndef _STATS_H
#define _STATS_H

#include "bb.h"
#include "dagset.h"

#include <filesystem>



#define BB_CUT 4
#define BB_BEST 8
#define BB_INIT 16
#define BB_EXACT 32
#define BB_PARENTCUT 64

extern int flag_globaldagcache;
extern int flag_savewhenimproved;

struct bbnode {
    int rtnumber;
    int type;
    long parent;
    COSTT cost;
    COSTT wrt_cost;
    double stime;
    bool visitedchild;
};


typedef struct BBStats
{
	double naivetime;
	long naivecnt;
	double dptime;
	long dpcnt;

	BBStats()
	{
		naivetime=0;
		naivecnt=0;
		dptime=0;
		dpcnt=0;
	}

	void add(struct BBStats &b);
	void print();	

} BBStats;

class BBTreeStats
{
	protected:		
		vector<bbnode> v;				

	public: 			
		
		BBStats stats;

		int minrtnumber;

		long init(int rtnumber, COSTT start_cost);
		long start(int rtnumber, int algtype, long parent);
		void stop(long id, COSTT cost);		
		void costcut(long id, COSTT wrt_cost);
		void bestupdated(long id, COSTT wrt_cost);
		void exactsolution(long id);		
		void parentcut(long parent, COSTT costcut);
		bool visitedchild(long parent);
		void savedot();
		void savetsv();

};



typedef struct ODTStats
{
	long displaytreecnt;
	double odtnaivetime; 
	long odtnaivecnt; // total
	long bbruns;

	double _naivestart;

	BBStats bbstats;
	
	ODTStats() {
		displaytreecnt=0;
		odtnaivetime=0;
		odtnaivecnt=0;		
		bbruns=0;
	}

	// After BB run
	void addbbstats(BBStats &b);

	// Merge
	void merge(ODTStats &s);
	void print();
	void startnaive();
	void stopnaive();


} ODTStats;

class ClimbStatsGlobal; 


class ClimbStatsBase
{
public:

	// Save to dat file
	void _savedat(ostream &odtf, bool labelled=false, bool partial=false);

	DagSet *bestdags;

	double optcost;
	long improvements;
	double climbtime;
	double mergetime;
	long steps; 
	long startingnets;
	int topnetworks; 
	int timeconsistency;
	unsigned int randseed;

	int networkclass;
		
	ODTStats odtstats;

	DagSet *visiteddagscache;

	vector<Dag*> eqdags;

	int _improvements;
	bool _newoptimal;
  	bool odtlabelled = false;


public:
 
	ClimbStatsBase(int networkclass, int timeconsistency, DagSet *visiteddagscache, unsigned int randseed);

	~ClimbStatsBase()
	{ 
    	if (bestdags)
      	delete bestdags;
	}

	DagSet *getbestdags() { return bestdags; }

	void setoutfiles(string _outdirectory, string _outfiles,  bool _odtlabelled);

	void step() { steps++; }

	Dag* addeq(Dag &n);

	bool inbestdags(Dag &n)
	{		
		return bestdags->contains(n);
	}

	Dag* addnewbest(Dag &n, double cost);

	// Set new optimal cost
	// Clear best dags cache
	void setcost(double cost) 
	{ 		
		improvements++;
		eqdags.clear();
		bestdags->clear(); 
		optcost = cost;
	}
	
	void finalize()
	{
		climbtime = gettime()-climbtime;
	}


	ODTStats &getodtstats() { return odtstats; }

	void start()
	{
		 climbtime = gettime();
	}

	// Print stats
	virtual void print();

	// Merge HC results
	int merge(ClimbStatsBase &climbstats, int printstats, bool fullmerge);

	bool alreadyvisited(Dag &n)
	{
		Dag *src;
		return visiteddagscache->add(&n, &src);
	}

	bool haseqdags()
	{				
		return eqdags.size()>0;
	}

	Dag *popeqdag()
	{			
		Dag *dag = eqdags.back();		
		eqdags.pop_back();
		return dag;
	}

	int gettimeconsistency()
	{ return timeconsistency; }

	unsigned int getrandseed()
	{ return randseed; }

	int getnetworkclass()
	{ return networkclass; }


	

	friend class ClimbStatsGlobal;

	void printnetworkinfo()
	{
		if (networkclass==NET_TREECHILD)
	      cout << " Class:TreeChild";

	    if (networkclass==NET_CLASS1RELAXED)
	      cout << " Class:Relaxed";

	    if (networkclass==NET_GENERAL)
	      cout << " Class:General";

	    cout << " TimeConsistency:" << timeconsistency;
	}
	

}; 

class ClimbStatsGlobal: public ClimbStatsBase
{

	protected:

		string outfile = "";  		  	
		string outdirectory = "";  		  	
		string curoutfile = "";
		double curoutfilecost = -100000;
		bool firstrun = true;

  		void _save(string filename, string s, bool printinfo, string sinfo="");
		virtual bool _checkoutfilename(bool finalfiles=false);

		bool savewhenimproved;

	public: 

	ClimbStatsGlobal(int networkclass, int timeconsistency, DagSet *visiteddagscache, unsigned int randseed, bool savewhenimproved) : 
	ClimbStatsBase(networkclass, timeconsistency, visiteddagscache, randseed), savewhenimproved(savewhenimproved) {}

	void globalstart()
	{		
		_improvements = 0;
		_newoptimal = false;
	}


	void savebestdags(bool printinfo = false, bool finalfiles=false);
	
	void savedat(bool printinfo = false);

	void addglobal(Dag *src, double cost);

	void setsavewhenimproved(bool _savewhenimproved) { savewhenimproved = _savewhenimproved; }

	int merge(ClimbStatsBase &nhc, int printstats, bool fullmerge);

	void savedatmerged(bool printinfo,  vector<ClimbStatsGlobal*> globalstatsarr, bool finalfiles);
	

	void setoutfiles(string _outdirectory, string _outfiles,  bool _odtlabelled);

	virtual float getsampling() { return 0; }

	virtual bool issampler() { return false; }

	virtual int testsampling(int reticulationcnt) { return 1 << reticulationcnt; }
	
	virtual ostream& info(ostream &os) { return os << "exact climbs"; }

	void print();


};

class ClimbStatsGlobalSampler: public ClimbStatsGlobal
{
	protected:
		float displaytreesampling;

		bool _checkoutfilename(bool finalfiles=false) { return false; }

	public:
		ClimbStatsGlobalSampler(int networkclass, int timeconsistency, DagSet *visiteddagscache, unsigned int randseed, float sampling) 
		: ClimbStatsGlobal(networkclass, timeconsistency, visiteddagscache, randseed, false), displaytreesampling(sampling)
	{}

	void saveglobal(bool printinfo = false) {}

	bool issampler() { return true; }

	float getsampling() { return displaytreesampling; }

	int testsampling(int reticulationcnt)
	{
		random_device rd;
		mt19937 gen(rd());

		exponential_distribution<> expdistr(displaytreesampling);

		int tid = 0;
		int cnt = 0;
		int maxtid = 1 << reticulationcnt;
		float sampledtid = 0;
		

    	while(tid<maxtid)
		{
			int nexttid=(int)floor(tid+expdistr(gen));
    		if (nexttid==tid) tid++;
    		else tid = nexttid;
			cnt++;
    	}
    	return cnt;

	}

	ostream& info(ostream &os) { return os << "sampler climbs " << displaytreesampling; }

	void print();


  	
};


class ClimbStats: public ClimbStatsBase
{

	protected:
		ClimbStatsGlobal *globalstats;
	
	public:
		ClimbStats(DagSet &_visiteddagscache, ClimbStatsGlobal*_globalstats) :  ClimbStatsBase(_globalstats->getnetworkclass(), _globalstats->gettimeconsistency(), visiteddagscache, _globalstats->getrandseed()), globalstats(_globalstats)
		{}		

		void start()
		{
			  ClimbStatsBase::start();  
  			  globalstats->globalstart();
		}

		int addnewbest(Dag &n, double cost) 
		{     
    		Dag *src = ClimbStatsBase::addnewbest(n, cost);
    		globalstats->addglobal(src, optcost);
    		return 1;     
		}

		int addeq(Dag &n)
		{
			Dag * src = ClimbStatsBase::addeq(n);
			globalstats->addglobal(src, optcost);
			return 1;
		}		

};


#endif
