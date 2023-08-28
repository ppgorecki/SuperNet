#ifndef _STATS_H
#define _STATS_H

#include "bb.h"
#include "dagset.h"

#define BB_CUT 4
#define BB_BEST 8
#define BB_INIT 16
#define BB_EXACT 32
#define BB_PARENTCUT 64

extern int flag_globaldagcache;
extern int flag_hcsavewhenimproved;

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



class NetworkHCStats
{

	DagSet *bestdags;

	double optcost;
	long improvements;
	double hctime;
	double mergetime;
	long steps; 
	long startingnets;
	int topnetworks; 
	int timeconsistency;
	unsigned int randseed;

	int networkclass;
		
	ODTStats odtstats;

	DagSet &visiteddags;

	vector<Dag*> eqdags;

	NetworkHCStats *globalstats;
	int _improvements;
	bool _newoptimal;
	string odtfile = "";
  	string datfile = "";
  	bool odtlabelled = false;



public:
 
	NetworkHCStats(int networkclass, int timeconsistency, DagSet &visiteddags, unsigned int randseed, NetworkHCStats *_globalstats=NULL);

	~NetworkHCStats();


	void setoutfiles(string _outfiles,  bool _odtlabelled) 
	{ 
		datfile = _outfiles+".dat";
		odtfile = _outfiles+".log";
		odtlabelled = _odtlabelled;
	}

	void step() { steps++; }

	int addeq(Dag &n) 
	{ 		
		Dag *src;
		if (bestdags->add(n, &src))	
		{						
			// new eq
			eqdags.push_back(src);
		}		
		topnetworks = bestdags->size();

		if (globalstats) 
			globalstats->addglobal(src, optcost);
		return 1; 		
	}

	bool inbestdags(Dag &n)
	{		
		return bestdags->contains(n);
	}

	int addnewbest(Dag &n, double cost) 
	{ 
		setcost(cost);		
		Dag *src;
		bestdags->add(n, &src);		
		topnetworks = bestdags->size();

		if (globalstats) 
			globalstats->addglobal(src, optcost);
		return 1; 		
	}

	void saveglobal(bool printinfo = false)
	{		
		if (odtfile.length()) 
		{ 
			save(odtfile);
			if (printinfo) 
			{
				cout << "Optimal networks saved: " << odtfile << endl;
			}
		}

		if (datfile.length()) 
		{ 
			savedat(datfile, odtlabelled);
			if (printinfo) 
				cout << "Stats data saved: " << datfile << endl;
		}		
	}


	void addglobal(Dag *src, double cost)
	{
		Dag *_;		
		if (!bestdags->size() || optcost > cost)
		{
			setcost(cost);
			bestdags->add(*src, &_);				
			_newoptimal = true;
			_improvements++;
			if (flag_hcsavewhenimproved)
				saveglobal();
		} else if (optcost == cost)
		{
			if (bestdags->add(*src, &_)) _improvements++;						
			if (flag_hcsavewhenimproved)
				saveglobal();
		}
	}

	// Set new optimal cost
	// Clear best dags cache
	void setcost(double cost) 
	{ 		
		improvements++;
		eqdags.clear();
		bestdags->clear(); 
		optcost = cost;
	}

	void save(string file)
	{
		bestdags->save(file);
	}
	
	void finalize()
	{
		hctime = gettime()-hctime;
	}

	void start()
	{
		hctime = gettime();
		if (globalstats)
		{
			globalstats->globalstart();
		}
	}

	void globalstart()
	{		
		_improvements = 0;
		_newoptimal = false;
	}

	ODTStats &getodtstats() { return odtstats; }

	// Save dat file
	void savedat(string file, bool labelled=false);

	// Print stats
	void print(bool global=false);

	// Merge HC results
	int merge(NetworkHCStats &nhc, int printstats, bool fullmerge);

	bool alreadyvisited(Dag &n)
	{
		Dag *src;
		return visiteddags.add(&n, &src);
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
	

}; 

#endif
