#ifndef _STATS_H
#define _STATS_H

#include "bb.h"
#include "dagset.h"

#define BB_CUT 4
#define BB_BEST 8
#define BB_INIT 16
#define BB_EXACT 32
#define BB_PARENTCUT 64


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

	DagSet *dagset;

	double optcost;
	long improvements;
	double hctime;
	double mergetime;
	long steps; 
	long startingnets;
	int topnetworks; 
	int timeconsistency;

	int networkclass;
		
	ODTStats odtstats;

public:
 
	NetworkHCStats(int networkclass, int timeconsistency);
	~NetworkHCStats();

	void step() { steps++; }

	int add(Dag &n) { 		
		dagset->add(n);		
		topnetworks = dagset->size();
		return 1; 		
	}

	void setcost(double cost) { 
		improvements++;
		dagset->clear(); 
		optcost = cost;
	}

	void save(string file)
	{
		dagset->save(file);
	}
	

	void finalize()
	{
		hctime = gettime()-hctime;
	}

	void start()
	{
		hctime = gettime();
	}

	ODTStats &getodtstats() { return odtstats; }

	// Save dat file
	void savedat(string file, bool labelled=false);

	// Print stats
	void print(bool global=false);

	// Merge HC results
	int merge(NetworkHCStats &nhc, int printstats);
	

}; 

#endif