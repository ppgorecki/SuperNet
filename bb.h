

#ifndef _BB_H_
#define _BB_H_

#if !defined USE_PRIORITY_QUEUE_MINRT && !defined USE_PRIORITY_QUEUE_MINLB && !defined USE_QUEUE_BFS
#define USE_PRIORITY_QUEUE_MINLB
#endif


#include "tools.h"

#define ALG_DP 1
#define ALG_NAIVE 2

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

class BBTreeStats
{
	protected:		
		vector<bbnode> v;				

	public: 			
		
		double algnaivetime;
		double algdptime;
		long algnaivecnt;
		long algdpcnt;
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


#endif
