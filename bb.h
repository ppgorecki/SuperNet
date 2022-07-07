

#ifndef _BB_H_
#define _BB_H_

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
};

class AdaptiveBB
{

	public: 
		AdaptiveBB() {}
		virtual long init(int rtnumber, COSTT start_cost) { return -1; }
		virtual void branch(long id, SPID rtid) {}

		virtual long start(int rtnumber, int algtype, long parent) { return 0; }
		virtual void stop(long id, COSTT cost) {}
		virtual void costcut(long id, COSTT wrt_cost) {}
		virtual void bestupdated(long id, COSTT wrt_cost) {}		
		virtual void exactsolution(long id) {}	
		virtual void parentcut(long parent, COSTT costcut) {}
		virtual ~AdaptiveBB() {} 
};

 

class AdaptiveBBTree : public AdaptiveBB
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
		void savedot();
		void savetsv();

};


#endif
