

#ifndef _DP_H
#define _DP_H

#define DELTA 1
#define DELTAUP0 2
#define DELTAUP1 4

#define min3(x,y,z) x < y ? (x < z ? x : z) : (y < z ? y : z)
#define naddr(n) ((n)*(genetree).nn)

class DP_DCE
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


public:
	DP_DCE(RootedTree &g, Network &n);

	void preprocess();
	void clean();

	void print(ostream &c);

	void alloc();

	~DP_DCE();

	inline RETUSAGE _deltaretusage(SPID g, SPID n);

	COSTT mindeltaroot(RETUSAGE &retusage);

	inline COSTT _delta(SPID g, SPID n);
		
	inline COSTT _deltaup1(SPID g, SPID n);
		
	inline COSTT _deltaup0(SPID g, SPID n);

};

#endif
