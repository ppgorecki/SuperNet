
#include "stats.h"

void ODTStats::addbbstats(BBStats &b)
{	
	bbstats.add(b);		
	bbruns++;
}

// Merge
void ODTStats::merge(ODTStats &s)
{	
	bbstats.add(s.bbstats);		
	odtnaivecnt+=s.odtnaivecnt;				
	bbruns+=s.bbruns;
	displaytreecnt+=s.displaytreecnt;
	odtnaivetime+=s.odtnaivetime;
}

void ODTStats::print()
{
	if (bbruns)
		cout << " BB:" << bbruns;     		
	if (odtnaivecnt)
		cout << " Naive:" << odtnaivecnt
		<< " Naivetime:" << odtnaivetime;
	if (displaytreecnt)
		cout << " DTcnt:" << displaytreecnt;
	bbstats.print();		
}

void ODTStats::startnaive()
{
	_naivestart = gettime();    	
}

void ODTStats::stopnaive()
{
	odtnaivecnt++;
	odtnaivetime += gettime() - _naivestart;
}



