#ifndef _DAG_SET_H
#define _DAG_SET_H

#include "rtree.h"
#include "network.h"
#include "tools.h"


// TODO: better implementation using hash arrays

class DagSet
{
	std::vector<Dag*> *v;
	bool usecounts = 0;
	bool maplabels = 0;

public:
	DagSet(bool counts=false, bool dagshapes=false) 
	{
		maplabels = !dagshapes;
		usecounts = counts;
		v = new std::vector<Dag*>;
	}


	Dag* contains(Dag &n) 
	{		
		for (size_t i=0; i<v->size(); i++)
			if ((*v)[i]->eqdags(&n,maplabels)) return (*v)[i];
		return NULL;
	}

	void insert(Dag *d)
	{
		d->count=1;
		v->push_back(d);
	}

	// insert a copy dag; 
	int add(Dag &n)
	{
		Dag *src = contains(n);
		if (src) 
		{ 
			if (usecounts) src->count++;
			return 0;			
		}		
		insert(new Dag(n));
		return 1;
	}


	// insert a dag (pointer)
	int add(Dag *n)
	{
		Dag *src = contains(*n);
		if (src) 
		{ 
			if (usecounts) src->count++;
			return 0;			
		}
		insert(n);
		return 1;
	}


	size_t size()
	{
		return v->size();
	}

	void clear() {	v->clear();  }

	void save(string file)
	{		
		std::ofstream odtf;
		odtf.open (file, std::ofstream::out);
		odtf << *this; 
		odtf.close();
	}

	// add pointers
	// TODO: counters summing
	size_t merge(DagSet &s) 
	{
		int cnt = 0;
	    for( size_t i = 0; i < s.v->size(); i++ )
	    {
	    	Dag *d = (*(s.v))[i];
	    	if (add(d)) cnt++;
		}
	   	return cnt;
	}

	    
	friend std::ostream &operator<<(std::ostream &os, DagSet const &ds) 
	{ 
		for( size_t i = 0; i < ds.v->size(); i++ )
		{
			Dag *d = (*ds.v)[i];
			if (ds.usecounts) os << d->getcount() << "\t";
	     	os << *d << endl;
	    }
		return os;
	}

};



#endif