#ifndef _DAG_SET_H
#define _DAG_SET_H

#include "rtree.h"
#include "network.h"
#include "tools.h"

// TODO: better implementation using hash arrays
// TODO: usagecount improve

class DagSet
{

protected:

	std::vector<Dag*> *v;
	bool usecounts = 0;  // --uniquedags (see)
	bool dagshapes = 0;

public:

	DagSet(bool counts=false, bool _dagshapes=false) 
	{
		dagshapes = _dagshapes;
		usecounts = counts;
		v = new std::vector<Dag*>;
	}

	std::vector<Dag*>::iterator begin() 
	{
    	return v->begin();
	}

	std::vector<Dag*>::iterator end() 
	{
    	return v->end();
	}


	Dag* contains(Dag &n) 
	{		
		for (size_t i=0; i<v->size(); i++)
			if ((*v)[i]->eqdags(&n, !dagshapes)) return (*v)[i];
		return NULL;
	}

	void insert(Dag *d)
	{
		d->usagecount++;
		v->push_back(d);
	}

	// insert a copy dag; 
	int add(Dag &n, Dag **retsrc)
	{
		Dag *src = contains(n);
		if (src) 
		{ 
			*retsrc = src;
			if (usecounts) src->usagecount++;
			return 0;			
		}		
		Dag *d = new Dag(n);
		*retsrc = d;
		d->usagecount = 0;
		insert(d);
		return 1;
	}


	// insert a dag (pointer)
	int add(Dag *n, Dag **retsrc)
	{
		Dag *src = contains(*n);
		if (src) 
		{ 
			*retsrc = src;
			if (usecounts) src->usagecount++;
			return 0;			
		}
		*retsrc = n;
		insert(n);
		return 1;
	}


	size_t size()
	{
		return v->size();
	}

	void clear() {	

		for( size_t i = 0; i < v->size(); i++ )
		{
	    	delete (*v)[i];
		}	
		v->clear();  
	}

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
	    	Dag *src;
	    	if (add(d, &src)) cnt++;
		}
	   	return cnt;
	}

	    
	friend std::ostream &operator<<(std::ostream &os, DagSet const &ds) 
	{ 
		for( size_t i = 0; i < ds.v->size(); i++ )
		{
			Dag *d = (*ds.v)[i];

			d->verifychildparent();
			if (ds.usecounts) os << d->getusagecount() << "\t";
	     	os << *d << endl;
	    }
		return os;
	}

};



#endif