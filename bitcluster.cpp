
#include "bitcluster.h"

#if MAXSPECIES==64 || MAXSPECIES==128

#define _2p64  (bitcluster)0x8000000000000000*2

bitcluster bcsingleton[MAXSPECIES] = {
	0x1,0x2,0x4,0x8,
	0x10,0x20,0x40,0x80,
	0x100,0x200,0x400,0x800,
	0x1000,0x2000,0x4000,0x8000,
	0x10000,0x20000,0x40000,0x80000,
	0x100000,0x200000,0x400000,0x800000,
	0x1000000,0x2000000,0x4000000,0x8000000,
	0x10000000,0x20000000,0x40000000,0x80000000,
	0x100000000,0x200000000,0x400000000,0x800000000,
	0x1000000000,0x2000000000,0x4000000000,0x8000000000,
	0x10000000000,0x20000000000,0x40000000000,0x80000000000,
	0x100000000000,0x200000000000,0x400000000000,0x800000000000,
	0x1000000000000,0x2000000000000,0x4000000000000,0x8000000000000,
	0x10000000000000,0x20000000000000,0x40000000000000,0x80000000000000,
	0x100000000000000,0x200000000000000,0x400000000000000,0x800000000000000,
	0x1000000000000000,0x2000000000000000,0x4000000000000000,0x8000000000000000
#if MAXSPECIES==128
	,
	_2p64*0x1,_2p64*0x2,_2p64*0x4,_2p64*0x8,
	_2p64*0x10,_2p64*0x20,_2p64*0x40,_2p64*0x80,
	_2p64*0x100,_2p64*0x200,_2p64*0x400,_2p64*0x800,
	_2p64*0x1000,_2p64*0x2000,_2p64*0x4000,_2p64*0x8000,
	_2p64*0x10000,_2p64*0x20000,_2p64*0x40000,_2p64*0x80000,
	_2p64*0x100000,_2p64*0x200000,_2p64*0x400000,_2p64*0x800000,
	_2p64*0x1000000,_2p64*0x2000000,_2p64*0x4000000,_2p64*0x8000000,
	_2p64*0x10000000,_2p64*0x20000000,_2p64*0x40000000,_2p64*0x80000000,
	_2p64*0x100000000,_2p64*0x200000000,_2p64*0x400000000,_2p64*0x800000000,
	_2p64*0x1000000000,_2p64*0x2000000000,_2p64*0x4000000000,_2p64*0x8000000000,
	_2p64*0x10000000000,_2p64*0x20000000000,_2p64*0x40000000000,_2p64*0x80000000000,
	_2p64*0x100000000000,_2p64*0x200000000000,_2p64*0x400000000000,_2p64*0x800000000000,
	_2p64*0x1000000000000,_2p64*0x2000000000000,_2p64*0x4000000000000,_2p64*0x8000000000000,
	_2p64*0x10000000000000,_2p64*0x20000000000000,_2p64*0x40000000000000,_2p64*0x80000000000000,
	_2p64*0x100000000000000,_2p64*0x200000000000000,_2p64*0x400000000000000,_2p64*0x800000000000000,
	_2p64*0x1000000000000000,_2p64*0x2000000000000000,_2p64*0x4000000000000000,
	_2p64*0x8000000000000000
	
#endif

};


#endif

std::ostream &ppspcluster(std::ostream &os, bitcluster c)
{	
	bool first=true;
	for (int i=MAXSPECIES;i--;i>=0)
	{
		if (SETHAS(c,i))
		{	
			os << '1';
			first=false;
		}
		else 
			os << '.';
	}
	return os;

}


std::ostream &ppbitclusterspecies(std::ostream &os, bitcluster c)
{	
	for (int i=0;i<MAXSPECIES;i++)
	
		if (SETHAS(c,i))		
			os << specnames[i];					
	
	return os;	
}

void clustergraphtester()
{
	for (int i=0; i<MAXSPECIES; i++)
	{
		bitcluster c = EMPTYSET();
		c = SETADD(c, i);	
		cout << i << ":";
		ppspcluster(cout, c);
		cout << endl;
	}

	for (int i=0; i<MAXSPECIES; i++)
		for (int j=i+1; j<MAXSPECIES; j++)
	{
		bitcluster c = EMPTYSET();
		c = SETADD(c, i);	
		c = SETADD(c, j);			
		ppspcluster(cout, c);
		cout << endl;		
	}
	
}