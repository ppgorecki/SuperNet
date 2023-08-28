

#include "stats.h"


void NetworkHCStats::print(bool global)
{
    cout << "Cost:" << optcost    
     << " Steps:" << steps 
     << " Climbs:" << improvements 
     << " TopNetworks:" << topnetworks;     

     if (networkclass==NET_TREECHILD)
     cout << " Class:TreeChild";

    if (networkclass==NET_CLASS1RELAXED)
     cout << " Class:Relaxed";

    if (networkclass==NET_GENERAL)
     cout << " Class:General";

    cout << " TimeConsistency:" << timeconsistency;

    if (global) 
     {     		
     		cout 
	     		<< " HCruns:" << startingnets 
	     		<< " HCTime:" << hctime;
	     		// << " MergeTime:" << mergetime;	    // usually low
		   	odtstats.print();
	   }
	   
     cout << endl;

}

// Merge stats. Returns 
//    2 - new networks with the same cost
//    1 - new cost
//    0 - no improvement

// printstats = 0 - no print
// printstats = 1 - improvements
// printstats = 2 - improvements and new nets
// printstats = 3 - always

int NetworkHCStats::merge(NetworkHCStats &nhc, int printstats, bool fullmerge) 
{
	double mtime = gettime();
	int res = 0;

  startingnets++;
  improvements += nhc.improvements;
  steps += nhc.steps;    
  hctime += nhc.hctime;

  odtstats.merge(nhc.getodtstats());

#define PNET  { cout << startingnets << ". ";  nhc.print();  printed=1; }

  bool printed = false;

  if (fullmerge)
  {

    if (!bestdags->size() || nhc.optcost < optcost)
    {
      delete bestdags;
      optcost = nhc.optcost;
      bestdags = nhc.bestdags;
      nhc.bestdags = NULL;
      if (printstats)   
      {    	
        PNET;
        cout << " New optimal cost: " << nhc.optcost << endl;             
      }
      res = 1;
    }
    else 
    	if (nhc.optcost == optcost)
      {      
        int insnets = bestdags->merge(*nhc.bestdags);

        if ((insnets && printstats==2) || printstats==3)
        {
        	PNET;
          cout << " New optimal networks:" << insnets << " " << "Total:" << bestdags->size() << endl;
        }
        if (insnets) res=2;  // new opt nets.            

      }
      else 
      	if (printstats == 3)
  	      PNET;

  }
  else
  {
    if (_newoptimal)
    {    
      if (printstats) 
      {     
        PNET;
        cout << " New optimal cost: " << nhc.optcost << endl;             
      }      
    }
    else {

        if ((_improvements && printstats==2) || printstats==3)
        {
          PNET;
          cout << " New optimal networks:" << _improvements << " " << "Total:" << bestdags->size() << endl;
        }
        if (_improvements) res=2;  // new opt nets.   
    }

  }

  mergetime += gettime()-mtime;
	topnetworks = bestdags->size();

	nhc.improvements = 0;
  nhc.steps = 0;  
  nhc.hctime = 0;

  return res; 
   
}


void NetworkHCStats::savedat(string file, bool labelled)
{
    std::ofstream odtf;
    odtf.open ( file, std::ofstream::out);
    
    if (labelled) odtf << "optcost=";
    odtf << optcost << endl; // cost

  	if (labelled) odtf << "time=";
    odtf << (hctime + mergetime) << endl; // time

  	if (labelled) odtf << "hctime=";
    odtf << hctime << endl; // hill climbing time 

  	if (labelled) odtf << "mergetime=";
    odtf << mergetime << endl; // merge time

  	if (labelled) odtf << "topnets=";
    odtf << topnetworks << endl; // networks

    if (labelled) odtf << "class=";
    odtf << networkclass << endl; // network type

    if (labelled) odtf << "timeconsistency=";
    odtf << timeconsistency << endl; 

  	if (labelled) odtf << "improvements=";
    odtf << improvements << endl; // improvements

  	if (labelled) odtf << "steps=";
    odtf << steps << endl; // steps

    if (labelled) odtf << "bbruns=";
    odtf << odtstats.bbruns << endl; // bbruns

  	if (labelled) odtf << "startnets=";
    odtf << startingnets << endl; // networks merged

    if (labelled) odtf << "memoryMB=";
    odtf << get_memory_size() << endl; // memory size in MB

    if (labelled) odtf << "dtcnt=";
    odtf << odtstats.displaytreecnt << endl; // number of generated display trees
    
    if (labelled) odtf << "naivetime=";
    odtf << odtstats.odtnaivetime << endl; 

    if (labelled) odtf << "naivecnt=";
    odtf << odtstats.odtnaivecnt << endl; 

    if (labelled) odtf << "bbnaivecnt=";
    odtf << odtstats.bbstats.naivecnt << endl; 
    
    if (labelled) odtf << "bbnaivetime=";
    odtf << odtstats.bbstats.naivetime << endl; 
    
    if (labelled) odtf << "bbdpcnt=";
    odtf << odtstats.bbstats.dpcnt << endl; 

    if (labelled) odtf << "bbdptime=";
    odtf << odtstats.bbstats.dptime << endl; 

    if (labelled) odtf << "randseed=";
    odtf << randseed << endl; 

    odtf.close();
}

NetworkHCStats::NetworkHCStats(int _networkclass, int _timeconsistency, DagSet &_visiteddags, unsigned int _randseed, NetworkHCStats *_globalstats) : visiteddags(_visiteddags), randseed(_randseed), globalstats(_globalstats)
{ 
    networkclass= _networkclass;
    timeconsistency = _timeconsistency;

    bestdags = new DagSet();
    improvements = 0; 
    hctime = 0;
    steps = 0;
    startingnets = 0;
    mergetime = 0;
    topnetworks = 0;    
}

NetworkHCStats::~NetworkHCStats()
{ 
    if (bestdags)
      delete bestdags;
}



void NetworkHCStats::setoutfiles(string _outfiles,  bool _odtlabelled)
{ 
    datfile = _outfiles+".dat";
    odtfile = _outfiles+".log";
    odtlabelled = _odtlabelled;
}


void NetworkHCStats::addglobal(Dag *src, double cost)
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


void NetworkHCStats::saveglobal(bool printinfo)
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
      { 
        cout << "Stats data saved: " << datfile << endl;
      }
    }   
}