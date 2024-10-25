

#include <ostream>
#include <cstdio>
#include "stats.h"
#include "hillclimb.h"


void ClimbStatsBase::print()
{
    cout << " Steps:" << steps 
     << " Climbs:" << improvements;     
}


void ClimbStatsGlobal::print()
{
    if (topnetworks>0)
    {
      cout << "Cost:" << optcost
           << " TopNets:" << topnetworks;          

    ClimbStatsBase::print();

    cout 
        << " CRuns:" << startingnets 
        << " CTime:" << climbtime;
        // << " MergeTime:" << mergetime;     // usually low
        odtstats.print();
    }
    else
    { 
        cout << " TopNets:" << topnetworks;          
    }

}

void ClimbStatsGlobalSampler::print()
{
    ClimbStatsBase::print();

    cout 
        << " CRuns:" << startingnets 
        << " CTime:" << climbtime;
        // << " MergeTime:" << mergetime;     // usually low
        odtstats.print();
}

// Merge stats. Returns 
//    2 - new networks with the same cost
//    1 - new cost
//    0 - no improvement

// printstats = 0 - no print
// printstats = 1 - improvements
// printstats = 2 - improvements and new nets
// printstats = 3 - always

int ClimbStatsBase::merge(ClimbStatsBase &nhc, int printstats, bool fullmerge) 
{
	double mtime = gettime();
	int res = 0;

  startingnets++;
  improvements += nhc.improvements;
  steps += nhc.steps;    
  climbtime += nhc.climbtime;

  odtstats.merge(nhc.getodtstats());

#define PNET  { cout << startingnets << ". ";  nhc.print();  printed=1; }

  bool printed = false;

  if (fullmerge)
  {

    if (!bestdags->size() || nhc.optcost < optcost)
    {
      // New best cost 
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

        // Eq best cost
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
  
  mergetime += gettime()-mtime;
	topnetworks = bestdags->size();

	nhc.improvements = 0;
  nhc.steps = 0;  
  nhc.climbtime = 0;

  return res; 
   
}


int ClimbStatsGlobal::merge(ClimbStatsBase &nhc, int printstats, bool fullmerge) 
{

  int res = ClimbStatsBase::merge(nhc, printstats, fullmerge);

  if (!fullmerge)
  {
    if (_newoptimal)
    {    
      if (printstats) 
      {     
        cout << startingnets << ". ";  
        nhc.print(); 
        cout << " New optimal cost: " << nhc.optcost << endl;             
      }      
    }
    else {

        if ((_improvements && printstats==2) || printstats==3)
        {
          cout << startingnets << ". ";  
          nhc.print(); 
          cout << " New optimal networks:" << _improvements << " " << "Total:" << bestdags->size() << endl;
        }
        if (_improvements) res=2;  // new opt nets.   
    }

  }

  return res;
}

void ClimbStatsBase::_savedat(ostream &odtf, bool labelled, bool partial)
{
    
    if (labelled) odtf << "optcost=";
    odtf << optcost << endl; // cost

  	if (labelled) odtf << "time=";
    odtf << (climbtime + mergetime) << endl; // time

  	if (labelled) odtf << "climbtime=";
    odtf << climbtime << endl; // hill climbing time 

  	if (labelled) odtf << "mergetime=";
    odtf << mergetime << endl; // merge time

  	if (labelled) odtf << "topnets=";
    odtf << topnetworks << endl; // networks

    if (!partial)
    {
      if (labelled) odtf << "class=";
      odtf << networkclass << endl; // network type


      if (labelled) odtf << "timeconsistency=";
      odtf << timeconsistency << endl; 
    }

  	if (labelled) odtf << "improvements=";
    odtf << improvements << endl; // improvements

  	if (labelled) odtf << "steps=";
    odtf << steps << endl; // steps

    if (labelled) odtf << "bbruns=";
    odtf << odtstats.bbruns << endl; // bbruns

  	if (labelled) odtf << "startnets=";
    odtf << startingnets << endl; // networks merged

    if (!partial)
    {
      if (labelled) odtf << "memoryMB=";
      odtf << get_memory_size() << endl; // memory size in MB
    }

    if (labelled) odtf << "dtcnt=";
    odtf << odtstats.displaytreecnt << endl; // number of generated display trees
    
    if (labelled) odtf << "naivetime=";
    odtf << odtstats.odtnaivetime << endl; 

    if (labelled) odtf << "naivecnt=";
    odtf << odtstats.odtnaivecnt << endl; 

    if (!partial)
    {

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
    }
    
}

ClimbStatsBase::ClimbStatsBase(int _networkclass, int _timeconsistency, DagSet *_visiteddagscache, unsigned int _randseed) : visiteddagscache(_visiteddagscache), randseed(_randseed)
{ 
    networkclass= _networkclass;
    timeconsistency = _timeconsistency;
    bestdags = new DagSet();
    improvements = 0; 
    climbtime = 0;
    steps = 0;
    startingnets = 0;
    mergetime = 0;
    topnetworks = 0;    
}


extern int flag_autooutfiles;

void ClimbStatsGlobal::setoutfiles(string _outdirectory, string _outfiles,  bool _odtlabelled)
{ 

    // cout << _outfiles << "---" << _outdirectory << endl;

    if (flag_autooutfiles && _outfiles == ODTBASENAME)
      _outfiles = "";

    if (_outdirectory.length())
    {
      std::filesystem::create_directories(_outdirectory);    
      outdirectory = _outdirectory + filesystem::path::preferred_separator;    
    }
    else
    {
      outdirectory = "";
    }  

    outfile = _outfiles;    

    odtlabelled = _odtlabelled;
    curoutfile = "";    
    firstrun = true;     
}

extern bool flag_hcignorecostgeq;
extern int flag_hcsavefinalodt;
extern float opt_hcignorecostgeq;

// Set new based on cost
bool ClimbStatsGlobal::_checkoutfilename(bool finalfiles)
{
  
  if (flag_hcignorecostgeq && (optcost>=opt_hcignorecostgeq))
  {    
    if (finalfiles==false || flag_hcsavefinalodt==false)
      return false; // ignore
  }
  
  if (flag_autooutfiles)
  { 
    string newoutfile = curoutfile;
    
    char buf[1000];

    if (firstrun || (curoutfilecost != optcost))
    {
      firstrun = false;
      curoutfilecost = optcost;
      sprintf(buf,"%g", optcost);      
      int cnt = 1;

      newoutfile = outdirectory + outfile + buf;

      // find new file name
      while (filesystem::exists(newoutfile+".dat"))
      {
        newoutfile = outdirectory + outfile + buf  + string(".") + to_string(cnt++);            
      }      
    }

    if ((newoutfile != curoutfile) && curoutfile.length())
    {            
      remove((curoutfile+".dat").c_str());

      if (flag_saveextnewick)
        remove((curoutfile+".tre").c_str());
      else remove((curoutfile+".log").c_str());
      
    }

    curoutfile = newoutfile;    
  }
  else
  {      
    curoutfile = outdirectory + outfile;        
  }  

  return true;
}


void ClimbStatsGlobal::addglobal(Dag *src, double cost)
{
    Dag *_;   
    if (!bestdags->size() || optcost > cost)
    {
      setcost(cost);
      bestdags->add(*src, &_);        
      _newoptimal = true;
      _improvements++;      
      if (savewhenimproved)
      {
        if (_checkoutfilename())
        { 
          savedat();
          savebestdags();
        }
      }
    } 
    else if (optcost == cost)
    {
      if (bestdags->add(*src, &_)) _improvements++;           
      if (savewhenimproved)
      {
        if (_checkoutfilename())
        {
          savedat();
          savebestdags();
        }
      }
    }
}


void ClimbStatsGlobal::savebestdags(bool printinfo, bool finalfiles)
{   

    if (_checkoutfilename(finalfiles))      
    { 
      
      ostringstream ss;     
      ss << *bestdags;    
      if (flag_saveextnewick)  
        _save(curoutfile + ".tre", ss.str(), printinfo, "Best networks");      
      else
        _save(curoutfile + ".log", ss.str(), printinfo, "Best networks");      

    }    
}

void ClimbStatsGlobal::savedat(bool printinfo)
{       
    if (curoutfile.length()) 
    { 
      ostringstream ss;
      _savedat(ss, odtlabelled);    
      _save(curoutfile + ".dat", ss.str(), printinfo, "Stats");
    }
}

void ClimbStatsGlobal::_save(string filename, string s, bool printinfo, string sinfo)
{
    ofstream f;
    f.open(filename, ofstream::out);  
    if(!f.is_open() || f.fail())
    {
        cout << "Could not open file " << filename << " to write. Redirecting to cout." << endl;          
        cout << s;
    }
    else
    {
      f << s;    
      f.close();

      if (printinfo) 
      { 
        cout << sinfo << " saved to " << filename << endl;
      }
    }   
}

extern int flag_hcdetailedsummarydat;

void ClimbStatsGlobal::savedatmerged(bool printinfo,  vector<ClimbStatsGlobal*> globalstatsarr, bool finalfiles)
{   

    bool savefiles = _checkoutfilename(finalfiles);

    if (!globalstatsarr.size())
    {
        savedat(printinfo);
        return;
    }

    if (outfile.length()) 
    {

      if (flag_hcdetailedsummarydat)
      {
        ostringstream ss;
        int i = 0;
        for (auto n: globalstatsarr) 
        {
          ss << "Sampler" << i++ << ":" << endl; // class
          ss << "displaytreesampling=" << n->getsampling() << endl; 
          n->_savedat(ss, odtlabelled, true);        
          ss << endl;
        }
        ss << "ExactClimb:" << endl;
        _savedat(ss, odtlabelled, true);              
        ss << endl;


        if (savefiles)
        {
          ofstream datf;
          datf.open(curoutfile+".detailed.dat", std::ofstream::out);
          datf << ss.str();
          datf.close();


          if (printinfo) 
          { 
            cout << "Detailed stats data saved: " << curoutfile << ".detailed.dat" << endl;
          }
        }

      }

      ostringstream ss;
      // Merge 
      for (auto n: globalstatsarr)       
          merge(*n, false, false);
            
      _savedat(ss, odtlabelled);        

      if (savefiles)
      {
        // Merge
        ofstream datf;
        datf.open(curoutfile+".dat", std::ofstream::out);
        datf << ss.str();
        datf.close();

      
        if (printinfo) 
        { 
          cout << "Stats data saved: " << curoutfile << ".dat" << endl;
        }
      }

      // cout << ss;
    }   
}


Dag* ClimbStatsBase::addeq(Dag &n)
{
    Dag *src;
    if (bestdags->add(n, &src)) 
    {           
      // new eq
      eqdags.push_back(src);
    }   
    topnetworks = bestdags->size();
    return src;     
}

Dag* ClimbStatsBase::addnewbest(Dag &n, double cost) 
{ 
    setcost(cost);        
    Dag *src;
    bestdags->add(n, &src);   
    topnetworks = bestdags->size();    
    return src;     
}


