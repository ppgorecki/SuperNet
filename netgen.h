#ifndef _RANDNETS_
#define _RANDNETS_

#include "network.h"

Network *addrandreticulations(int reticulations, Network *n, int networktype, int timeconsistency, bool uniform, Clusters *guideclusters, Clusters *guidetree, bool keepsource=false);

class NetGen
{
  protected:
    long current;
    long count;
    virtual Network *get() { return NULL; }
  
  public:
    NetGen(long count=-1) : current(-1), count(count) { reset(); }

    virtual void reset()
    { current=-1; }

    Network *next()
    {
       current++;
       if (count<0) return get(); // always
       if (current<count) return get();
       return NULL;
    }

    long getcount() { return count; }


    long pos() { return current; }

    friend class NetRetGenerator;
};


class NetVecGenerator : public NetGen
{
  protected:
    VecNetwork &netvec;
    Network *get() {       
      return new Network(*(netvec[current])); }
    
  public:
    NetVecGenerator(VecNetwork &netvec) : netvec(netvec), NetGen(netvec.size()) 
    { reset(); }    
    
};


class QuasiConsTreeGenerator : public NetGen
{
  protected:
    
    Clusters *genetreeclusters;
    RootedTree *preserverootst;
    Clusters *guideclusters;
    Clusters *guidetree;
    float genetreessimilarity;

    Network *get() { 
      string r = genetreeclusters->genrootedquasiconsensus(preserverootst, guideclusters, guidetree, genetreessimilarity);
      if (!r.length())
      {
        cerr << "Cannot create initial species tree" << endl;
        exit(-1);
      }   
      return new Network(r); 
    }

  public:

    QuasiConsTreeGenerator(long count, RootedTree *preserverootst, Clusters *guideclusters, Clusters *guidetree, Clusters *genetreeclusters, float genetreessimilarity) :
       NetGen(count),
      genetreeclusters(genetreeclusters),
      preserverootst(preserverootst),
      guideclusters(guideclusters),
      genetreessimilarity(genetreessimilarity),
      guidetree(guidetree) {}

};

class RandTreeGenerator : public QuasiConsTreeGenerator
{
    
  public:

    RandTreeGenerator(long count, RootedTree *preserverootst, Clusters *guideclusters, Clusters *guidetree) :  QuasiConsTreeGenerator(count, preserverootst, guideclusters, guidetree, new Clusters(), 1.0) {}

    ~RandTreeGenerator() { delete genetreeclusters; }

};



class NetRetGenerator : public NetGen
{
 protected:
   //return addrandreticulations(reticulations, new Network(r), networkclass, timeconsistency, false, guideclusters, guidetree);

    NetGen *basegen;
    int reticulations;
    int networkclass;
    int timeconsistency;
    bool uniform;
    Clusters *guideclusters;
    Clusters *guidetree;    
    int samplespernet;

    Network *last;

    Network *get() 
    { 
      if (!last || (current%samplespernet)==0)
        last = basegen->next();      
      if (!last) return NULL;

      Network *k = addrandreticulations(reticulations, last, networkclass, timeconsistency, uniform, guideclusters, guidetree, true);            
      return k;
    }

    void reset()
    { 
      last = NULL;
      NetGen::reset();
      basegen->reset(); 
    }

 public: 
   NetRetGenerator(NetGen *basegen, int reticulations, int networkclass, int timeconsistency, bool uniform, long count = -1, int samplespernet=1, Clusters *guideclusters=NULL, Clusters *guidetree=NULL) :
      NetGen(count),
      basegen(basegen),
      reticulations(reticulations),
      networkclass(networkclass),
      timeconsistency(timeconsistency),
      uniform(uniform),
      samplespernet(samplespernet),
      guideclusters(guideclusters),
      guidetree(guidetree),
      last(NULL)
   { reset(); }       
};

class NetGenCollect : public NetGen
{

protected:
  vector<NetGen*> v;
  int currentgen;

  Network *get() 
  { 
    while (currentgen<v.size())
    {                
        Network *n = v[currentgen]->next();
        if (n) return n;

        currentgen++;    
        if (currentgen<v.size())
          v[currentgen]->reset();
        else break;
    }
    return NULL;    
  }


public:
  NetGenCollect(int count=-1) : NetGen(count) {}

  void reset()
  {    
    NetGen::reset();
    currentgen = 0;
    for (auto i: v)
    {
      i->reset();
    }
  }

  void add(NetGen *netgen)
  {
     v.push_back(netgen);
     netgen->reset();
  }


};


// NetGenerator 
//   - returns all networks from netvec
//   - returns quasiconsensus nets (-1 infinity), i.e., quasi cons trees +  reticulations random reticulations
//   - returns randomnetworks (-1 infinity) nets with reticulations reticulations
//   - returns all networks from next generator
//   - TODO: replace it 
class NetGenerator : public NetGen
{  
  int randomnetworks; 
  int quasiconsensusnetworks;
  VecNetwork *netvec;
  int reticulations;  
  int timeconsistency;
  int randnetuniform;
  RootedTree *preserverootst;
  Clusters *gtc;
  Clusters *guideclusters;
  Clusters *guidetree;
  int networkclass;
  NetGenerator *nextgenerator;
  int initial_randomnetworks;
  int initial_quasiconsensusnetworks;

public:

  NetGenerator(VecNetwork *netvec) : NetGenerator(netvec,0,0,NULL,NULL,0,0,0,0,NULL,NULL) {}

  NetGenerator(
    VecNetwork *netvec, 
    int randomnetworks, 
    int quasiconsensusnetworks, 
    Clusters *gtc,
    RootedTree *preserverootst,
    int reticulations, 
    int networkclass, 
    int timeconsistency, 
    bool randnetuniform,
    Clusters *guideclusters,
    Clusters *guidetree,
    NetGenerator *nextgenerator=NULL)
  : 
    NetGen(),
    netvec(netvec),
    initial_randomnetworks(randomnetworks),
    initial_quasiconsensusnetworks(quasiconsensusnetworks),
    reticulations(reticulations),    
    timeconsistency(timeconsistency),
    randnetuniform(randnetuniform),
    preserverootst(preserverootst),
    guideclusters(guideclusters),
    guidetree(guidetree),
    networkclass(networkclass),
    nextgenerator(nextgenerator),
    gtc(gtc) 

  {        
    reset();
  }

  Network *next();
  

  void reset() { 
    NetGen::reset(); 

    if (nextgenerator) nextgenerator->reset();

    quasiconsensusnetworks = initial_quasiconsensusnetworks;
    randomnetworks = initial_randomnetworks;
  }

};

Network *randquasiconsnetwork(int reticulations, int networktype, int timeconsistency, Clusters *genetreeclusters, RootedTree *preserverootst, 
  Clusters *guideclusters,
  Clusters *guidetree);


Network *randnetwork(int reticulations, int networktype, int timeconsistency, bool uniform);




#endif