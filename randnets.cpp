
#include "network.h"

// How many tries to obrain time consistent network
#define RANDCNTREPEAT 100 

Network *addrandreticulations(int reticulationcnt, Network *n, int networkclass, int timeconsistency, bool uniform)
{   

  Network *src = n;

  for (auto k=0; k<RANDCNTREPEAT; k++)
  {

    if (timeconsistency)
    {
      bool tc = n->istimeconsistent();

      if (timeconsistency==TIMECONSISTENT && !tc)
      {
        cerr << "Network " << *n << " is not time consistent" << endl;
        exit(-1);
      }            
    }

    for (auto i=0; i<reticulationcnt; i++)
    {
      Network *prev = n;
      int addrandcnt = RANDCNTREPEAT;
      while (addrandcnt)
      {

        prev = n;

        n = n->addrandreticulation("", networkclass, uniform);
        if (!n)
        {
          cerr << "Cannot insert random " << (i+1) << "-th reticulation into " << *prev << endl;
          exit(-1);
        }                

        if (!timeconsistency) 
          // OK found
          break;


        if (timeconsistency==TIMECONSISTENT)
        {
            bool tc = n->istimeconsistent();
            if (tc)
            {
                // OK net is found            
                break;
            }
            else
            {         
                //TC condition not satisfied: revert
                delete n;
                n = prev;
            }        

            addrandcnt--;

            if (!addrandcnt)
            {
              cerr << "Cannot insert reticulation without violating time consistency condition" << endl;
              exit(-1);
            }            
        }       
        else break;

      } // while randcnt

      if (prev!=src)
        delete prev;      
    } // for 

    if ((timeconsistency==NOTIMECONSISTENT) && n->istimeconsistent())      
    {
        //revert  
        n = src;
    }
    else 
    {
      if (src != n)
        delete src;

      return n;
    }
  }

  cerr << "Cannot generate network without time consistency" << endl;
  exit(-1);

}


Network *randnetwork(int reticulationcnt, int networkclass, int timeconsistency, bool uniform)
{

    string r = randspeciestreestr();
    if (!r.length())
    {
      cerr << "Cannot create initial random species tree" << endl;
      exit(-1);
    }        
    return addrandreticulations(reticulationcnt, new Network(r), networkclass, timeconsistency, uniform);        
}

Network *randquasiconsnetwork(int reticulationcnt, int networkclass, int timeconsistency, TreeClusters *gtc, RootedTree *preserverootst)
{
    string r = gtc->genrootedquasiconsensus(preserverootst);
    if (!r.length())
    {
      cerr << "Cannot create initial quasi consensus species tree" << endl;
      exit(-1);
    }      

    return addrandreticulations(reticulationcnt, new Network(r), networkclass, timeconsistency, false);

}


// interator over networks
// i should be initialized with -1
Network* netiterator(long int &i, VecNetwork &netvec, int &randomnetworkscnt, int &quasiconsensuscnt, 
  TreeClusters *gtc,
  RootedTree *preserverootst,
  int reticulationcnt, int networkclass, int timeconsistency, bool randnetuniform)
{
  
  i++;

  if (netvec.size()>0 && i<netvec.size())  
      return netvec[i];    
    
  if (randomnetworkscnt!=0)   // with -1 infitite 
  { 
    if (randomnetworkscnt>0)  randomnetworkscnt--;             

    return randnetwork(reticulationcnt, networkclass, timeconsistency, randnetuniform);

  }


  if (quasiconsensuscnt!=0)  // with -1 infitite 
  { 
      if (quasiconsensuscnt>0) quasiconsensuscnt--;              
      return randquasiconsnetwork(reticulationcnt, networkclass, timeconsistency, gtc, preserverootst);                     
  }

  return NULL;

}