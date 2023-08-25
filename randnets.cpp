
#include "network.h"

// How many tries to obrain time consistent network
#define RANDCNTREPEAT 100 

Network *addrandreticulations(int reticulationcnt, Network *n, int networkclass, int timeconsistency, bool uniform, Clusters *guidetreeclusters)
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

        if (guidetreeclusters)
        {
          if (!n->hasclusters(guidetreeclusters))
          {
              delete n; //revert
              n = prev;
              addrandcnt--;

              if (!addrandcnt)
              {
                cerr << "Cannot insert reticulation without violating guide tree :(" << endl;
                exit(-1);
              }            
              continue;
          }

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

            //TC condition not satisfied: revert
            delete n;
            n = prev;
                   
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
      {
        delete prev;      
      }

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
    return addrandreticulations(reticulationcnt, new Network(r), networkclass, timeconsistency, uniform, NULL);        
}

Network *randquasiconsnetwork(int reticulationcnt, int networkclass, int timeconsistency, Clusters *genetreeclusters, RootedTree *preserverootst, Clusters *guidetreeclusters)
{
    string r = genetreeclusters->genrootedquasiconsensus(preserverootst, guidetreeclusters);
    if (!r.length())
    {
      cerr << "Cannot create initial quasi consensus species tree" << endl;
      exit(-1);
    }      

    return addrandreticulations(reticulationcnt, new Network(r), networkclass, timeconsistency, false, guidetreeclusters);

}


// interator over networks
// i should be initialized with -1
Network* netiterator(long int &i, VecNetwork &netvec, int &randomnetworkscnt, int &quasiconsensuscnt, 
  Clusters *gtc,
  RootedTree *preserverootst,
  int reticulationcnt, int networkclass, int timeconsistency, 
  bool randnetuniform,
  Clusters *guidetreeclusters)
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
      return randquasiconsnetwork(reticulationcnt, networkclass, timeconsistency, gtc, preserverootst, guidetreeclusters);                     
  }

  return NULL;

}