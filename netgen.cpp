
#include "network.h"
#include "netgen.h"

// How many tries to obtain time consistent network
#define RANDCNTREPEAT 100 

Network *addrandreticulations(int reticulations, Network *n, int networkclass, int timeconsistency, bool uniform, Clusters *guideclusters, Clusters *guidetree, bool keepsource)
{   

  Network *src = n;

  for (auto k=0; k<RANDCNTREPEAT; k++)
  {

    if (!n->checktimeconsistency(timeconsistency, true))
    {      
        cerr << "Time consistency condition is not satisfied" << endl;
        exit(-1);                
    }

    for (auto i=0; i<reticulations; i++)
    {
      Network *prev = n;
      int addrandcnt = RANDCNTREPEAT;
      while (addrandcnt)
      {

        prev = n;

        n = n->addrandreticulation("", networkclass, uniform);
        if (!n)
        {
          cerr << "Cannot insert random " << (i+prev->rtcount()+1) << "-th reticulation into " << *prev << endl;
          exit(-1);
        }                

        if (guideclusters)
        {
          if (!n->hasclusters(guideclusters))
          {
              delete n; //revert
              n = prev;
              addrandcnt--;

              if (!addrandcnt)
              {
                cerr << "Cannot insert reticulation without violating guide clusters :(" << endl;
                exit(-1);
              }            
              continue;
          }
        }


        if (guidetree)
        {
          if (!n->hastreeclusters(guidetree))
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


        if (timeconsistency==NET_TIMECONSISTENT)
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

    if ((timeconsistency==NET_NOTIMECONSISTENT) && n->istimeconsistent())      
    {
        //revert  
        n = src;
    }
    else 
    {
      if (src != n && keepsource==false)
        delete src;

      return n;
    }
  }

  cerr << "Cannot generate network without time consistency" << endl;
  exit(-1);

}


Network *randnetwork(int reticulations, int networkclass, int timeconsistency, bool uniform)
{

    string r = randspeciestreestr();
    if (!r.length())
    {
      cerr << "Cannot create initial random species tree" << endl;
      exit(-1);
    }        
    Network *n = new Network(r);
    return addrandreticulations(reticulations, n, networkclass, timeconsistency, uniform, NULL, NULL);        
}

extern float opt_genetreessimilarity;

Network *randquasiconsnetwork(int reticulations, int networkclass, int timeconsistency, Clusters *genetreeclusters, RootedTree *preserverootst, 
  Clusters *guideclusters, Clusters *guidetree)
{
    string r = genetreeclusters->genrootedquasiconsensus(preserverootst, guideclusters, guidetree, opt_genetreessimilarity);
    if (!r.length())
    {
      cerr << "Cannot create initial quasi consensus species tree" << endl;
      exit(-1);
    }      

    return addrandreticulations(reticulations, new Network(r), networkclass, timeconsistency, false, guideclusters, guidetree);

}


// interator over networks
Network* NetGenerator::next()
{  
  current++;

  if (netvec && (netvec->size()>0 && current<netvec->size()))
      return new Network(*(*netvec)[current]);    

  if (nextgenerator)
  {
      Network *n  = nextgenerator->next();
      if (n) return n;
  }
  
  if (randomnetworks!=0)   // with -1 infitite 
  { 
    if (randomnetworks>0)  randomnetworks--;             

    return randnetwork(reticulations, networkclass, timeconsistency, randnetuniform);
  }
  
  if (quasiconsensusnetworks!=0)  // with -1 infitite 
  { 
      if (quasiconsensusnetworks>0) quasiconsensusnetworks--;              
      return randquasiconsnetwork(reticulations, networkclass, timeconsistency, gtc, preserverootst, guideclusters, guidetree);                     
  }

  return NULL;

}