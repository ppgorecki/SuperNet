
#include "testers.h"
#include "network.h"
#include "contrnet.h"
#include "treespace.h"
#include "randnets.h"
#include "neteditop.h"

extern TreeSpace *globaltreespace;

void testcomparedags(int reticulationcnt_R, int timeconsistency, int randnetuniform, int networkclass)
{

    int cnt = 0;
    while (1) 
    {
      cnt++;
      string r1 = randspeciestreestr();
      string r2 = randspeciestreestr();
      int EE = 0;
      int e1, e2, e3;
      if (!r1.length() || !r2.length()) {
        cerr << "Cannot create initial random species tree" << endl;
        exit(-1);
      }
      Network *n1 =
          addrandreticulations(reticulationcnt_R, new Network(r1), networkclass,
                               timeconsistency, randnetuniform, NULL, NULL);
      Network *n2 =
          addrandreticulations(reticulationcnt_R, new Network(r2), networkclass,
                               timeconsistency, randnetuniform, NULL, NULL);

      e1 = n1->eqdags(n2);
      e2 = n1->eqdagsbypermutations(n2);

      cout << (*n1) << "\t" << (*n2) << "\tE1=" << e1 << "\tE2=" << e2
           << "\tError=" << (e1 != e2) << endl;

      cerr << (*n1) << "\t" << (*n2) << "\tE1=" << e1 << "\tE2=" << e2
           << "\tError=" << (e1 != e2) << endl;

      delete n1;
      delete n2;

      if (e1 != e2)
        break; // stop
    }
    exit(0);
}

void testtreerepr(VecNetwork &netvec)
{
  for (auto &ntpos: netvec) 
    {
      DISPLAYTREEID tid = 0;
      SNode *t = NULL;
      cout << *ntpos << endl;
      while ((t = ntpos->gendisplaytree2(tid, t, globaltreespace)) != NULL) {
        ppSNode(cout, t) << endl;
        tid++;
      }
    }
}

void testcontract(int argc, char **argv, VecRootedTree &gtvec,  VecNetwork &netvec, CostFun *costfun)
{

    // Last two args: leftretusage rightretusage
    // ./supnet -g "(a,(b,c))" -n "((b)#1,(((#2)#3,((#1)#2,(c,#3))),a))" -eXng 2
    // 1
    // ./dot -Tpdf contr.dot -o c.pdf && evince c.pdf

    RETUSAGE retusage;
    emptyretusage(retusage);

    // ugly
    if (optind + 1 <= argc) {
      long us = atol(argv[optind]);

      int rid = 0;
      while (us) {
        cout << "A" << us << " " << rid << endl;
        if (us & 1)
          addleftretusage(retusage, rid);
        us = us >> 1;
        rid++;
      }
    }

    if (optind + 2 <= argc) {
      long us = atol(argv[optind + 1]);
      int rid = 0;
      while (us) {
        cout << "B" << us << " " << rid << endl;
        if (us & 1)
          addrightretusage(retusage, rid);
        us = us >> 1;
        rid++;
      }
    }

    cout << "R" << retusage << endl;
    cout << "Conflicted " << conflicted(retusage) << endl;
    
    for (int i = 0; i < netvec.size(); i++) 
    {
      Network *n1 = netvec[i];
      std::stringstream ss;
      n1->print(ss);
      cout << ss.str() << endl;

      ContractedNetwork *c = new ContractedNetwork(ss.str());

      c->contract(retusage);

      cout << "Rtcount " << c->rtcount() << " "
           << "rt=" << n1->rtcount() << endl;
      ;

      for (auto &gtpos: gtvec)
      {
        cout << "retmindc:" << (c->approxmindce(*gtpos, *costfun)) << endl;
      }

      ofstream s("contr.dot");
      std::ofstream sf;
      sf.open("contr.dot", std::ofstream::out);
      sf << "digraph SN {" << endl;
      sf << " inp [label=\"InRT=" << retusage << "\"]" << endl;
      c->gendot(sf);
      c->gendotcontracted(sf);
      sf << "}" << endl;
      sf.close();
      cout << c->newickrepr() << endl;

      RootedTree *t = NULL;
      DISPLAYTREEID tid = 0; // id of display tree

      while ((t = c->gendisplaytree(tid, t)) != NULL) {
        cout << tid << " " << (*t) << endl;
        tid++;
      }
    }
}

void testeditnni(VecNetwork &netvec)
{

   for (auto &ntpos: netvec) 
    {
      std::ofstream dotf;
      dotf.open("edit.dot", std::ofstream::out);

      int cnt = 0;
      dotf << "digraph SN {" << endl;

      NNI nni;
      nni.init(ntpos);

      cout << *ntpos << " verify==" << ntpos->verifychildparent() << endl;
      ntpos->printdeb(cout, 2);
      ntpos->printdot(dotf, cnt++);

      while (nni.next()) 
      {
        cout << "========" << endl << endl;
        int err = ntpos->verifychildparent();

        cout << *ntpos << " verify==" << err << endl;
        ntpos->printdeb(cout, 2);
        ntpos->printdot(dotf, cnt++);

        if (err)
        {
          break;
        }
      }

      dotf << "}" << endl;
      dotf.close();
    }
}

void testedittailmove(VecNetwork &netvec)
{
  for (auto &ntpos: netvec) 
{
      std::ofstream dotf;

      dotf.open("tm.dot", std::ofstream::out);

      int cnt = 0;
      dotf << "digraph SN {" << endl;

      TailMove tailmove(0);
      tailmove.init(ntpos);

      cout << *ntpos << " verify==" << ntpos->verifychildparent() << endl;
      // (*ntpos)->printdeb(cout,2);
      ntpos->printdot(dotf, cnt++);

      while (tailmove.next()) 
      {
        cout << "========" << endl << endl;
        int err = ntpos->verifychildparent();

        cout << ntpos << " verify==" << err << endl;
        // (*ntpos)->printdeb(cout,2);
        ntpos->printdot(dotf, cnt++);

        if (err)
        {
          break;
        }
      }

      dotf << "}" << endl;
      dotf.close();
    }
}



