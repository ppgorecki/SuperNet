
#include "tools.h"
#include "clusters.h"
#include "rtree.h"
#include "costs.h"



NODEID** RootedTree::getspclusterrepr()
{

    NODEID **res=new NODEID*[nn];

    for (NODEID i=0; i<lf; i++) 
    	res[i]=spec2gtcluster[lab[i]]->spcluster;

    for (NODEID i=lf; i<nn; i++) 
    	res[i]=joinspclusters(res[leftchild[i]],res[rightchild[i]]);

    return res;
  	
}

bitcluster* RootedTree::getbitclusterrepr()
{

    bitcluster *res=new bitcluster[nn];

    for (NODEID i=0; i<lf; i++) 
        res[i]=bcsingleton[lab[i]];

    for (NODEID i=lf; i<nn; i++) 
        res[i]=UNION(res[leftchild[i]],res[rightchild[i]]);

    return res;
    
}


// // Set levels of nodes, i.e., the height of the subtree rooted at i
// void RootedTree::setlevel(NODEID i)
// {
// 	if (i<lf) level[i] = 0;
// 	else
// 	{
//   		setlevel(leftchild[i]);
//   		setlevel(rightchild[i]);
//   		level[i] = 1 + max(level[leftchild[i]], level[rightchild[i]]);
// 	}
// }


string randspeciestreestr()
{

  if (!specnames.size())
  {
        cerr << "No species is defined (e.g., use -A NUM opt.)" << endl;
        exit(-1);
  }
  NODEID tr[specnames.size()];
  for (int i = 0; i < (int)specnames.size(); i++) tr[i] = i;
  return genrandomtree(tr, specnames.size());
}


RootedTree *randspeciestree()
{  
  return new RootedTree(randspeciestreestr());
}


double RootedTree::cost(RootedTree &speciestree, CostFun &cost)
{
    speciestree.initlca();
    speciestree.initdepth();    

    NODEID *lcamap = getlcamapping(speciestree);

    double v = cost.compute(*this, speciestree, lcamap);

    if (lcamap)
        delete[] lcamap;

    return v;

}

// double RootedTree::_cost(RootedTree &speciestree, NODEID* lcamap, CostFun* costfunc)
// {
//     double ccost = 0;

//     switch (costfunc)
//     {
//         case COSTDUPLICATION: 
//             ccost = costduplication(speciestree,lcamap);  
//             break;
//         case COSTLOSS: 
//             ccost = costloss(speciestree,lcamap); 
//             break;
//         case COSTDUPLICATIONLOSS: 
//             ccost = costduplicationloss(speciestree,lcamap); 
//             break;
//         case COSTDEEPCOAL: 
//             ccost = costdeepcoalx(speciestree,lcamap); 
//             break;
//         case COSTROBINSONFOULDS: 
//             ccost = costrobinsonfoulds(speciestree,lcamap); 
//             break;
//         default: 
//             cerr << "Unknown costfunc type?" << endl;
//             exit(-1);
//     }

//     return ccost;    
// }


// long RootedTree::costloss(RootedTree &speciestree, NODEID *lcamap)
// {
//     return costdeepcoalx(speciestree,lcamap)+2*costduplication(speciestree,lcamap);
// }

// long RootedTree::costrobinsonfoulds(RootedTree &speciestree, NODEID *lcamap)
// {
//     //TODO
//     cerr << "RF is not implemented yet" << endl;
//     exit(-1);    
// }

// long RootedTree::costduplicationloss(RootedTree &speciestree, NODEID *lcamap)
// {
//   return costdeepcoalx(speciestree,lcamap)+3*costduplication(speciestree,lcamap);
// }


// long RootedTree::costduplication(RootedTree &speciestree, NODEID *lcamap)
// {      
//   long s = 0;
//   for (NODEID i=lf; i<nn; i++)        
//     if (lcamap[leftchild[i]]==lcamap[i] || lcamap[rightchild[i]]==lcamap[i]) 
//             s++; 
//   return s;
// }

// long RootedTree::costdeepcoalx(RootedTree &speciestree, NODEID *lcamap)
// {      
//   long s = 0;            
//   for (NODEID i=0; i<nn; i++) 
//         if (i!=root) 
//             s+=speciestree.depth[lcamap[i]]-speciestree.depth[lcamap[parent[i]]]-1; 
//   return s;
// }

void RootedTree::inittreestr()
{
    depth = new NODEID[nn];
#ifdef LCATAB
    lcatab = NULL;
#endif
}

// Build caterpillar tree
RootedTree::RootedTree(int _lf, NODEID* labels, double weight): Dag(_lf,labels,weight), depthinitialized(false)
{   
    inittreestr();
}

RootedTree::RootedTree(char *s, double weight): Dag(s,weight), depthinitialized(false) 
{
    if (rt>0)
    {
        cerr << "A tree cannot have reticulations!" << endl;
        exit(-1);
    }       
    inittreestr();      
}


// Set depth of nodes in a subtree 
void RootedTree::_setdepth(NODEID i, int dpt)
{
    depth[i] = dpt;
    if (i<lf) return;
    _setdepth(leftchild[i], dpt + 1);
    _setdepth(rightchild[i], dpt + 1);
}

void RootedTree::initdepth() 
{ 
    if (!depthinitialized)
    {
        _setdepth(root, 0);       
        depthinitialized = true;
    }
}

NODEID* RootedTree::getlcamapping(RootedTree &speciestree, NODEID* lcamap)
{
    if (nn<1) return NULL;
    if (!lcamap) lcamap = new NODEID[nn];
    for (NODEID i=0; i<lf; i++) 
        lcamap[i] = speciestree.findlab(lab[i]);    
    for (NODEID i=lf; i<nn; i++)
        lcamap[i] = speciestree.lca(lcamap[leftchild[i]], lcamap[rightchild[i]]);
    return lcamap;
}

void RootedTree::initlca()
{
   
    initdepth();

#ifdef LCATAB
    if (!lcatabinitialized)
    {
        if (!lcatab) 
        {
            lcatab = new NODEID*[nn];
            for (NODEID i = 1; i < nn; i++) lcatab[i] = new NODEID[i];
        }
        for (NODEID i = 1; i < nn; i++)
            for (NODEID j = 0; j < i; j++)
                lcatab[i][j] = _lca(i, j);

        lcatabinitialized = true;
    }
#endif
  }


#ifdef LCATAB

NODEID RootedTree::lca(NODEID a, NODEID b)
{
    //cout << a << " " << b <<  " " << depth[a] << " " << depth[b] << endl;
    if (a == b) return a;
    if (!lcatabinitialized) initlca();
    if (a < b) return lcatab[b][a];
    return lcatab[a][b];
}

#endif


#ifdef LCATAB
NODEID RootedTree::_lca(NODEID a, NODEID b) 
#else
NODEID RootedTree::lca(NODEID a, NODEID b) 
#endif
{
    //cout << a << " " << b <<  " " << depth[a] << " " << depth[b] << endl;
    while (depth[a]<depth[b]) b = parent[b];
    while (depth[b]<depth[a]) a = parent[a];

    // equal dephts; iterate until lca is reached
    while (a!=b)
    {            
      a = parent[a];
      b = parent[b];      
    }
    return a;
}



// NODEID* RootedTree::repr2(NODEID *t)
// {
//     if (!t) t = new NODEID[lf*3-2];
//     NODEID mm[nn];
//     _setmm_repr(root, mm);
//     int cpos = 0;
//     _repr(t, root, cpos, mm);
//     return t;
// }

// Unique representation 
// Generate unique representation of a tree as a vector
// (sorted representation)
// If t==NULL new array will be allocated
SPID* RootedTree::repr(SPID *t)
{
    if (!t) t = new SPID[lf*3];
    SPID mm[nn];
    t[0]=lf*3-1;
    t[1]=_setmm_repr(root, mm);
    int cpos = 0;
    _repr(t+2, root, cpos, mm);
    return t;
}

// (sorted representation)
// Used in sorted representation of tree to compute 
// a node in a subtree with min spid
SPID RootedTree::_setmm_repr(NODEID i, SPID* mm)
{
    if (i<lf) mm[i] = lab[i];
    else
    {
      NODEID lm = _setmm_repr(leftchild[i], mm);
      NODEID rm = _setmm_repr(rightchild[i], mm);
      mm[i] = min(lm, rm);
    }
    return mm[i];
}

// Used in sorted representation of tree
void RootedTree::_repr(SPID *t, NODEID i, int &pos, SPID* mm)
{
    if (i<lf)
    {
      t[pos++] = lab[i];
      return;
    }
    t[pos++] = REPROPEN; // (

    if (mm[leftchild[i]] == mm[i])
    {
      _repr(t, leftchild[i], pos, mm);
      // t[pos - 1] = -t[pos - 1]; //,
      _repr(t, rightchild[i], pos, mm);
    }
    else
    {
      _repr(t, rightchild[i], pos, mm);
      // t[pos - 1] = -t[pos - 1]; //,
      _repr(t, leftchild[i], pos, mm);
    }
    t[pos++] = REPRCLOSE; // )
  }

ostream& RootedTree::printdeb(ostream&s, int gse, string tn)
{
    Dag::printdeb(s,gse,tn);
    if (depthinitialized)
    {
        s << "depth=[";
        for (NODEID i=0; i<nn; i++)
            s << (int)depth[i] << " ";
        return s << "]";
    }
    else
    {
      s << " depth=uninitialized";
    }
    return s << endl;
}