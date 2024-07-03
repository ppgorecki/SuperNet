#include <iostream>
#include <algorithm>
#include <regex>
#include <string>
#include <cctype>
#include <sstream> 
#include <vector> 

#include "algotok.h"
#include "netgen.h"
#include "hillclimb.h"

#define TLABEL 1
#define TNUM 2
#define TRNUM 3

using namespace std;

const std::string AlgTokenizer::DELIMITERS(" \t\n\r;");

string types2names[] = { "NULL", "Label", "Int", "Float", "String", 
                    "Networks", "List", "Vector of Values", "MAIN", "Assignment",
                    "repeat", "", "", "", "", 
                    "",
                    "Call", "hcalgorithm", "italgorithm", "storebest", "print", "printenv", "read", "netgenerator"
                };



ExprList* parseExprList(AlgTokenizer &s, char sep, char termtoken=')');


#define ELAST ENETS

string ptokenname(int tokentype)
{
    if (tokentype>=0 && tokentype <= ELAST)
        return types2names[tokentype];
    return string(1,(char)tokentype);
}

AlgTokenizer::AlgTokenizer(const std::string& s) :
    m_prev(""),
    m_string(s), 
    m_offset(0), 
    m_delimiters(DELIMITERS),
    linenum(0) {}

AlgTokenizer::AlgTokenizer(const std::string& s, const std::string& delimiters) :
    m_prev(""),
    m_string(s), 
    m_offset(0), 
    m_delimiters(delimiters),
    linenum(0) {}

int AlgTokenizer::expect(const std::string s)
{
    if (m_token==s) return true;
    cerr << "Line " << linenum << ". Expected " << s << " found " << m_token << " after " <<  m_prev << endl;
    exit(-1);
}    

int AlgTokenizer::next() 
{
    return next(m_delimiters);
}

bool islabel(char c)
{
    return isalpha(c) || c=='_';
}

bool isnum(char c)
{
    // ugly
    return isdigit(c) || c=='.' || c=='-' || c=='+';
}

int AlgTokenizer::_next(const std::string& delimiters) 
{    

    m_prev = m_token;

    size_t i = m_string.find_first_not_of(delimiters, m_offset);    

    if (std::string::npos == i) 
    {
        m_offset = m_string.length();
        return 0;
    }

    for (int p = m_offset; p<i; p++)
        if (m_string[p]=='\n') linenum++;

    if (m_string[i]=='#')
    {
        size_t i = m_string.find('\n', m_offset);
        if (std::string::npos == i) 
        {
            m_offset = m_string.length();
            return 0;
        }
        m_offset = i+1;        
        linenum++;
        return next(delimiters);
    }


    if (strchr("=(){},",m_string[i]))
    {
        m_token = m_string.substr(i,1);        
        m_offset = i+1;        
        return  m_string[i];
    }

    if (m_string[i]=='"')
    {        

        if (m_string.find("\"\"\"", i)==i)
        {
            // Multiline comment
            i+=3;
            int j;
            if ((j = m_string.find("\"\"\"", i+3))==std::string::npos)
            {
                cerr << "Unterminated \"\"\" comment" << endl;
                exit(-1);
            }
            for (;i<j;i++)
                if (m_string[i]=='\n') linenum++;
            m_offset = j+1;                                
            return next(delimiters);

        }

        size_t j = m_string.find_first_of("\"", i+1);  // TODO: escape \"
        if (std::string::npos == j)    
        {
            cerr << "Line " << linenum << ". Unterminated string: " << m_string.substr(i,10) << endl;
            exit(-1); 
        }
        m_token = m_string.substr(i+1, j-i-1);
        m_offset = j+1;

        return ESTR;
    }

    if (islabel(m_string[i])) 
    {
        int j=i;
        while (islabel(m_string[j])) j++;        
        m_token = m_string.substr(i, j-i);        
        m_offset = j;        
                                    
        if (m_token=="randtrees")
        {
            return ERANDTREES;
        } 
        else if (m_token=="netlist")
        {
            return ENETLIST;
        }  
        else if (m_token=="quasiconsensustrees")
        {
            return EQUASICONSTREES;
        }                           
        else if (m_token=="addreticulations")
        {
            return EADDRETICULATIONS;
        }                                   
        else if (m_token=="climb")
        {
            return ECLIMB;
        }        
        else if (m_token=="storebest")
        {
            return ESTOREBEST;
        }
        else if (m_token=="print")
        {
            return EPRINT;
        }
        if (m_token=="println")
        {
            return EPRINTLN;
        }
        if (m_token=="printsp")
        {
            return EPRINTSP;
        }
        if (m_token=="showstats")
        {
            return ESHOWSTATS;
        }
        else if (m_token=="for")
        {
            return EFOR;
        }
        else if (m_token=="to")
        {
            return ETO;
        }
        else if (m_token=="printenv")
        {
            return EPRINTENV;
        }
        else if (m_token=="read")
        {
            return EREAD;
        }
        else 
        { 
            return ELABEL;
        }
    }

    if (isnum(m_string[i])) 
    {
        int j=i;
        int d=0;
        int dot=0;

        if (m_string[j]=='-' || m_string[j]=='+') j++;    
        while (isdigit(m_string[j])) { j++; d++; } 
        if (m_string[j]=='.') { 
            dot=true;
            j++; 
            while (isdigit(m_string[j])) {j++; d++;}
        }

        if (!d)
        {
            cerr << "Line " << linenum << ". Unknown token: " << m_string.substr(i,10) << endl;
            exit(-1);
        }

        m_token = m_string.substr(i, j-i);        
        m_offset = j;

        std::size_t pos; 

        if (dot)
        {
            floatval = stof(m_token, &pos);            
            return ERNUM;
        }

        intval = stoi(m_token, &pos);
        floatval = intval;
        return ENUM;
                    
    }

    cerr << "Line " << linenum << ". Unknown token at " << m_string.substr(i,10) << "..." << endl;
    exit(-1);
    
    return 0;
}



Expr* parseValue(AlgTokenizer &s)
{   
    int tokentype = s.gettokentype();
#ifdef TOKPARSEDEBUG        
    cout << "parseValue: " << tokentype << " " << s.get() << endl;
#endif    
    Expr *res = NULL;
    if (tokentype==ERNUM)
    {
        res = new VFloat(s.floatval);
    }
    else if (tokentype==ENUM)
    {
        res =  new VInt(s.intval);
    }
    else if (tokentype==ESTR)
    {
        res =  new VStr(s.get());
    }
    else if (tokentype==ELABEL)
    {
        res =  new ELabel(s.get()); // TODO: arithmetic expressions?
    }
    else if (tokentype & ECALL)
    {        
        s.next(); // eat label
        s.expect("("); 
        s.next(); // eat (
        res = new ExprCall(tokentype, parseExprList(s,','));        
        s.expect(")");    
    }    
    if (res) s.next(); // eat
    else
    {
        cerr << "Expected value. Found " << ptokenname(tokentype) << endl;
        exit(-1);
    }
    return res;
}


void ExprList::append(Expr *e)
{        
    t.push_back(e); 
}
//#define TOKPARSEDEBUG

int cx = 0;
ExprList* parseExprList(AlgTokenizer &s, char sep, char termtoken)
{      
      ExprList *el = new ExprList();      
      int tokentype = s.gettokentype();
      int command = 0;
#ifdef TOKPARSEDEBUG
      cout << "parseExprList:" << s.get() << " sep=" << sep << " termtoken=" << termtoken << endl;
#endif

      while (tokentype != 0)
      {             
          
          if (tokentype == termtoken)
          {
                break;
          }
          else if (tokentype == ELABEL)
          {
                string label = s.get();                
                tokentype = s.next();                
                
                if (tokentype=='=')                    
                {                                                 
                    s.next();
                    el->append(new ExprAssignment(label, parseValue(s)));             
                }    
                else if (tokentype == sep || tokentype ==')' || tokentype==0)
                {                                        
                    el->append(new ELabel(label));                                
                }                
                else
                {
                    cerr << "Line " << s.linenum << ". Expected label, value or assignment. Found: " << s.get() << " after " << s.prev() << endl;  
                    exit(-1);
                }
                tokentype = s.gettokentype();

                if (sep && tokentype == sep) 
                    s.next();
          }          
          else if (tokentype == EFOR)
          {              
              ExprAssignment *ea = NULL;
              s.next(); // eat for

              string label = s.get();                
              tokentype = s.next();                
                
              if (tokentype=='=')                    
              {                                                 
                  s.next();
                  ea = new ExprAssignment(label, parseValue(s));             
              }   
              else
              {
                  cerr << "Line " << s.linenum << ". Expected label, value or assignment. Found: " << s.get() << " after " << s.prev() << endl;  
              } 
              
              if ((tokentype = s.gettokentype())!=ETO)
              {
                 cerr << "Line " << s.linenum << ". Expected to in for loop. Found: " << s.get() << " after " << s.prev() << endl;  
                  exit(-1); 
              }
              s.next(); // eat to
              Expr *et = parseValue(s);                            

              if ((tokentype = s.gettokentype())!='{')
              {
                  cerr << "Line " << s.linenum << ". Expected block {...} in for. Found: " << s.get() << " after " << s.prev() << endl;  
                    exit(-1);
              }

              s.next(); // eat {

              ExprList *block = parseExprList(s,0,'}');
              tokentype = s.gettokentype();

              if (tokentype != '}')
              {
                  cerr << "Line " << s.linenum << ". Expected block } in repeat. Found: " << s.get() << " after " << s.prev() << endl;  
                    exit(-1);
              }
              s.next(); // eat {

              el->append(new ExprFor(ea, et, block));   
          } 
          else
          {                        
            
              el->append(parseValue(s)); // take directly call            
              tokentype = s.gettokentype();
              if (sep && tokentype == sep) 
                    s.next();              
          }

          tokentype = s.gettokentype();
          
      }
#ifdef TOKPARSEDEBUG
      cout << "[Done] parseExprList:" << s.get() << " sep=" << sep << " termtoken=" << termtoken << endl;
#endif      
      return el;
}


ExprList* parsescript(string script)
{
      AlgTokenizer s(script);
      s.next();
      return parseExprList(s, 0, 0);
}


ExprValue* ExprFor::eval(Env &env)
{
    string label = init->label;
    ExprValue *prev = env.get(label);
    
    ExprValue *var = init->eval(env);    

    int start = env.getint(string("for loop"), var);
    int term = env.getint(string("for loop"), to->eval(env));

    VInt it(start);
    env.setvar(label, &it);

    if (start<=term) 
    {
    
        int step=1;    
        int i = start - step;
        do 
        {
            i+=step;
            it.v = i;
            block->eval(env);        

        } while (i!=term);
    }

    if (prev)
    {        
        env.setvar(label, prev); // revert
    }
    else
    {
        env.delvar(label);
    }
    
    return var;
}

VNetGen* applyvec(NetGen &netgen, bool cleanunused, bool vectorize)
{
    if (vectorize)
    {
        VecNetwork *nets = new VecNetwork();
        Network *n; 
        while ((n = netgen.next())!=NULL)      
            nets->push_back(n);
        return new VNetGen(nets, cleanunused);

    }
    return new VNetGen(netgen, cleanunused);
    
}

ExprValue *ExprCall::eval(Env &env)
{   
    Env localenv(env);  
    localenv.setparentenv(env);

    EValueVec *a = (EValueVec*)args->eval(localenv);
    VecNetwork *nets = new VecNetwork();

    auto preserverootst = localenv.preserverootst; // Todo: via variable
    auto guideclusters = localenv.guideclusters; // Todo: via variable
    auto guidetree = localenv.guidetree; // Todo: via variable

    int stopatcostdefined = localenv.has("stopatcost");
    int stopatcost = localenv.getint("stopatcost");
    int vectorize = localenv.getint("vectorize");

    if (type==ERANDTREES)
    {        
        int count = localenv.getint("count", 1);
        bool cleanunused = localenv.getint("cleanunused", false);    
        auto netgenerator = new RandTreeGenerator(
                  count,                                    
                  preserverootst,                   
                  guideclusters,
                  guidetree);

        return applyvec(*netgenerator, cleanunused, vectorize);                    

    }


    if (type==ENETLIST)
    {        
        bool cleanunused = localenv.getint("cleanunused", false);        
        int count = localenv.getint("count", -1); 
        NetGenCollect *netcol = new NetGenCollect(count);
        string file = "";

        for (int i=0; i<a->size(); i++)
        {
            auto v = a->get(i);

            // Check if variable is save and extract file name
            if (args->t[i]->type==EASSIGNMENT && ((ExprAssignment*)args->t[i])->label=="save")
            {
                file = v->getstr("file");
                continue;
            }

            if (v->type==ENETS)
            {
                netcol->add( &((VNetGen*)v)->netgen );
            }            
            else if (v->type==EVALUEVEC)
            {
                VecNetwork *nets = new VecNetwork();
                EValueVec *evv = (EValueVec*)v;
                for (int i=0; i< evv->size(); i++)
                {                    
                    nets->push_back(new Network(evv->get(i)->getstr("netlist")));
                }
                netcol->add( new NetVecGenerator(*nets) );
            }
            else if (v->type==ESTR)
            {                
                VecNetwork *nets = new VecNetwork();
                stringstream ss(v->getstr("netlist")); 
                string strtree; 

                while (getline(ss, strtree, ';')) { 
                    nets->push_back(new Network(strtree));                        
                }
                netcol->add( new NetVecGenerator(*nets) );                                 
            }
        }
        auto res = applyvec(*netcol, cleanunused, vectorize);                     
        //new VNetGen(*netcol, cleanunused);  
        
        if (file!="")
        {
            res->save(file);
        }

        return res;                           
    }

    if (type==EQUASICONSTREES)
    {        
        int count = localenv.getint("count", 1);

        bool cleanunused = localenv.getint("cleanunused", false);
        
        float genetreessimilarity = localenv.getfloat("genetreessimilarity");

        auto genetreeclusters = new Clusters();
        for (auto & gtpos: localenv.genetreesv) 
            genetreeclusters->adddag(gtpos);
        
        auto netgenerator = new QuasiConsTreeGenerator(
                  count,                  
                  preserverootst,                   
                  guideclusters,
                  guidetree,
                  genetreeclusters,
                  genetreessimilarity);

        return applyvec(*netgenerator, cleanunused, vectorize);                  
    }

    if (type==EADDRETICULATIONS)
    {        
        int count = localenv.getint("count", -1); 
        int networkclass = localenv.getint("networkclass"); 
        int reticulations = localenv.getint("reticulations");
        
        VNetGen *source = localenv.getsource("addreticulations");
            
        int randnetuniform = localenv.getint("uniform");
        int timeconsistency = localenv.getint("timeconsistency");       
        bool cleanunused = localenv.getint("cleanunused", false);
        int samplespernet = localenv.getint("samplespernet", 1);
            
        auto netgenerator = new NetRetGenerator(
                  &source->netgen,
                  reticulations,                  
                  networkclass,
                  timeconsistency,
                  randnetuniform,
                  count,
                  samplespernet,

                  guideclusters,
                  guidetree);                
        return applyvec(*netgenerator, cleanunused, vectorize);     
        

    }

    
    
    if (type==EREAD)
    {           
        if (!a->size())
        {
            cerr << "At last one arg expected in read" << endl;               
        }
        vector<ExprValue*> evec;
        for (auto i=0; i<a->size(); i++)
        {
            ExprValue  *arg = a->get(i);

            if (arg->type!=ESTR)
            {
                cerr << "Expected string in read function. Found " << arg->type << endl;                               
            }
            std::ifstream f(((VStr*)arg)->s);
            std::vector<string> v;
            string line;
            while (f >> line) evec.push_back(new VStr(line));
        }        
        return new EValueVec(evec);
    }

    if (type==ESTOREBEST)
    {
        
        if (!localenv.genetreesv.size())
        {
            cerr << "No gene trees defined..." << endl;
            exit(-1);
        }

        VNetGen *source = localenv.getsource("storebest");
        
        int count = localenv.getint("count", 1);                
        float displaytreesampling = 0; // TODO: sampling
        int verbosehccost = localenv.getint("verbosehccost");
        int stopinit = localenv.getint("stopinit");
        

        DagSet visiteddags;       
        ClimbStats climbstats(visiteddags, localenv.getglobalstats());  
        climbstats.start();

        DagScoredQueue dgs(count);
        COSTT optcost;
        bool first = true;

        long int hccnt = -1;
        int lastimprovement = 0;
        int improvementscnt = 0;
        int noimprovementstep = 0;
        int cnt = 0;

        for (auto net: *source)
        {
            if (stopinit && (noimprovementstep > stopinit))
                break; // stop

            double curcost = net->odtcost(
                localenv.genetreesv, 
                *localenv.costfun, 
                localenv.getint("usenaive"), 
                localenv.getint("runnaiveleqrt"), 
                climbstats.getodtstats(), 
                displaytreesampling);

            // cout << noimprovementstep << " " << curcost << *net << endl;

            cnt++;
            hccnt++;

            Dag *n = dgs.push(net, curcost);

            if (first)
            {
                climbstats.addnewbest(*net, curcost);
                optcost = curcost;
                first = false;

            }
            else if (curcost<optcost)
            {                
                optcost = curcost;                      

                if (verbosehccost>=1)
                {
                    cout << "   >: " << *net << " cost=" << optcost << endl;    
                }

                // new optimal; forget old   
                climbstats.addnewbest(*net, optcost);                   
            
                improvementscnt++;
                noimprovementstep = 0; // reset counter

                if (stopatcostdefined && stopatcost >= curcost)
                {
                    break;
                }
                
                // reset 
                // lasthcimprovementtime = gettime();
                    
            } else
            {
                noimprovementstep++;
                if (source->cleanunused && n)
                {
                    delete n;
                }
            }



        }

        climbstats.finalize();

        while (!dgs.empty())
        {
            DagScored ds = dgs.top();
            // cout << ds.cost << " " << *ds.dag << endl;
            nets->push_back((Network*)ds.dag);
            dgs.pop();            
        }

        env.globalstats->merge(climbstats, env.getint("printstats"), false);

        // return applyvec(*netgenerator, cleanunused, vectorize);     
        return new VNetGen(nets);
    }

    if (type==ECLIMB)
    {
        
        int flag_hcedit_nni = localenv.getint("hceditnni");
        int networkclass = localenv.getint("networkclass"); 

        
        EditOp *editop;
        if (flag_hcedit_nni)
            editop = new NNI();
        else
            editop = new TailMove(networkclass, guideclusters, guidetree);

        DagScoredQueue dgs(localenv.getint("return",1));

        vector<ClimbStatsGlobal*> globalstatsarr;
         globalstatsarr.push_back(env.globalstats);         

        supnetheuristic(   
            localenv.genetreesv,       
            &localenv.getsource("hcalgorithm")->netgen,
            editop,
            localenv.costfun,
            env.getint("printstats"),       
            localenv.getint("stopinit"),
            localenv.getint("stopclimb"),            
            localenv.getint("usenaive"),            
            localenv.getint("runnaiveleqrt"),
            localenv.getint("maximprovements"),
            globalstatsarr,
            localenv.getint("cutwhendtimproved"),
            localenv.getint("multipleoptima", 0),
            stopatcostdefined, 
            stopatcost,
            &dgs
        );

        delete editop;


        if (localenv.getint("hcdetailedsummary"))
        {
           // print summary stats     
           for (auto nhc: globalstatsarr) 
           {        
              cout << "HC stats: "; 
              nhc->info(cout);
              cout << " ";        
              nhc->print();        
              cout << endl;
           }
        }

        globalstatsarr.pop_back(); 

        int verbosealg = localenv.getint("verbosealg");

        // merge all data and save odt/dat file(s); optional     
        localenv.globalstats->savedatmerged(verbosealg >= 4, globalstatsarr, true);     

        // best dags to file
        localenv.globalstats->savebestdags(verbosealg >= 4, true);

        // print summary
        // localenv.globalstats->print();
        // localenv.globalstats->printnetworkinfo();     
        // cout << endl;

        while (!dgs.empty())
        {
            DagScored ds = dgs.top();
            // cout << ds.cost << " " << *ds.dag << endl;
            nets->push_back((Network*)ds.dag);
            dgs.pop();            
        }


        return new VNetGen(nets);

    }

    if (type==EPRINT || type==EPRINTLN || type==EPRINTSP)
    {

        if (type==EPRINT)
            a->print(cout, ' ');
        else if (type==EPRINTLN)
            a->print(cout, '\n');
        else a->print(cout);
        
        if (type!=EPRINTLN) 
            cout << endl;

        return NULL;
    }

    if (type==ESHOWSTATS)
    {
        ClimbStatsGlobal *globalstats = env.getglobalstats();
        globalstats->print();
        cout  << endl;
        return NULL;
    }

    if (type==EPRINTENV)
    {
        cout << env << endl;
        return NULL;
    }

    cerr << "Warning: unknown call " << type << endl;
    return NULL;
}

ExprValue *ELabel::eval(Env &env)
{    
    if (env.vars.find(s)==env.vars.end())
    {
        if (s=="geneclusters")        
        {
            // Special variable -> return string
            auto genetreeclusters = new Clusters();
            for (auto & gtpos: env.genetreesv) 
                genetreeclusters->adddag(gtpos);


            stringstream ss; 
            ss << *genetreeclusters << endl;

            return new VStr(ss.str());

    

        }
        cerr << "Unknown variable <<" << s << ">>" << endl;
        exit(-1);
    }   
    return env.vars[s];        
}

ExprValue* ExprList::eval(Env &env)
{    
    vector<ExprValue*> res;
    for (auto p: t)
    {
        ExprValue *ev = p->eval(env);
        res.push_back(ev);

#ifdef TOKEVALDEBUG        
        cout << "evalExprList:" << p << " " << p->type << endl;
        if (ev)
             cout << "elinsert:" << *ev << endl;
#endif        
    }
    return new EValueVec(res);
}

extern TreeSpace *globaltreespace;
extern int maxdisplaytreecachesize;

ExprValue *ExprAssignment::eval(Env &env)
{    
    ExprValue *ev = value->eval(env);        

    if (label=="genetrees")
    {
        VecRootedTree *genetreesv = env.getgenetrees();

        int startid = genetreesv->size();

        if (ev->type==EVALUEVEC)
        {
            EValueVec *evv = (EValueVec*)ev;
            for (int i=0; i< evv->size(); i++)
            {
                RootedTree *gtree = new RootedTree(evv->get(i)->getstr("genetrees"));
                gtree->setid(startid++);
                genetreesv->push_back(gtree);
            }
        }
        else if (ev->type==ESTR)
        {
            
            stringstream ss(ev->getstr("genetrees")); 
            string strtree; 

            while (getline(ss, strtree, ';')) { 
                RootedTree *gtree = new RootedTree(strtree);
                gtree->setid(startid++);
                genetreesv->push_back(gtree);
            }
                
        }
        else        
        {
            cerr << "Vector of values or a string expected in genetrees assignment. Found " << ev->type << endl;
            exit(-1);
        }     

        // Initialize treespace
        if (globaltreespace->used())
        {

            cerr << "Cannot initialize tree space with gene trees (already used). Set all gene trees before cost computations" << endl;
            exit(-1);
        }
        globaltreespace = new TreeSpace(env.genetreesv, maxdisplaytreecachesize);


    }
    else if (label=="setspecies")
    {
        setspecies(ev->getint("setspecies"));
    }
    else if (label=="savewhenimproved")
    {
        
        env.getglobalstats()->setsavewhenimproved(ev->getint("setspecies"));
        env.setvar(label, ev);         
    }
    else
    {
        env.setvar(label, ev);      
    }
    return ev;
}


ostream& operator<<(ostream& os, const Env& e)
{
    os << "ENV" << endl;    
    os << ":genetreesv:" << e.genetreesv.size() << endl;
    // os << ":printstats:" << e.printstats << endl;
    // os << ":usenaive:" << e.usenaive << endl;
    // os << ":runnaiveleqrt:" << e.runnaiveleqrt<< endl;
    // os << ":cutwhendtimproved:" << e.cutwhendtimproved << endl;
    // os << ":networkclass:" << e.networkclass << endl;
    // os << ":timeconsistency:" << e.timeconsistency << endl;
    // os << ":reticulations:" << e.reticulations << endl;
    // os << ":randnetuniform:" << e.randnetuniform << endl;
    // os << ":randomnetworks:" << e.randomnetworks << endl;
    // os << ":quasiconsensusnetworks:" << e.quasiconsensusnetworks << endl;
    if (e.guidetree) os << ":guidetree:present" << endl;
    if (e.guideclusters) os << ":guideclustes:present" << endl;
    for (auto const& [key, val] : e.vars)
    {
        os << ":::" << key << ':'  << *val << std::endl;
    }
    return os;
}

