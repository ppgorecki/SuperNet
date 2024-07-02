

#include <iostream>
#include <regex>
#include <string>
#include <cctype>
#include "netgen.h"
#include "network.h"

#define ELABEL 1
#define ENUM 2
#define ERNUM 3
#define ESTR 4
#define ENETS 5
#define ELIST 6
#define EVALUEVEC 7
#define EMAIN 8
#define EASSIGNMENT 9
#define EFOR 10
#define ETO 11
#define EGTREES 12

#define ECALL 16
#define ECLIMB (1|ECALL)
#define EITALGORITHM (2|ECALL)
#define ESTOREBEST (3|ECALL)
#define EPRINT (4|ECALL)
#define EPRINTENV (5|ECALL)
#define EREAD (6|ECALL)
#define ERANDTREES (8|ECALL)
#define EQUASICONSTREES (9|ECALL)
#define EADDRETICULATIONS (10|ECALL)
#define EPRINTLN (12|ECALL)
#define ESHOWSTATS (13|ECALL)
#define ENETLIST (14|ECALL)

using namespace std;

struct Env;

class AlgTokenizer 
{

    public:

        int linenum;
        float floatval;
        int intval;
        int m_tokentype;

        static const std::string DELIMITERS;
        AlgTokenizer(const std::string& str);
        AlgTokenizer(const std::string& str, const std::string& delimiters);
        int next();
        int next(const std::string& delimiters) {
            m_tokentype = _next(delimiters);
#ifdef _DEBUG_ALGTOK_
            cout << "NEXT:" << m_tokentype << " off=" << m_offset << " m_token=<" << m_token << ">" << endl;
#endif
            return m_tokentype;
        }   
        int expect(const std::string);
        int _next(const std::string& delimiters);
        int gettokentype() { return m_tokentype; }
        string get() {return m_token; }
        const std::string get() const;
        std::string prev() { return m_prev; }
        void Reset();

    protected:
        size_t m_offset;
        const std::string m_string;
        std::string m_prev;
        std::string m_token;
        std::string m_delimiters;
};


class ExprValue;

class Expr
{    
public:
    int type;        
    Expr(int t) : type(t) {}
    virtual bool isvalue() { return false; }
    virtual ExprValue *eval(Env &env) { return NULL; }
};

class ELabel; 

class ExprList : public Expr
{    
public:    
    vector<Expr*> t;
    ExprList() : Expr(ELIST), t() {}    
    ExprValue *eval(Env &env);
    void append(Expr *e);
};

class ExprValue:
 public Expr
{

public: 
    virtual void print(ostream &os) const {}
    ExprValue(int type) : Expr(type) {}
    bool isvalue() { return true; }
    virtual ExprValue *eval(Env &env) { return this; }
    friend ostream& operator<<(ostream& os, ExprValue const& e) { e.print(os); return os; } 
    virtual string getstr(string comment)
    {
        cerr << "String expected in " << comment << ". Found " << type << endl;
        exit(-1);
    }
    virtual int getint(string comment)
    {
        cerr << "Int expected in " << comment << ". Found " << type << endl;
        exit(-1);
    }
};

class VInt : public ExprValue
{
public: 
    int v;
    VInt(int _v) : ExprValue(ENUM), v(_v) {}  
    int getint(string comment) { return v; } 
    void print(ostream &os) const { os << v; };
};

class VFloat : public ExprValue
{

public:
    float v;
    VFloat(float _f) : ExprValue(ERNUM), v(_f) {}    
    void print(ostream &os) const { os << v; };
};

class EValueVec : public ExprValue
{
protected:
    vector<ExprValue*> t;

public:
    void print(ostream &os, bool sepln=false) const 
    {
        bool prev=false;    
        if (!sepln) os << "[";
        for (auto v: t)
        {                    
            if (!sepln && prev) os << ", ";
            if (v)
            {
                os << *v;                
            }
            else os << "NULL";
            if (sepln) os << endl;
            prev = true;              
        }      
        if (!sepln) 
            os << "]";        
    };
    EValueVec(vector<ExprValue*> &r) : ExprValue(EVALUEVEC), t(r) {} 
    ExprValue *get(int i) { return t[i]; }
    int size() { return t.size(); }
};


class VStr : public ExprValue
{

public:
    void print(ostream &os) const { os << s; };
    string s;
    VStr(const string &_s) : ExprValue(ESTR), s(_s) {}
    friend ostream& operator<<(ostream& os, const VStr& e) { return os << e.s; }
    string getstr(string comment) { return s; }
};

class ELabel: public Expr
{
public:
    string s;
    ELabel(const string &_s) : Expr(ELABEL), s(_s) {}    
    ExprValue *eval(Env &env);    
};  

typedef Network* PNetwork;

class VGeneTrees: public ExprValue
{

    public:    
      void print(ostream &os) const { for (auto gt: *gtrees) os << (*gt) << endl; };

      VecRootedTree *gtrees;    
      VGeneTrees(VecRootedTree *gtrees) : ExprValue(EGTREES), gtrees(gtrees) {}
};


class VNetGen: public ExprValue
{
protected: 
public:    
    void print(ostream &os) const {         
        Network *n; 
        while ((n = netgen.next())!=NULL)      
            os << (*n) << endl;         
        netgen.reset();
    }

    struct Iterator 
    {
        using iterator_category = std::forward_iterator_tag;
        using difference_type   = std::ptrdiff_t;
        using value_type        = PNetwork;
        using pointer           = PNetwork*;
        using reference         = PNetwork&;

        Iterator(NetGen &netgen) : netgen(netgen) { current = netgen.next(); }
        Iterator(NetGen &netgen, int) : netgen(netgen) { current = NULL; }
        
        reference operator*() { return current; }
        pointer operator->() { return &current; }
        Iterator& operator++() { current = netgen.next(); return *this; }  
        Iterator operator++(int) { Iterator tmp = *this; ++(*this); return tmp; }
        friend bool operator== (const Iterator& a, const Iterator& b) { return a.current == b.current; };
        friend bool operator!= (const Iterator& a, const Iterator& b) { return a.current != b.current; };  

    private:
        PNetwork current;
        NetGen &netgen;    
    };

    Iterator begin() { netgen.reset(); return Iterator(netgen); }
    Iterator end()   { return Iterator(netgen, 0); }

    NetGen &netgen;    
    bool cleanunused;

    VNetGen(NetGen &netgen, bool cleanunused=false) : ExprValue(ENETS), netgen(netgen), cleanunused(cleanunused) {}

    VNetGen(VecNetwork *nets, bool cleanunused=false) : VNetGen(*(new NetVecGenerator(*nets)), cleanunused) {}
    
    // VNetGen(VecNetwork *netvec) : ExprValue(ENETGENERATOR), cleanunused(false), netgen(*(new NetGen(netvec))) {}
};


class ExprAssignment : public Expr
{
public:
    string label;
    Expr *value;    
    ExprAssignment(string elabel, Expr *evalue) : 
        Expr(EASSIGNMENT), label(elabel), value(evalue) {}
    ExprValue *eval(Env &env);
};

class ExprFor : public Expr
{
public:
    ExprAssignment *init;    
    Expr *to;
    ExprList *block;
    ExprFor(ExprAssignment *init, Expr *to, ExprList *block) : 
        Expr(EFOR), init(init), to(to), block(block) {}
    ExprValue *eval(Env &env);
};

class ExprCall : public Expr
{    
public:
    ExprList *args;
    ExprCall(int fun, ExprList *eargs) : Expr(fun), args(eargs) {}
    ExprValue *eval(Env &env);    
};


ExprList* parsescript(string script);

struct Env
{
    VecRootedTree &genetreesv;                  
    std::map<std::string, struct ExprValue*> vars;
    CostFun *costfun;       
    ClimbStatsGlobal *globalstats;

    Clusters *guideclusters;
    Clusters *guidetree;
    RootedTree *preserverootst;       

    Env *parent;
    bool rootenv;

    Env(VecRootedTree &genetreesv, 
        VecNetwork &networksv,
        CostFun *costfun,        
        Clusters *guideclusters,
        Clusters *guidetree,
        RootedTree *preserverootst        
        ) :
    genetreesv(genetreesv),    
    costfun(costfun),    
    guideclusters(guideclusters),
    guidetree(guidetree),
    preserverootst(preserverootst),
    parent(NULL),
    globalstats(NULL)
    {
        vars["NET_GENERAL"]=new VInt(NET_GENERAL);
        vars["NET_CLASS1RELAXED"]=new VInt(NET_CLASS1RELAXED);
        vars["NET_TREECHILD"]=new VInt(NET_TREECHILD);
        vars["NET_ANY"]=new VInt(NET_ANY);
        vars["NET_TIMECONSISTENT"]=new VInt(NET_TIMECONSISTENT);
        vars["NET_NOTIMECONSISTENT"]=new VInt(NET_NOTIMECONSISTENT);
        vars["networks"] = new VNetGen(&networksv);
        vars["genetrees"] = new VGeneTrees(&genetreesv);
    }

    ExprValue* get(string label, ExprValue *def=NULL)
    {
        if (vars.find(label)==vars.end())
        {
            return def; 
        }
        return vars[label];
    }

    void setparentenv(Env &env)
    {
        parent = &env;
    }

    int getint(string info, ExprValue *v)
    {
        if (v->type!=ENUM)
        {
            cerr << "Int value expected in " << info << ". Found " << v->type << endl;
            exit(-1);
        }
        return ((VInt*)v)->v;
    }

    void set(string label, int v)
    {
        setvar(label, new VInt(v));
    }

    void set(string label, string v)
    {
        setvar(label, new VStr(v));
    }

    bool has(string label)
    {
        return vars.find(label)!=vars.end();
    }

    VecRootedTree *getgenetrees()
    {
        VecRootedTree *vg = NULL;        
        if (has("genetrees"))
        {
            VecRootedTree *vg = ((VGeneTrees*)vars["genetrees"])->gtrees;
            if (vg!=&genetreesv)
            {
                cerr << "Warning: gene trees vectors are inconsistent in Env" << endl;                
            }            
            return vg;
        }
        return vg;
    }


    void setvar(string label, ExprValue *v)
    {
        vars[label] = v; // Always add
    }

    void delvar(string label)
    {
        vars.erase(label);   
    }

    int getint(string label, int def=0)
    {
        if (vars.find(label)==vars.end())
        {
            return def; 
        }
        ExprValue *e = vars[label];
        if (e->type != ENUM)
        {
            cerr << "Expected int. Found type " << e->type << endl;
        }
        return ((VInt*)e)->v;
    }

    int getfloat(string label, float def=0)
    {
        if (vars.find(label)==vars.end())
        {
            return def; 
        }
        ExprValue *e = vars[label];
        if (e->type == ENUM)
        {
            return ((VInt*)e)->v;
        }
        if (e->type == ERNUM)
        {
            return ((VFloat*)e)->v;
        }
        cerr << "Expected float or int. Found type " << e->type << endl;
        exit(-1);
    }


    string getstr(string label)
    {
        if (vars.find(label)==vars.end())
        {
            cerr << "Variable " << label << " undefined" << endl;
            exit(-1);
        }
        ExprValue *e = vars[label];
        if (e->type != ESTR)
        {
            cerr << "Expected string. Found type " << e->type << endl;
            exit(-1);
        }
        return ((VStr*)e)->s;
    }


    VNetGen* getsource(string info, bool allownull=false)
    {
        ExprValue *src = get("source", NULL);
        if (allownull && !src) return NULL;

        if (!src)
        {
            cerr << "source arg missing in" << info << endl;
            exit(-1);
        }        
        if (src->type==ENETS)
        {
            return (VNetGen*)src;
        }
        cerr << "Networks expected in source variable in " << info << endl;
        exit(-1);

    }

    ClimbStatsGlobal* getglobalstats()
    {
        if (parent) return parent->getglobalstats();

        int networkclass = getint("networkclass");
        unsigned int randseed = getint("randseed");
        int timeconsistency = getint("timeconsistency");
        int savewhenimproved = getint("savewhenimproved");
        
        if (globalstats==NULL)
        {
            // at main level
            string outfiles = getstr("outfiles");
            string outdirectory = getstr("outdirectory");
            int odtlabelled = getint("odtlabelled");
            globalstats = new ClimbStatsGlobal(networkclass, timeconsistency, new DagSet, randseed, savewhenimproved);
            if (outfiles.length()) 
            {
                globalstats->setoutfiles(outdirectory, outfiles, odtlabelled);
            }
        }

        if (networkclass!=globalstats->getnetworkclass())
        {
            cerr << "Warning. Networks class is inconsistent with the global setting" << endl;
        }        
        if (timeconsistency!=globalstats->gettimeconsistency())
        {
            cerr << "Warning. Time consistency is inconsistent with the global setting" << endl;
        }        
        return globalstats;
    }

    void finalizestats()
    {
        if (globalstats)
        {            
            int verbosealg = getint("verbosealg");
            vector<ClimbStatsGlobal*> globalstatsarr; 
            // merge all data and save odt/dat file(s); optional     
            globalstats->savedatmerged(verbosealg >= 4, globalstatsarr, true);     

            // best dags to file
            globalstats->savebestdags(verbosealg >= 4, true);

            // print summary
            globalstats->print();
            globalstats->printnetworkinfo();     
            cout << endl;
            // globalstats = NULL;
        }
    }

    friend ostream& operator<<(ostream& os, const Env& e);

  
};


