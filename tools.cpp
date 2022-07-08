

/************************************************************************
 SuperNetwork Inference - based on FastUrec project

 *************************************************************************/

#include <ctype.h>
#include <iostream>
#include <cstdlib>
#include <cstring>
#include <stdlib.h>

using namespace std;

#include "tools.h"

vector<string> specnames;
vector<int> specorder;
map<string, int> specnames2id;
map<SPID, SPID*> spec2spcluster;

int linenum = 0;
int charnum = 0;
char *inputfilename = (char*)NULL;


void printspcluster(ostream&s, SPID *a) {
  s << a[0] << ":";
  for (int i = 1; i <= a[0]; i++) s << a[i] << " ";
}

SPID topspcluster[2] = {0, 0};
 
// Join two sets of species
// Each set is an array where the first item (at index 0) is the size of the set
// Species ids are sorted
// Special case: topspcluster where all species names are present

SPID* joinspclusters(SPID* a, SPID* b, SPID *dest)
{
  if (!a[0]) return a; // topcluster
  if (!b[0]) return b; // topcluster

  int as = a[0], bs = b[0], cs = 1;  
  
  SPID r[as+bs+1];

  int i = 1, j = 1;

  while (i <= as && j <= bs)
    if (a[i] == b[j]) { r[cs++] = a[i];  i++; j++; }
    else if (a[i] < b[j]) r[cs++] = a[i++];
    else if (a[i] > b[j]) r[cs++] = b[j++];

  while (i <= as) r[cs++] = a[i++];
  while (j <= bs) r[cs++] = b[j++];

  if (cs == 2) 
    return spec2spcluster[r[1]];

  if (cs - 1 == specnames.size())  // all species -> topcluster   
    return topspcluster;

  SPID *rx = dest;
  if (!dest) rx = new SPID[cs];

  for (i=1; i<cs; i++) rx[i] = r[i];
  rx[0] = cs - 1;  // store size

  return rx;
}


bool spsubseteq(SPID *a, SPID *b)
{
  if (b==topspcluster) return 1;
  if (a==topspcluster) return b==topspcluster;
  if (spsize(a)>spsize(b)) return false;

  int jp=1;
  for (int i = 1; i <= a[0]; i++)
  {
    while (a[i]>b[jp]) 
      { jp++; 
        if (jp>b[0]) return false;  //last
      } // not found
    if (a[i]<b[jp]) return false;
  }
  return true;
}

void initlinenuminfo(char *fname)
{
  linenum = 1;
  charnum = 0;
  inputfilename = fname;
}

SPID* spidcopy(SPID* a, int size)
{
  SPID *res = new SPID[size];
  for (int i = 0; i < size; i++) res[i] = a[i];
  return res;
}

void printlinepos()
{
  if (inputfilename) cerr << inputfilename << ":";
  cerr << linenum << ":" << charnum << ": ";
}

void expectTok(const char *tok, char *s, int p)
{
  printlinepos();
  cerr << "parse error - expected #" << tok << "#" << " after " << endl;
  int x = p;
  if (x < 70) x = 0;
  else { x -= 70; cerr << "..."; }
  while (s[x] && x < p + 70)
  {
    if (x == p) cerr << endl << "Found: " << endl;
    cerr << s[x++];
  }
  if (s[x]) cerr << "...";
  cerr << endl;
  exit(-1);
}

char* getTok(char *s, int &p, int num)
{
#define inctok { p++; if (s[p]=='\n') { linenum++; charnum=0; } charnum++; }
  while (isspace(s[p])) inctok;

  char *cur = s + p;
  if (isalnum(s[p]) || (s[p] == '_') || (s[p] == '#'))
  {
    if (s[p] == '#') inctok;
    while (isalnum(s[p]) || (s[p] == '_') || (s[p] == '.') || (s[p] == '-')) inctok;
    return cur;
  }
  if ((s[p] == '(')  || (s[p] == ')') || (s[p] == ',') || s[p] == ':')
  {
    inctok;
    return cur;
  }
  if (num >= 0) cerr << "Tree nr " << num << ". ";
  printlinepos();
  cerr << " parse error - expected expected LABEL, '(', ')' or ','; found: '" << s[p]
       << "'. current string: <<" << cur << ">>" << endl;
  exit(-1);
}


void checkparsing(char *s)
{
  int p = 0;
  while (s[p])
  {
    if (s[p] == ';') { inctok; continue; }
    if (!isspace(s[p]))
    {
      printlinepos();
      cerr << "EOLN, EOF or ';' expected, found '" << (s + p) << "'" << endl;
      exit(-1);
    }
    inctok;
  }
}


string getstringn(char *s, int len)
{
  char old;  
  if (len)
  {
    old = s[len];
    s[len] = 0;
  }
  string t = s;
  
  if (len) s[len] = old;
  return t;
}

int getspecies(char *_s, int len)
{
  string s = getstringn(_s,len);

  if (specnames2id.count(s)) return specnames2id[s];

  int num;

  specnames2id[s] = num = specnames.size();
  specorder.push_back(num);
  specnames.push_back(s);    
  SPID* speccluster = new SPID[2];
  speccluster[0] = 1;
  speccluster[1] = num;
  spec2spcluster[num] = speccluster;
  return num;
  
}

void cleanspecies()
{
    int i;
    for (SPID i = 0; i < specnames.size(); i++)
      delete[] spec2spcluster[i];
}


char *mstrndup(const char *s, size_t n)
{
  char *result;
  size_t len = strlen (s);
  if (n < len) len = n;
  result = (char *) malloc (len + 1);
  if (!result) return 0;
  result[len] = '\0';
  return (char *) memcpy (result, s, len);
}


void randomizearr(SPID *t, int size)
{
  for (int i = 0; i < size; i++)
  {
    int np = rand() % size;
    int tmp = t[i];
    t[i] = t[np];
    t[np] = tmp;
  }
}



double gettime()
{
    struct timespec now;
    clock_gettime(CLOCK_REALTIME, &now);
    return now.tv_sec + now.tv_nsec*1e-9;
}



extern const char *SUPNET;

int usage(int argc, char **argv) {

  cout <<
       "SuperNetwork " << SUPNET << "\n"
       "Software homepage: http:...\n"
       "Usage: " << argv[0] << " [options]\n\n"
       "INPUT OPTIONS:\n"
       "  -g STR: gene trees from a string (separator is semicolon)\n"
       "  -s species tree\n"
       "  -n network\n"
       "  -G filename - gene trees from a file; EOLN separated; - stdin\n"              
       "  -S filename - species trees from a file; EOLN separated; - stdin\n"
       "  -N filename - networks from a file; EOLN separated; - stdin\n"

       "  -l [aDELIM|bDELIM|pPOS] - rules for matching gene tree names with species tree names;\n"       
       "       a - after delimiter; b - before delimiter; p - from position \n"
       
       "\n"
       
       "Random and quasi consensus generators\n"
       "  -r NUM - gen NUM random species trees/networks\n"       
       "  -q NUM - gen NUM quasi-consensus species trees/networks\n"              
       "  -R NUM - inserts NUM reticulation nodes in all networks/trees from -q, -r, -n, -N (default 0)\n"
       "  -A SPECIESNUM - define SPECIESNUM species a,b,...\n"       
       
       "  By default generator produces tree-child networks. Use -e1, -e2 for other classes (see below)\n"
       "\n"


       "ADDITIONAL OPTIONS\n"
       "  -v NUM - verbose level 1 (print the name of outputfile only),\n"
       "       2 (simple output), 3 (detailed output)\n"       
       "  -e [gsrD...]+ - multiple options:\n"
       "     g - print a gene tree\n"
       "     s - print a species tree\n"
       "     n - print a network\n"
       "     t - print all display trees\n"
       "     T - print all display trees with their ids\n"
       "     o - min total cost of all gene trees vs network using naive approach\n"
       "     r - preserve root when searching the species tree space and in quasi-consensus\n"              
       "     D - detailed tree info\n"       
       "     i - list species dictionary\n"       
       "     c - print cost\n"       
       "     C - print cost with trees\n"              
       "     d - run DP algorithm only (approx DC) per each (G,N)\n"
       "     b - run BB algorithm only (exact DC); see jJ subopt.\n"
       "     1 - class 1 networks in rand generator (int node has at most one reticulation child); see -R\n"       
       "     2 - general phylogenetic networks (no restrictions); see -R\n"
       "     p - pairwise comparison of all networks\n"
       "     u - print unique networks; summary stats printed on stderr\n"
       "     U - similar to u plus counts of networks\n"
       "     x - two networks are equal if their shapes are isomorphic in u,U,p (i.e., ignore leaf labels)\n"
       "     L - for each v in V(N), print the number of nodes reachable from v (only from networks)\n"
       "     l - for each v in V(N), print the number of leaves reachable from v (only from networks)\n"       
       "\n"               
       "COST SETTING OPTIONS\n"
       // " TODO: -D dupweight  - set weight of gene duplications (def. 1.0)\n"
       // " TODO: -L lossweight - set weight of gene losses (def. 1.0)\n" 
       "  -C COST - set cost function from {DL,D,L,DC,RF}\n"
       "\n"
       
       "BB algorithm for DC (with -eb)\n"      
       "  -t THRESHOLD - run naive odt computation, when the number of reticulations is < THRESHOLD; otherwise run DP;"
       "  -ebj - gen bb.tsv with stats\n"
       "  -ebJ - gen bb.dot with bb tree search\n"
       "  -ebk - print time per each pair\n"
       "\n"
       "\n"

       "ODT HEURISTIC SEARCH\n"

       "  -o [TNt123sq]+ - run hill climbing heuristic using cost function and print optimal cost, non TC networks are allowed, tail/nni moves, all optimal networks are written in odt.log file; summary stats are save to odt.dat (opt cost, total time in sec., time of hill climbing, time of merge step, number of networks, improvements, steps and the number of starting networks)\n"
       "   By default each HC step starts from the set of predefined networks (if -n or -N are provided); otherwise by using quasi-random networks\n"
       "       T - use TailMoves (default)\n"
       "       N - use NNI instead of TailMoves\n"      
       "       t - TailMoves limited to tree-child (with 3a-3c conditions)\n"
       "       1 - print visited network after each improvement (strictly)\n"
       "       2 - print visited networks if the cost is equal to the current or improved\n"
       "       3 - print all visited networks\n"
       "       S - print stats after locating optimal networks after each HC runs\n"
       "       s - print extended stats after each HC run\n"
       "       q - do not save odt.log (and odt.dat) with optimal networks\n"
       

       "  -K NUM - stopping criterion: stop when there is no new network after NUM HC runs\n"

       "\n"
       "  -O ODTFILE - change the name of odt.log and odt.dat files\n"
              
       "\n"
       // "*Starting tree options:\n"
       // "  -q NUM - generate NUM quasi-consensus trees sampled from gene tree clusters;\n"
       // "       (def. is 1 if no user species tree is defined)\n"
       // "  -r NUM - generate NUM random trees (def. 0)\n"
       // "  -Q p - defines probability for sampling gene trees (use with -q, default 0.2)\n"
       // "  See also -qNUM, -ei and -er.\n"

       "Examples: \n\n"

       "Print 10 quasi consensus trees\n"
       " supnet -g '(a,((b,c),d));(a,(b,d))' -q10 -en\n"
       
       "Print 10 quasi consensus trees with preserved split of the root (-er)\n"
       " supnet -g '(a,((b,c),d));(a,(b,e))' -s'((a,b),(c,(d,e)))' -q10 -enr\n"

       "Detailed tree/network info\n"
       "  supnet -g '(a,((b,a),c))' -eD\n"

       "Print DC cost\n"
       "  supnet -g '(b,(a,c))' -s '(a,(b,c))' -CDC -ec \n"

       "Print 10 random species trees over a..e (5)\n"
       "  supnet -A5 -r10 -en\n"

       "Print display trees (based on reticulation switching; trees maybe non-unique)\n"
       "  supnet -n '((((c)#B,b))#A,(#A,(#B,a)))' -et\n"

       "Print display trees with ids\n"
       "  supnet -n '((((c)#B,b))#A,(#A,(#B,a)))' -eT\n"

       "Print min total cost 10 random gene trees vs random tree-child network with 5 reticulations over 8 species; print the initial network\n"
       "  supnet -r10 -A8 -en  | supnet -G- -r1 -A8 -R5 -eon"
       //"  supnet -r10 -A10 -en | supnet -G- -n $( embnet.py -n \"rand:10:5\" -pn ) -eo\n"

       "Printing and visualizing\n"
       "  supnet -n '((a)#A,(#A,(c,b)))' -enD\n"       
       "  supnet -n '((((c)#B,b))#A,(#A,(#B,a)))' -enD\n"       
       "  supnet -n '(((b)#A,a),(#A,c))' -d\n"
       "  supnet -n '((a)#A,(#A,(c,b)))' -d | dot -Tpdf > n.pdf\n"
       "  supnet -n '((((c)#B,b))#A,(#A,(#B,a)))' -d | dot -Tpdf > n.pdf\n"
       "  supnet  -r1 -R10  -A20 -d | dot -Tpdf > n.pdf\n"
       "  N=$( embnet.py -n 'rand:3:1' -pn ); echo $N; supnet -n$N -d | dot -Tpdf > n.pdf\n"

       "\nHILL CLIMBING HEURISTIC (-o...)\n"

       " Minimalistic run: print cost; result in odt.log; one random initial network (-r1, i.e., one HC climb) with 2 (-R2) reticulations\n"
       "  supnet -g \"(a,(b,(c,d))); ((a,b),(c,d))\" -r1 -R2 -oT\n"
       
       " Print summary stats; 10 HC runs with random networks\n"
       "  supnet -g \"(a,(b,(c,(d,e)))); ((a,b),(c,(e,a))); ((b,c),(d,a))\" -r10 -R3 -o2S\n"

       " Print cost, improvements and stats (s); tree-child search (t); quasi consensus initial network (-q) with 3 random reticulations\n"
       "  supnet -g \"(a,(b,(c,(d,e)))); ((a,b),(c,(e,a))); ((b,c),(d,a))\" -q1 -R3  -o2st\n"

       " Larger instance; tree-child search:\n"
       "  supnet -g \"((i,c),((a,d),(f,(b,((g,(e,j)),h)))));((i,h),((c,f),((a,(d,e)),(j,(g,b)))));(((f,(b,d)),((j,g),(e,i))),(c,(h,a)));(((f,(i,(j,d))),(c,b)),(((h,a),g),e));((((h,e),((c,f),a)),(d,b)),((g,i),j))\" -R8 -q1 -o2st\n"

       " Print only improvements:\n"
       "  supnet -g \"(a,(b,(c,d))); ((a,b),(c,d))\" -r1 -R3 -o1\n"
       
       " Print only improvements and equal cost networks:\n"
       "  supnet -g \"(a,(b,(c,d))); ((a,b),(c,d))\" -r1 -R3 -o2\n"

       " Print improvements and equal cost networks; NNI moves:\n"
       "  supnet -g \"(a,(b,(c,d))); ((a,b),(c,d))\" -r1 -R3 -o3N\n"
       
       " Print improvements; skip odt.log:\n"
       "  supnet -g \"(a,(b,(c,d))); ((a,b),(c,d))\" -r1 -R3 -o3Nq\n"

       " Recommended with large HC-runs using quasi-consensus rand networks\n"
       "  supnet -g \"(a,(b,(c,d))); ((a,b),(c,d))\" -q1000 -R3 -oS\n"

       " Using stoppning criterion (-K1000) with quasi-consensus rand networks\n"
       "  supnet -g \"(a,(b,(c,d))); ((a,b),(c,d))\" -q-1 -R3 -oS -K1000\n"

       " Display trees usage stats:\n"
       "  supnet -N odt.log -et | sort | uniq -c | sort -k1 -n\n"

       " Insert 2 reticulations into a network (tree-child output)\n"
       "  supnet -R2 -n '(a,((d)#1,(b,(c,#1))))' -en\n"

       " Insert 10 reticulations into a network (general network)\n"
       "  supnet -R10 -n '(a,((d)#1,(b,(c,#1))))' -en2\n"

       " Pairwise comparison of random dags (general dags)\n"
       "  supnet -r10 -R1 -A2 -en2p\n"

       " Print unique networks from odt.log (-eu)\n"
       "  supnet -g '(a,(b,(c,d))); ((a,b),(c,d))' -r1 -R3 -o2 && supnet -eu -N odt.log\n"

       " Print unique random networks with counts\n"
       "  supnet -r100000 -R1 -A3 -eU\n"

       " Print unique random shapes of networks with counts\n"
       "  supnet -r100000 -R1 -A3 -eUx\n"

       ;


       

  exit(-1);
}


// return random tree string from given vector of species
string genrandomtree(SPID *sp, int len)
{
  vector<string> v;
  int i;
  for (i = 0; i < len; i++)
    v.push_back(species(sp[i]));

  while (v.size() > 1)
  {
    int a = rand() % v.size();
    string as = v[a];
    v.erase(v.begin() + a);

    int b = rand() % v.size();
    string bs = v[b];
    v.erase(v.begin() + b);

    v.push_back("(" + as + "," + bs + ")");
  }
  return v[0];
}

int strcount(char *s, char c)
{  
  int cnt = 0;  
  while (*s) 
  {
    if (*s == c) cnt++;
    s++;
  }
  return cnt;
}


void shuffle(SPID *a, size_t n)
{
    if (n > 1) 
    {
        size_t i;
        for (i = 0; i < n - 1; i++) 
        {
          size_t j = i + rand() / (RAND_MAX / (n - i) + 1);
          SPID t = a[j];
          a[j] = a[i];
          a[i] = t;
        }
    }
}
