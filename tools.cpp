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
#include "bb.h"

vector<string> specnames;
vector<int> specorder;
map<string, int> specnames2id;
map<NODEID, NODEID*> spec2spcluster;

int linenum = 0;
int charnum = 0;
char *inputfilename = (char*)NULL;

void printspcluster(ostream&s, NODEID *a) {
  s << a[0] << ":";
  for (int i = 1; i <= a[0]; i++) s << a[i] << " ";
}

NODEID topspcluster[2] = {0, 0};
 
// Join two sets of species
// Each set is an array where the first item (at index 0) is the size of the set
// Species ids are sorted
// Special case: topspcluster where all species names are present

NODEID* joinspclusters(NODEID* a, NODEID* b, NODEID *dest)
{
  if (!a[0]) return a; // topcluster
  if (!b[0]) return b; // topcluster

  int as = a[0], bs = b[0], cs = 1;  
  
  NODEID r[as+bs+1];

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

  NODEID *rx = dest;
  if (!dest) rx = new NODEID[cs];

  for (i=1; i<cs; i++) rx[i] = r[i];
  rx[0] = cs - 1;  // store size

  return rx;
}


bool spsubseteq(NODEID *a, NODEID *b)
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

NODEID* spidcopy(NODEID* a, int size)
{
  NODEID *res = new NODEID[size];
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
  NODEID* speccluster = new NODEID[2];
  speccluster[0] = 1;
  speccluster[1] = num;
  spec2spcluster[num] = speccluster;
  return num;
  
}

void cleanspecies()
{
    int i;
    for (NODEID i = 0; i < specnames.size(); i++)
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


void randomizearr(NODEID *t, int size)
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
       
       "SupNet: Phylogenetic networks visualisation, evaluation and inference\n"

       "SuperNetwork " << SUPNET << " with "
#if defined USE_PRIORITY_QUEUE_MINRT
       "MinRT-strategy\n"
#endif
#if defined USE_PRIORITY_QUEUE_MINLB
       "MinLB-strategy\n"
#endif 
#if defined USE_QUEUE_BFS
      "BFS-strategy\n"
#endif      

       "Software homepage: https://github.com/ppgorecki/SuperNet\n"
       
       "Usage: " << argv[0] << " [options]\n\n"

       "\n"

       "Input:\n"
       " --genetrees, -g STR: gene trees from a string (separator is semicolon)\n"
       " --speciestrees, -s STR: species trees from a string\n"
       " --networks, -n STR: networks from a string\n"

       " --genetreesfile, -G FILE: gene trees from a file; EOLN separated; - stdin\n"              
       " --speciestreesfile, -S FILE: species trees from a file; EOLN separated; - stdin\n"
       " --networksfile, -N filename: networks from a file; EOLN separated; - stdin\n"

       "\n"

       "Print:\n"
       " --pgenetrees: print gene trees\n"
       " --pspeciestrees: print species trees\n"
       " --pnetworks: print networks\n"
       " --pdisplaytrees: print all display trees\n"
       " --pdisplaytreesext: print all display trees with their ids\n"
       " --dot: print dot representation of networks\n"

       "Matching gene leaf labels with species leaf labels:\n"

       " [TODO] --matchparam, -l [NUM|STR]: parameter for the occurrence of species name in gene label; see --match* options below\n"
       " [TODO] --matchafter: species name is after --matchparam (STR) in gene label\n" 
       " [TODO] --matchbefore: species name is before --matchparam (STR) in gene label\n" 
       " [TODO] --matchatpos: species name is at fixed --matchparam (NUM) position in gene label; negative NUM positions from the end of label\n" 

       //"  -l [aDELIM|bDELIM|pPOS] - rules for matching gene tree names with species tree names;\n"       
       // "       a - after delimiter; b - before delimiter; p - from position \n"
       
       "\n"

       "Binary network classes in random generator and HC heuristic:\n"
       " --treechild: tree-child networks (default); a non-leaf node has at least one tree node child \n"
       " --timeconsistent: time-consistent networks \n"              
       " --notimeconsistent: no time-consistent networks \n"              
       " --relaxed: a node has at most one reticulation child; see -R\n"       
       " --general: general phylogenetic networks (not applicable with DP algorithm) \n"
       " --detectclass: for each network N print 't tc r N', where t is 1 if N is time consistent, tc is 1 if N is tree-child network, r is 1 if N is relaxed network.\n"
       
       "\n"       
       "Random and quasi consensus tree/network generators\n"
       " -r NUM: generate NUM random networks; if NUM=-1 generator is infinite\n"       
       " -q NUM: generate NUM quasi-consensus networks; if NUM=-1 generator is infinite\n"              
       " -R NUM: insert NUM reticulation nodes into networks from -q, -r, -n, -N; default is 0\n"       
       " -A SPECIESNUM: define SPECIESNUM species a,b,... [TODO: expand label range]\n"       
       " --uniformedgesampling: draw uniformly a pair of edges to create random reticulation; the default: first draw the source edge, then the uniformly the destination edge given the source\n" 
       " --preserveroot: preserve root in quasi-consensus and in HC algorithm\n"              
       " --randseed, -z SEED: set seed for random generator (srand)\n"       
       
       // " S - print species trees using compressed representation\n" 

       "\n"
       "Tree/network print:"
       " --pgenetrees: print gene trees\n"
       " --pspeciestrees: print species trees\n"
       " --pnetworks: print networks\n"
       " --pdisplaytrees: print all display trees\n"
       " --pdisplaytreesext: print all display trees with their ids\n"
       " --dot: print dot representation of networks\n"

       "\n"

       "Cost functions:\n"
       " --cost, -C COST: set cost function from {DL,D,L,DC,RF,DCE}; default is DC\n"
       " [TODO] -D dupweight: set weight of gene duplications (def. 1.0)\n"
       " [TODO] -L lossweight: set weight of gene losses (def. 1.0)\n" 

       "\n"

       "Tree(s) vs tree(s) cost computation:\n"
       " --ptreecost: print cost between two trees (G,S)\n"       
       " --ptreecostext: print cost between two trees with trees (G,S)\n"              
       
       "\n"

       "DP algorithm to compute lower bound of DC between a gene tree and a network by dynamic programming\n"
       "--DP: run DP and print the bound\n"            

       "\n"

       "ODT-naive algorithm via enumeration of display trees to compute exact DC cost between gene trees and a network:\n"
       " --odtnaivecost: run ODT-naive and print the cost\n"
       " [TODO] --odtnaivesampling: set sampling in ODT naive (exponential distribution)\n"       
       " --ptreesinodtnaive: print additional details\n"

       "\n"

       "Branch and bound (BB) algorithm to compute exact DCE cost by multiple calls of ODT-naive and/or DP\n"
       " --BB: compute exact DCE by BB\n"      
       " -t THRESHOLD, --runnaiveleqrt: run ODT-naive cost computation, when the number of reticulations is <= THRESHOLD; otherwise run DP;\n"
       " --bbtsvstats: gen bb.tsv with stats\n"
       " --bbtreesearch: gen bb.dot with bb tree search\n"
       " --bbtimestats: print time per each pair\n"
       " --bbstartscore=FLOAT: define initial score in BB (testing only)\n"
                          
       "\n"

       "Hill Climbing (HC) heuristic for DC cost (only) by using BB and/or ODT-naive \n"

       " --HC: run hill climbing heuristic for each initial network using cost function and print optimal cost, all optimal networks are written in odt.log file; summary stats are saved to odt.dat;\n"
       "  By default each HC step starts from the set of predefined networks (if -n or -N are provided); otherwise by using quasi-random networks\n"              
       " --hcnnimove: use NNI instead of TailMoves\n"      
       " --hcrunstats: print improvements stats after each HC runs\n"
       " --hcrunstatsext: print improvements stats after locating optimal networks after each HC runs\n"
       " --hcrunstatsalways: print stats after each HC run\n"
       " --noodtfiles: do not save odt.log and odt.dat with optimal networks\n"
       " --hcusenaive: use only ODT naive algorithm in HCH (no BB)\n"
       " --hcstopinit=NUM: stop when there is no new optimal network after NUM HC runs\n"
       " --hcstopclimb=NUM: stop in a single HC when there is no improvements after NUM steps\n"
       " --hcmaximprovements=NUM: stop after NUM improvements in HC climb\n"       
       " -O ODTFILE: the name of odt.log \n"
       " -D DATFILE: the name of odt.dat; see also --odtlabelled and --noodtfiles \n"
       " --noodtfiles: do not generate odt and dat files\n"
       " --odtlabelled: odt and dat file in labelled format\n"

       " By default HC is limited to tree-child and non time consistent networks. Use --general, --relaxed or/and --timeconsistent to change the search space.\n"

       "\n"

       "\nVerbose:"
       " --verbose, -v [123][34]: where\n"
       "  0: quiet mode in HC, BB and ODT-naive\n"
       "  1: print visited network in HC after each improvement\n"
       "  2: print visited networks in HC if the cost is equal to the current or improved\n"
       "  3: print all visited networks in HC \n"
       "  4: print basic info in HC, BB, ODT-naive \n"
       "  5: print detailed info in HC, BB, ODT-naive \n"

       "\n"

       "Networks/dags comparing and aggregating\n"
       " --comparedags: pairwise comparison of all networks or their shapes\n"
       " --uniquedags: print unique networks; summary stats printed on stderr\n"
       " --uniquedagscnts: as above plus counts of networks\n"
       " --dagshapes: two networks are equal if their shapes are isomorphic in --*dags* options (i.e., ignore leaf labels)\n"

       "\n"

       "Other options, print, debug, test, etc.\n"
       " --pspeciesdictionary: list species dictionary\n"      
       " --pdetailed: print debug tree and network structures\n"
       " --reachablenodescnt: for each v in V(N), print the number of nodes reachable from v\n"
       " --reachableleafvescnt: for each v in V(N), print the number of leaves reachable from v\n"         
       " --pstsubtrees: print species tree subtrees\n"
       " --maxdisplaytreecachesize: set limit for the size of the cache of display tree nodes (only if DTCACHE is enabled)\n"

                       
       "\n"
       // "*Starting tree options:\n"
       // "  -q NUM - generate NUM quasi-consensus trees sampled from gene tree clusters;\n"
       // "       (def. is 1 if no user species tree is defined)\n"
       // "  -r NUM - generate NUM random trees (def. 0)\n"
       // "  -Q p - defines probability for sampling gene trees (use with -q, default 0.2)\n"
       // "  See also -qNUM, -ei and -er.\n"

       "Examples: \n\n"

       "Print 10 quasi consensus trees\n"
       "  supnet -g '(a,((b,c),d));(a,(b,d))' -q10 --pnetworks\n"
       
       "Print 10 quasi consensus trees with preserved split of the root (-er)\n"
       "  supnet -g '(a,((b,c),d));(a,(b,e))' -s'((a,b),(c,(d,e)))' -q10 --pnetworks --preserveroot\n"

       "Print DC cost\n"
       "  supnet -g '(b,(a,c))' -s '(a,(b,c))' -CDC --ptreecost\n"

       "Print 10 random species trees (i.e., networks with no reticulation) over a..e\n"
       "  supnet -A5 -r10 --pnetworks\n"

       "Print display trees based on reticulation switching; trees maybe non-unique\n"
       "  supnet -n '((((c)#B,b))#A,(#A,(#B,a)))' --pdisplaytrees\n"

       "Print the cost of 10 random gene trees vs random tree-child network with 5 reticulations over 8 species\n"
       "  supnet -r10 -A8 --pnetworks | supnet -G- -r1 -A8 -R5 --odtnaivecost --pnetworks\n"

       "Printing and visualizing\n"                     
       "  supnet -n '((a)#A,(#A,(c,b)))' --dot | dot -Tpdf > n.pdf\n"
       "  supnet -n '((((c)#B,b))#A,(#A,(#B,a)))' --dot | dot -Tpdf > n.pdf\n"
       "  supnet -r1 -R10  -A20 --dot | dot -Tpdf > n.pdf\n"
       //"  N=$( embnet.py -n 'rand:3:1' -pn ); echo $N; supnet -n$N -d | dot -Tpdf > n.pdf\n"

       "\nDynamic programming --DP (approximate DC cost)\n"
       "  supnet -g '(a,(b,c))' -n '((a)#A,(#A,(c,b)))' --DP\n"

       "\nBranch and bound --BB (exact DCE cost)\n"
       "  supnet -g '(a,(b,c))' -n '((a)#A,(#A,(c,b)))' --BB\n"

       "\nHill climbing heuristic (--HC ...)\n"

       " Minimalistic run: print cost; result in odt.log; one random initial network (-r1, i.e., one HC climb) with 2 (-R2) reticulations\n"
       "  supnet -g '(a,(b,(c,d))); ((a,b),(c,d))' -r1 -R2 --HC\n"
       
       " Print summary stats; 10 HC runs with random networks\n"
       "  supnet -g '(a,(b,(c,(d,e)))); ((a,b),(c,(e,a))); ((b,c),(d,a))' -r10 -R3 --HC -v2 --hcpstatsext\n"

       " Print cost, improvements and stats; tree-child search; quasi consensus initial network (-q) with 3 random reticulations\n"
       "  supnet -g '(a,(b,(c,(d,e)))); ((a,b),(c,(e,a))); ((b,c),(d,a))' -q1 -R3 --HC -v2 --hcpstats\n"

       " Larger instance; relaxed and time consistent network search:\n"
       "  supnet -g '((i,c),((a,d),(f,(b,((g,(e,j)),h)))));((i,h),((c,f),((a,(d,e)),(j,(g,b)))));(((f,(b,d)),((j,g),(e,i))),(c,(h,a)));(((f,(i,(j,d))),(c,b)),(((h,a),g),e));((((h,e),((c,f),a)),(d,b)),((g,i),j))' -R8 -q1 --HC -v2 --hcpstats --relaxed --timeconsistent\n"

       " Print only improvements, general network search:\n"
       "  supnet -g '(a,(b,(c,d))); ((a,b),(c,d))' -r1 -R3 --HC -v1 --general \n"
       
       " Print only improvements and equal cost networks:\n"
       "  supnet -g '(a,(b,(c,d))); ((a,b),(c,d))' -r1 -R3 --HC -v2\n"
       
       " Print improvements; skip odt.log and odt.dat:\n"
       "  supnet -g '(a,(b,(c,d))); ((a,b),(c,d))' -r1 -R3 -v1 --HC --noodtfiles \n"

       " Recommended with large HC-runs using quasi-consensus rand networks\n"
       "  supnet -g '(a,(b,(c,d))); ((a,b),(c,d))' -q1000 -R3 --HC --hcpstatsext \n"

       " Using stopping criterion (1000) with quasi-consensus networks\n"
       "  supnet -g '(a,(b,(c,d))); ((a,b),(c,d))' -q-1 -R3 --HC --hcpstatsext --hcstopinit=1000\n"

       " Display trees usage stats:\n"
       "  supnet -N odt.log --pdisplaytrees | sort | uniq -c | sort -k1 -n\n"

       " Insert 2 reticulations into a network (tree-child output)\n"
       "  supnet -R2 -n '(a,((d)#1,(b,(c,#1))))' --pnetworks \n"

       " Insert 10 reticulations into a general network\n"
       "  supnet -R10 -n '(a,((d)#1,(b,(c,#1))))' --pnetworks --general\n"

       " Pairwise comparison of random dags (general networks)\n"
       "  supnet -r10 -R1 -A2 --pnetworks --general --comparedags\n"

       " Print unique networks from odt.log\n"
       "  supnet -g '(a,(b,(c,d))); ((a,b),(c,d))' -r1 -R3 --HC -v2 && supnet --uniquedags -N odt.log\n"

       " Print unique random networks with counts\n"
       "  supnet -r100000 -R1 -A3 --uniquedagscnts\n"

       " Print unique random networks with counts under uniform pairs model\n"
       "  supnet -r100000 -R1 -A3 --uniquedagscnts --uniformedgesampling\n"

       " Print unique random shapes of networks with counts\n"
       "  supnet -r100000 -R1 -A3 --dagshapes --uniquedagscnts \n"

       " Gen BB-tree search\n" 
       "  supnet -r1 -A15 --pnetworks  | supnet -G- -r1 -A15 -R12 --BB --bbtreesearch; dot bb.dot -Tpdf > bb.pdf\n"

       ;

/*
 To test correctness of the above commands run:
grep '^[ \t]*"  supnet' tools.cpp | sed 's#[\]"#\"#g;  s/^[[:space:]]*"/( /g; s/\\n"/ ) || exit -1; echo ==================;sleep 1;/g' | bash -x 
*/

  exit(-1);
}


// return random tree string from given vector of species
string genrandomtree(NODEID *sp, int len)
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


void shuffle(NODEID *a, size_t n)
{
    if (n > 1) 
    {
        size_t i;
        for (i = 0; i < n - 1; i++) 
        {
          size_t j = i + rand() / (RAND_MAX / (n - i) + 1);
          NODEID t = a[j];
          a[j] = a[i];
          a[i] = t;
        }
    }
}


//  repr:
// 0 - lf
// 1 - rt
// rest: tree rt+1 encodings

ostream& pprepr(SPID *t, ostream&s) 
{
  int lf = t[0];
  int rt = t[1];
  // cout << "[[[";
  // for (int i=2; i<COMPRSIZE(lf,rt); i++)
  //   {
  //     cout << (int)t[i]  << ".";
  //   }
  // cout << "]]]" << endl;
  int lst=2;
  for (int i=0; i<=rt; i++)
    lst += ppreprsingle(t+lst, s);
  
  return s;
}


int ppreprsingle(SPID *t, ostream&s) 
{
    SPID prev = REPROPEN;
    int i = 0;
    int open = 0;
    int speciescnt = specnames.size(); 

    while(1)
    {
        SPID c = t[i];        
        if (c==REPROPEN) 
        { 
          open++;
          if (prev!=REPROPEN && prev!=MAXSPECIES) s << ",";
          s << "("; 
        }
        else if (c==REPRCLOSE) { 
          s << ")";       
          open--;
        }
        else { 
          if (prev!=REPROPEN) s << ",";
          if (c<specnames.size())
            s << specnames[c];                    
          else
            s << "^" << specnames[c-speciescnt];
        }
        prev = c;
        i++;
        if (!open) break;
    }    
    return i;
}


// Adopted from https://stackoverflow.com/questions/33010010/how-to-generate-random-64-bit-unsigned-integer-in-c
#define IMAX_BITS(m) ((m)/((m)%255+1) / 255%255*8 + 7-86/((m)%255+12))

#define RAND_MAX_WIDTH IMAX_BITS(RAND_MAX)

uint64_t rand64(void) {
  uint64_t r = 0;
  for (int i = 0; i < 64; i += RAND_MAX_WIDTH) {
    r <<= RAND_MAX_WIDTH;
    r ^= (unsigned) rand();
  }
  return r;
}



// Return the size of total used memory in MB
unsigned long get_memory_size()
{
  unsigned long dummy, size;

  const char* statm_path = "/proc/self/statm";

  FILE *f = fopen(statm_path,"r");
  if(!f){
    perror(statm_path);
    abort();
  }
  if(7 != fscanf(f,"%ld %ld %ld %ld %ld %ld %ld",
    &size,&dummy,&dummy,&dummy,&dummy,&dummy,&dummy))
  {
    perror(statm_path);
    abort();
  }
  fclose(f);

  return size * getpagesize() / (1024 * 1024 );

  //                 size       (1) total program size
  //                            (same as VmSize in /proc/[pid]/status)
  //                 resident   (2) resident set size
  //                            (inaccurate; same as VmRSS in /proc/[pid]/status)
  //                 shared     (3) number of resident shared pages
  //                            (i.e., backed by a file)
  //                            (inaccurate; same as RssFile+RssShmem in
  //                            /proc/[pid]/status)
  //                 text       (4) text (code)
  //                 lib        (5) library (unused since Linux 2.6; always 0)
  //                 data       (6) data + stack
  //                 dt         (7) dirty pages (unused since Linux 2.6; always 0)
}





