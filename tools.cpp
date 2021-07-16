

/************************************************************************
 SuperNetwork Inference - based on FastUrec project

 *************************************************************************/

#include <ctype.h>
#include <iostream>
#include <cstdlib>
#include <cstring>

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

       "  -l [aDELIM|bDELIM|pPOS] - rules for matching gene tree names with species tree names;\n"       
       "       a - after delimiter; b - before delimiter; p - from position \n"
       
       "\n"
       "ADDITIONAL OPTIONS\n"
       "  -v NUM - verbose level 1 (print the name of outputfile only),\n"
       "       2 (simple output), 3 (detailed output)\n"
       "  -w NUM - print NUM quasi-consensus species trees and exit\n"              
       "  -r NUM - print NUM random species trees and exit\n"
       "  -S SPECIESNUM - define SPECIESNUM species a,b,...\n"

       "  -e [gsrD...]+ - multiple options:\n"
       "     g - print a gene tree\n"
       "     s - print a species tree\n"
       "     n - print a network\n"
       "     d - print display trees\n"
       "     o - ODT cost using naive approach\n"
       "     r - preserve root when searching the species tree space and in quasi-consensus\n"              
       "     D - detailed tree info\n"       
       "     i - list species dictionary\n"       
       "     c - print cost\n"       
       "     C - print cost with trees\n"       

       "\n"       
       "COST SETTING OPTIONS\n"
       // " TODO: -D dupweight  - set weight of gene duplications (def. 1.0)\n"
       // " TODO: -L lossweight - set weight of gene losses (def. 1.0)\n" 
       "  -C COST - set cost function from {DL,D,L,DC,RF}\n"
       "\n"

       "HEURISTIC SEARCH\n"
              
       // "*Starting tree options:\n"
       // "  -q NUM - generate NUM quasi-consensus trees sampled from gene tree clusters;\n"
       // "       (def. is 1 if no user species tree is defined)\n"
       // "  -r NUM - generate NUM random trees (def. 0)\n"
       // "  -Q p - defines probability for sampling gene trees (use with -q, default 0.2)\n"
       // "  See also -wNUM, -ei and -er.\n"

       "Examples: \n\n"

       "Print 10 quasi consensus trees\n"
       " supnet -g '(a,((b,c),d));(a,(b,d))' -w10\n"
       
       "Print 10 quasi consensus trees with preserved split of the root (-er)\n"
       " supnet -g '(a,((b,c),d));(a,(b,e))' -s'((a,b),(c,(d,e)))' -w10 -er\n"

       "Detailed tree info\n"
       "  supnet -g '(a,((b,a),c))' -eD\n"

       "Print DC cost\n"
       "  supnet -g '(b,(a,c))' -s '(a,(b,c))' -CDC -ec \n"

       "Print 10 random species trees over a..e (5)\n"
       "  supnet -S5 -r10\n"

       "Print display trees\n"
       "  supnet -n '((((c)#B,b))#A,(#A,(#B,a)))' -ed\n"

       "Printing and visualizing\n"
       "  supnet -n '((a)#A,(#A,(c,b)))' -enD\n"       
       "  supnet -n '((((c)#B,b))#A,(#A,(#B,a)))' -enD\n"       
       "  supnet -n '(((b)#A,a),(#A,c))' -d\n"
       "  supnet -n '((a)#A,(#A,(c,b)))' -d | dot -Tpdf > n.pdf\n"
       "  supnet -n '((((c)#B,b))#A,(#A,(#B,a)))' -d | dot -Tpdf > n.pdf\n"
       "  supnet -n $( embnet.py -n 'rand:20:10' -pn ) -d | dot -Tpdf > n.pdf\n"
       "N=$( embnet.py -n 'rand:3:1' -pn ); echo $N; supnet -n$N -d | dot -Tpdf > n.pdf\n"

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




