/*
 * tools.h
 *
 *  Created on: 18-09-2010
 *      Author: gorecki
 */


#ifndef TOOLS_H_
#define TOOLS_H_

#include <fstream>
using namespace std;

typedef long COSTT;

unsigned long get_memory_size();

#define INFTY 1000000

#ifndef MAXSPECIES
  #define MAXSPECIES 64    
  #define MAXSPID 255
#endif

#if MAXSPECIES<=253
  #define SPID unsigned char
#elif 
  #define SPID unsigned short int
#endif

#define REPRCLOSE (MAXSPID-2)
#define REPROPEN (MAXSPID-1)

#if !defined(SPSMALL) && !defined(SPMED) 
#define SPMED
#endif

#ifdef SPSMALL
#define NODEID char
#define MAXNODEID 255
#endif

#ifdef SPMED
#define NODEID short
#define MAXNODEID 32000
#endif

typedef NODEID *NODEIDARR;

typedef NODEID NODEIDPair[2];

#define GSFULL 0
#define GSPOS 1
#define GSAFTER 2
#define GSBEFORE 3


#include <set>
using namespace std;

#include <vector>
#include <iomanip>
#include <map>
using namespace std;

#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <random>

#define ERR_NOSPECIESDEFINED "No species are defined."

class GTCluster;

extern vector<string> specnames;
extern vector<int> specorder;
extern map<string, int> specnames2id;
extern vector<GTCluster*> spec2gtcluster;

extern int rsort; // Important in nni - used to generated unique r|u species trees
extern int gcenter; // Useful for drawing gene trees
extern int costfromroot;
extern int printleaves;
extern int stprogress; // Print progress for species trees
extern int ppgraphviz; // Print graphviz output
extern NODEID topspcluster[2];
string getstringn(char *s, int len);
extern int flag_globaldagcache; 

// Count occurences of c in s
int strcount(char *s, char c);

#define species(i) specnames[i]

#define tostr(x) dynamic_cast< std::ostringstream & >(( std::ostringstream() << std::dec << x ) ).str()

string genrandomtree(NODEID *sp, int len);

NODEID* joinspclusters(NODEID* a, NODEID* b, NODEID *dest=NULL, bool compresstopcluster=true);
NODEID* spidcopy(NODEID* a, int size);
void printspcluster(ostream&s, NODEID *a);
void pprintspcluster(ostream&s, NODEID *a);

void initlinenuminfo(char *);
void printlinepos();
void checkparsing(char *);
void expectTok(const char *, char *, int );

#define spsize(x) ((((x)[0])==0)?specnames.size():((x)[0]))

char *mstrndup(const char *s, size_t n);

char* getTok(char *s, int &p, int num);
char* seeTok(char *s, int p, int num);
int getspecies(char *s, int len);
int getspecies(string s);
void setspecies(int spcnt);

void randomizearr(NODEID*, int);

int usage(int argc, char **argv);

int eqspclusters(NODEID *a, NODEID *b);
bool spsubseteq(NODEID *a, NODEID *b);

void cleanspecies();

double gettime();

int ppreprsingle(SPID *t, ostream&s = cout); 
ostream& pprepr(SPID *t, ostream&s = cout);

struct comparespids
{
  bool operator()(const NODEID* a, const NODEID* b) const
  {  
      int i = 1, al = a[0] + 1, bl = b[0] + 1;

      if (al<bl) { return true; }
      if (al>bl) { return false; }
      for ( ; i < al; i++ ) {
      if (a[i] < b[i]) { 
        return true;
      }
      if (a[i] > b[i]) { 
        return false; }
      }      
      return false;  
  }
};


void shuffle(NODEID *a, size_t n);



uint64_t rand64(void);

#define COMPRSIZE(LF,RT) (3*(LF)+4*(RT))

#endif /* TOOLS_H_ */
