/*
 * tools.h
 *
 *  Created on: 18-09-2010
 *      Author: gorecki
 */

#ifndef TOOLS_H_
#define TOOLS_H_

#define SPMED
//#define SPSMALL

#ifdef SPSMALL
#define SPID char
#define MAXSP 126
#endif

#ifdef SPMED
#define SPID short
#define MAXSP 32760
#endif

#ifdef SPLARGE
#define SPID int
#define MAXSP 1000000
#endif


typedef SPID *SPIDARR;

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

extern vector<string> specnames;
extern vector<int> specorder;
extern map<string, int> specnames2id;
extern map<SPID, SPID*> spec2spcluster;
extern int rsort; // Important in nni - used to generated unique r|u species trees
extern int gcenter; // Useful for drawing gene trees
extern int costfromroot;
extern int printleaves;
extern int stprogress; // Print progress for species trees
extern int ppgraphviz; // Print graphviz output
extern SPID topspcluster[2];
string getstringn(char *s, int len);

// Count occurences of c in s
int strcount(char *s, char c);

#define species(i) specnames[i]

#define tostr(x) dynamic_cast< std::ostringstream & >(( std::ostringstream() << std::dec << x ) ).str()

string genrandomtree(SPID *sp, int len);

SPID* joinspclusters(SPID* a, SPID* b, SPID *dest=NULL);
SPID* spidcopy(SPID* a, int size);
void printspcluster(ostream&s, SPID *a);

void initlinenuminfo(char *);
void printlinepos();
void checkparsing(char *);
void expectTok(const char *, char *, int );

#define spsize(x) ((((x)[0])==0)?specnames.size():((x)[0]))

char *mstrndup(const char *s, size_t n);

char* getTok(char *s, int &p, int num);
int getspecies(char *s, int len);

void randomizearr(SPID*, int);

int usage(int argc, char **argv);

int eqspclusters(SPID *a, SPID *b);
bool spsubseteq(SPID *a, SPID *b);

void cleanspecies();

struct comparespids
{
  bool operator()(const SPID* a, const SPID* b) const
{
      int i = 1, al = a[0] + 1, bl = b[0] + 1;

      if (al<bl) { //cout << "Tx";
  return true; }
      if (al>bl) { //cout << "Fx";
  return false; }
    for ( ; i < al; i++ ) {
      if (a[i] < b[i]) { //cout << "T";
  return true;}
      if (a[i] > b[i]) { //cout << "F";
  return false; }
    }
      //    cout <<"F";
    return false;
  
  }
};





#endif /* TOOLS_H_ */
