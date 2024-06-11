# SupNet

Phylogenetic networks visualisation, evaluation and inference\n"

## Usage

```
supnet ...
```

## Basic input and output 

### Reading trees and networks

 - `--genetrees, -g STR`: gene trees from a string; separator is semicolon
 - `--speciestrees, -s STR`: species trees from a string
 - `--networks, -n STR`: networks from a string
 - `--genetreesfile, -G FILE`: gene trees from a file; EOLN separated; - stdin
 - `--speciestreesfile, -S FILE`: species trees from a file; EOLN separated; - stdin
 - `--networksfile, -N filename`: networks from a file; EOLN separated; - stdin

### Printing trees and networks

 - `--pgenetrees`: print gene trees
 - `--pspeciestrees`: print species trees
 - `--pnetworks`: print networks
 - `--pdisplaytrees`: print all display trees
 - `--pdisplaytreesext`: print all display trees with their ids
 - `--pranddisplaytrees=NUM`: print NUM randomly chosen display trees
 - `--dot`: print dot representation of networks

### Gene tree(s)

Read (`-g`) and print gene trees (`--pgenetrees`).

```
> supnet -g '(a,(b,c));(a,(a,d))' --pgenetrees
(a,(b,c))
(a,(a,d))
```

Gene trees from a file `-G`. 

```
> supnet -G examples/gtrees.txt --pgenetrees
(a,(b,c))
(a,(a,d))
```

Use `-` for stdin.

```
> cat examples/gtrees.txt | supnet -G- --pgenetrees
(a,(b,c))
(a,(a,d))
```

### Species tree(s) 

Use `-s/-S` and `--pspeciestrees` to read and print species tree sets.

### Networks

Use `-n/-N` and `--pnetworks` to read and print networks.

### Display trees of networks

Print display trees based on reticulation switching `--pdisplaytrees`. Trees maybe non-unique.

```
>  supnet -n '((((c)#B,b))#A,(#A,(#B,a)))' --pdisplaytrees
((c,a),b)
((c,b),a)
((c,a),b)
((c,b),a)
```

Print all display trees with their ids '--pdisplaytreesext'.
```
>  supnet -n '((((c)#B,b))#A,(#A,(#B,a)))' --pdisplaytreesext
0 ((c,a),b)
1 ((c,b),a)
2 ((c,a),b)
3 ((c,b),a)
```

Print 2 random display trees with ids '--pdisplaytreesext' (optional).

```
> supnet -n '((((c)#B,b))#A,(#A,(#B,a)))' --pranddisplaytrees=2 --pdisplaytreesext
1 ((c,b),a)
1 ((c,b),a)
```

Display tree usage stats after HC run:
```  
> supnet --HC -g '(a,(b,c));(b,(a,c))' -R1 -q1; supnet -N odt.log --pdisplaytrees | sort | uniq -c | sort -k1 -n
   i: (a,(((c)#1,b),#1)) cost=1
   >: ((a,#1),((c)#1,b)) cost=0
Cost:0 TopNetworks:1 Steps:18 Climbs:4 HCruns:1 HCTime:4.22001e-05 Naive:19 Naivetime:9.05991e-06 DTcnt:38 Class:TreeChild TimeConsistency:0
      1 ((a,c),b)
      1 (a,(c,b))
```

### Visualization of networks using dot format and graphviz

To use dot visualizations install `graphviz` package.

```
supnet -n '(((b)#A,a),(#A,c))' --dot
```

```
>  supnet -n '((a)#A,(#A,(c,b)))' --dot | dot -Tpdf > n.pdf
```

```
>  supnet -n '((((c)#B,b))#A,(#A,(#B,a)))' --dot | dot -Tpdf > n.pdf
```

```
>  supnet  -r1 -R10  -A20 --dot | dot -Tpdf > n.pdf
```

## (Not implemented yet) Matching gene tree names with species names

 - `[TODO] --matchparam, -l [NUM|STR]`: parameter for the occurrence of species name in gene label; see `--match*` options below
 - `[TODO] --matchafter`: species name is after `--matchparam (STR)` in gene label
 - `[TODO] --matchbefore`: species name is before `--matchparam (STR)` in gene label
 - `[TODO] --matchatpos`: species name is at fixed `--matchparam (NUM)` position in gene label; negative `NUM` positions from the end of label

## Network classes

Binary network classes in random generator and HC heuristic: 
 - `--treechild`: tree-child networks (default); a non-leaf node has at least one non-reticulation child 
 - `--relaxed`: a node has at most one reticulation child; see -R
 - `--general`: general phylogenetic networks (not applicable with DP algorithm) 
 - `--timeconsistent`: time-consistent networks 
 - `--notimeconsistent`: no time-consistent networks 
 - `--detectclass`: for each network N print "t c r N", where t is 1 if N is time consistent, c is 1 if N is tree-child network, r is 1 if N is relaxed network. If network N does not belong to a given class the corresponding value is 0.

These values are reported in `[Cc]lass` and `TimeConsistency` field in output and odt.dat file. See HC heuristic.

Detecting classes:
```
> supnet -r10 -R2 -A5 --general --detectclass  --randseed 13
0 0 0 ((c,((d)#2,((b)#1,(a,(#2,#1))))),e)
0 1 1 ((((((c,((d)#1,b)),#1))#2,a),e),#2)
0 1 1 (((#1,((e)#1,a)),b),((c)#2,(d,#2)))
0 0 0 (((a)#1,#2),((((d,#1))#2,c),(b,e)))
0 1 1 ((c)#1,(#1,((e)#2,(((#2,d),a),b))))
1 1 1 ((c,#1),((#2,((a)#1,b)),((e)#2,d)))
0 1 1 (((#1,(((c)#2,a),#2)),d),((b)#1,e))
0 0 0 ((((e)#2,(#2)#1),((d,c),(b,#1))),a)
0 1 1 (#1,((((d,((c)#2,a)))#1,(e,b)),#2))
0 1 1 (((((#2,((b)#2,a)),d))#1,(c,#1)),e)
```

## Random trees and networks and quasi consensus generators

Random and quasi consensus tree/network generators
 - `-r NUM`: generate `NUM` random networks; if `NUM=-1` generator is infinite
 - `-q NUM`: generate `NUM` quasi-consensus networks; if `NUM=-1` generator is infinite
 - `-R NUM`: insert NUM reticulation nodes into networks from `-q`, `-r`, `-n`, `-N`; default is `0`
 - `-A SPECIESNUM`: define `SPECIESNUM` species a,b,... [TODO: expand label range]
 - `--uniformedgesampling`: draw uniformly a pair of edges to create random reticulation; the default: first draw the source edge, then the uniformly the - destination edge given the source
 - `--preserveroot`: preserve root in quasi-consensus and in HC algorithm
 - `--randseed, -z SEED`: set seed for random generator (srand)

If NUM is -1, the network generator will provide network when requested (e.g., when used with `-K/--hcstop`, see `--HC` option).

By default, the generator produces tree-child networks without time-consistency restriction. 

Print 2 quasi consensus trees with preserved split of the root.
```
> supnet -g '(a,((b,c),d));(a,(b,e))' -s'((a,b),(c,(d,e)))' -q2 --preserveroot --pnetworks
((d,(e,c)),(a,b))
((b,a),(d,(e,c)))
```

### Guide clusters: used in quasi consensus 

Guide clusters are defined as a list of multifurcated trees separated by `;`. All clusters from the trees must be present in all quasi consensus tree and networks in hill climbing. Note that random networks are not generated using guide clusters (i.e., with `-r NUM`).

```
> supnet --guideclusters "(a,b,c);(d,e,(f,g))" -q4 --pnetworks
(((b,c),a),((d,e),(g,f)))
(((a,c),b),((d,e),(f,g)))
(((g,f),(d,e)),((a,b),c))
(((e,d),(g,f)),(c,(a,b)))
```

### Guide tree: forcing tree structure in networks

Similarly to guide clusters guide trees are defined as a multifurcated tree. 
A cluster present in every guide tree must be present as a cluster of a node in a network. Additionally, every such a node is a root of a subtree, i.e., no reticulation is allowed.
Note that random networks are not not generated using guide trees (i.e., with `-r NUM`).

```
> supnet --guidetree "((a,b,c),(d,e,f))" -q4 --pnetworks
(((d,f),e),((a,c),b))
((d,(f,e)),(a,(b,c)))
(((f,d),e),((a,b),c))
(((e,d),f),(c,(a,b)))
```

With reticulations; guide trees are separated by `;`.
```
> supnet -A10 -q2  --pnetworks --guidetree '(a,b);(d,e,f)' -R3
((((#3,(((h)#2,g),c)))#1,(b,a)),((((#2,((f,d),e)))#3,j),(i,#1)))
((h,(g,(((((b,a))#2,j))#1,((c)#3,(#2,(e,(f,d))))))),((i,#1),#3))
```

### Random trees
    
Print two random networks with no reticulations (i.e., species trees) using species {a,b,c,d,e} (`-A5`). 
```
>  supnet -A5 -r2 --pnetworks
((c,((e,a),b)),d)
(((a,d),b),(c,e))
```

Print two quasi-consensus trees.
```
> supnet -g '(a,((b,c),d));(a,(b,d))' -q2 --pnetworks
((c,a),(b,d))
((b,d),(a,c))
```

### Random networks 

Print one random relaxed network with species {a,b,c}  and 5 reticulations.
```
>  supnet -A3 -R5 -r1 --pnetworks --relaxed
(((#1)#5,(#4,(#5,a))),(#2,((b)#3,((((#3)#1)#2)#4,c))))
```

Print one quasi-consensus network with two reticulations.
```
> supnet -G examples/gtrees.txt -q1 -R2 --pnetworks --pnetworkclusters
(((a)#2,((#1,d),((c)#1,b))),#2)
a 
b 
c 
d 
b c 
c d 
b c d 
a b c d 
```

Print one quasi-consensus network with one reticulations and having guide tree clusters, and print all network clusters.
```
> supnet -q1 -R2 --pnetworks --guideclusters '((a,b,c),(d,e,f))' --pnetworkclusters
(((f)#1,(#2,((b)#2,(c,a)))),(#1,(d,e)))
a 
b 
c 
d 
e 
f 
a c 
d e 
a b c 
d e f 
a b c f 
a b c d e f 
```


### Fixed seed for random generator

Use `-z SEED` or `--randseed SEED` to set seed for random generator (srand). 

### Inserting reticulations

Insert 2 reticulations into a network; tree-child network in output (default).
```  
> supnet -R2 -n '(a,((d)#1,(b,(c,#1))))' --pnetworks
(#2,((d)#1,(#3,((b)#3,((a)#2,(c,#1))))))
```

Insert 8 reticulations into a network; general network `--general`
```  
> supnet -R8 -n '(a,((d)#1,(b,(c,#1))))' --pnetworks --general 
((((#3)#2)#9,((d)#5,((#1)#8,#2))),(((#9)#7,(#5)#1),((((((#8)#4,#7),((a)#3,c)),#4))#6,(b,#6))))
```

## Cost functions 

Available cost function are {DL,D,L,DC,RF,DCE}; the default is DC. 
Limitations: 
- In DP only DC is available, HC cannot be run using PD via BB.
- DTCache has only DC cost implementation.

Cost functions:
 - `--cost, -C COST`: set cost function from {DL,D,L,DC,RF,DCE}; default is DC
 - `[TODO] -D dupweight`: set weight of gene duplications (def. 1.0)
 - `[TODO] -L lossweight`: set weight of gene losses (def. 1.0)

## Tree(s) vs tree(s) cost computation:

- `--ptreecost`: print cost between two trees (G,S)
- `--ptreecostext`: print cost between two trees with trees (G,S)

Print DC cost
```
>  supnet -g '(b,(a,c))' -s '(a,(b,c))' -CDC --ptreecost
1
```

Print DCE cost (DCE = DC + |V(E)|).
```
>  supnet -g '(b,(a,c))' -s '(a,(b,c))' -CDCE --ptreecost
5
```

Print RF cost [TODO]
```
>  supnet -g '(b,(a,c))' -s '(a,(b,c))' -CRF --ptreecost
RF is not implemented yet
```

## Optimal display tree cost: network(s) vs collections of gene tree: 

For DC cost:
```
> supnet -n '((c,#1),((#2,((a)#1,b)),((e)#2,d)))'  -g '(e,((a,b),(d,c)));(d,((c,a),(e,b)));((d,b),((a,c),e))'  --odtcost
6 ((c,#1),((#2,((a)#1,b)),((e)#2,d)))
```

## DP algorithm to compute lower bound of DC between a gene tree and a network by dynamic programming

`--DP`: run DP and print the bound


DP algorithm computes aproximately the DC cost. Use `--DP` to execute it directly.

```
> supnet -g '(a,(b,(c,d)))'  -n '((#1,d),(((c)#1,b),a))' --DP
7
```


## ODT-naive algorithm via enumeration of display trees to compute exact DC cost between gene trees and a network: 

 - `--odtnaivecost`: run ODT-naive and print the cost
 - `[TODO] --odtnaivesampling`: set sampling in ODT naive (exponential distribution)
 `--ptreesinodtnaive`: print additional details


## Branch and bound (BB) algorithm 

BB algorithm computes exactly DCE cost. It may call DP and ODT naive multiple times. Use `--BB` to execute it directly.

 - `--BB`: compute exact DCE by BB
 - `-t THRESHOLD, --runnaiveleqrt`: run ODT-naive cost computation, when the number of reticulations is `<= THRESHOLD`; otherwise run DP;
 - `--bbtsvstats`: gen bb.tsv with stats
 - `--bbtreesearch`: gen bb.dot with bb tree search
 - `--bbtimestats`: print time per each pair
 - `--bbstartscore=FLOAT`: define initial score in BB (testing only)


Generate the picture of BB-tree search:
```  
> supnet -r1 -A15 --pnetworks  | supnet -G- -r1 -A15 -R12 --BB --bbtreesearch; dot bb.dot -Tpdf > bb.pdf
-4
```

## Hill Climbing (HC) heuristic for DC cost (only) by using BB and/or ODT-naive 

 - `--HC`: run hill climbing heuristic for each initial network using cost function and print optimal cost, all optimal networks are written in odt.log file (default) or in odt.tre (if --saveextnewick is set); summary stats are saved to odt.dat;
 - ` --hcnnimove`: use NNI instead of TailMoves
 - ` --hcpstats`: print stats after locating optimal networks after each HC runs
 - ` --hcpstatsext`: print extended stats after each HC run

 - `--hcrunstats`: print improvements stats after each HC run
 - `--hcrunstatsext`: print improvements stats after locating optimal networks after each HC run
 - `--hcrunstatsalways`: print stats after each HC run
 - `--noodtfiles`: do not save odt.log/odt.tre and odt.dat with optimal networks
 - `--hcusenaive`: use only ODT naive algorithm in HCH (no BB)
 - `--hcstopinit=NUM`: stop when there is no new optimal network after HC with NUM initial networks
 - `--hcstopclimb=NUM`: stop in HC when there is no improvements after NUM steps
 - `--hcmaximprovements=NUM`: stop after NUM improvements in HC climb
 - `-O BASENAMEFILE, --outfiles`: base name of odt output files .log/.tre and .dat
 - `-d DIRECTORY, --outdirectory`: directory to save output files, default is the current dir
  - `--odtlabelled`: odt and dat file in labelled format

 By default HC is limited to tree-child and non time consistent networks. Use `--general`, `--relaxed` or/and `--timeconsistent` to change the search space.

 ODT is solved by a hill climbing heuristic algorithm that optimizes a cost function for networks. The algorithm will start from a set of predefined networks or generate random networks partially based on input gene trees if specified. The algorithm will make use of two moves: "tail" and "nni" to improve the networks.

The algorithm will keep track of various statistics. It will print the optimal cost found and write all optimal networks to an "odt.log" file. Summary statistics, including the optimal cost, total time taken, time taken for hill climbing, time taken for the merge step, number of networks, number of improvements made, number of steps taken, and the number of starting networks will be saved to an "odt.dat" file.

Sketch of the algorithm:

```
Input: a set of gene trees GSet and a cost function (now only DC).

For each initial (random, quasi or input) network N:
  Run hill climbing HC(GSet, N) and store the result if the resuling cost is not worst than the current one.
  If -K NUM is set, stop if there was no better network in the last NUM steps.
Return the best scoring set of networks

HC(GSet, N): 
  cost := ODTcost(GSet,N)
  while there is a neighbour N' of N:
     cost' := ODTcost(GSet, N')
     if cost' is better than cost: 
         clean stored networks
         cost := cost'
         N := N'
         continue
      if cost' = cost: 
         store N'
   if there is a stored and not visited network N':
      N := N'
      go to while loop 
   return N, cost

ODTcost(GSet, N): 
    if `--hcusenaive` is set or N has less than '-t THRESHOLD` reticulations:
       compute the cost by ODT naive algorithm, i.e., enumeration of all possible display trees)
    otherwise:        
       compute the cost by BB algorithm (using multiple calls of DP and/or ODT naive)
```

By default each HC step starts from the set of predefined networks (see `-n`, `-N` or random networks); otherwise by using quasi-random networks.

Minimalistic hill climbing (HC) run using Tail Moves (default) with a single random network (and one HC run) `-r1` with two reticulations `-R2`. Print cost and save result to odt.log. 

```
> supnet -g "(a,(b,(c,d))); ((a,b),(c,d))" -r1 -R2 --HC --randseed 13
   i: ((#1,d),(((((a)#2,c),#2))#1,b)) cost=4
   >: (((#1,d),#2),((((a)#2,c))#1,b)) cost=2
   >: (#2,((#1,d),((((a)#2,c))#1,b))) cost=1
   >: (#2,((#1,d),((a)#2,((c)#1,b)))) cost=0
Cost:0 TopNetworks:20 Steps:703 Climbs:8 HCruns:1 HCTime:0.00291181 Naive:704 Naivetime:0.000615597 DTcnt:2795 Class:TreeChild TimeConsistency:0
```

See more examples in the last section.

## Saving extended newick and inserting prefix in reticulation ids

```
> supnet -g "(a,(b,(c,d))); ((a,b),(c,d))" -r1 -R2 --HC --randseed 13 --saveextnewick --retidprefix=H
   i: ((#H1,d),(((((a)#H2,c),#H2))#H1,b)) cost=4
   >: (((#H1,d),#H2),((((a)#H2,c))#H1,b)) cost=2
   >: (#H2,((#H1,d),((((a)#H2,c))#H1,b))) cost=1
   >: (#H2,((#H1,d),((a)#H2,((c)#H1,b)))) cost=0
Cost:0 TopNetworks:20 Steps:703 Climbs:8 HCruns:1 HCTime:0.00276279 Naive:704 Naivetime:0.000634909 DTcnt:2795 Class:TreeChild TimeConsistency:0
```

Output in odt.tre (not odt.log).
```
> cat odt.tre
(#H2,((#H1,d),((a)#H2,((c)#H1,b))));
(#H2,(#H1,((a)#H2,(((c,d))#H1,b))));
(#H2,(((c)#H1,(#H1,d)),((a)#H2,b)));
(#H2,((#H1,d),((c)#H1,((a)#H2,b))));
(#H2,((c)#H1,((#H1,d),((a)#H2,b))));
(#H2,((#H1,d),((((c)#H1,a))#H2,b)));
(#H2,(#H1,(((((c,d))#H1,a))#H2,b)));
(#H2,((c,#H1),((((d)#H1,a))#H2,b)));
(#H2,(#H1,(((c,d))#H1,((a)#H2,b))));
((#H2,(((c,d))#H1,((a)#H2,b))),#H1);
(#H2,((c,#H1),((d)#H1,((a)#H2,b))));
(#H2,(#H1,((c,(d)#H1),((a)#H2,b))));
(#H2,(((#H1,c),(d)#H1),((a)#H2,b)));
(#H2,((c,(d)#H1),((a)#H2,(#H1,b))));
((#H2,((c,(d)#H1),((a)#H2,b))),#H1);
((#H2,(((c)#H1,d),((a)#H2,b))),#H1);
((#H2,((a)#H2,(((c,d))#H1,b))),#H1);
((#H2,(((((a)#H2,c),d))#H1,b)),#H1);
((#H2,(((c,((a)#H2,d)))#H1,b)),#H1);
((#H2,((((a)#H2,(c,d)))#H1,b)),#H1);
```

## Verbose in algorithms 

 Use `--verbose, -v [123][34]`, where 

  - `0`: quiet mode in HC, BB and ODT-naive
  - `1`: print visited network in HC after each improvement
  - `2`: print visited networks in HC if the cost is equal to the current or improved
  - `3`: print all visited networks in HC 
  - `4`: print basic info in HC, BB, ODT-naive 
  - `5`: print detailed info in HC, BB, ODT-naive 


Print only improvements in HC `-v1`:
```  
> supnet -g "(a,(b,(c,d))); ((a,b),(c,d)); ((a,c),(b,d))" -r1 -R1 --HC -v1 --randseed 13
   i: ((#1,d),(((c,a))#1,b)) cost=5
   >: (((c,#1),d),((a)#1,b)) cost=2
Cost:2 TopNetworks:3 Steps:96 Climbs:4 HCruns:1 HCTime:0.000170946 Naive:97 Naivetime:5.34058e-05 DTcnt:194 Class:TreeChild TimeConsistency:0
```

Print improvements and equal cost networks `-v2`:
```  
> supnet -g "(a,(b,(c,d))); ((a,b),(c,d)); ((a,c),(b,d))" -r1 -R1 --HC -v2 --randseed 13
   i: ((#1,d),(((c,a))#1,b)) cost=5
   >: (((c,#1),d),((a)#1,b)) cost=2
   =: ((a)#1,(((c,#1),d),b)) cost=2
   =: ((a)#1,((c,d),(b,#1))) cost=2
Cost:2 TopNetworks:3 Steps:96 Climbs:4 HCruns:1 HCTime:0.000129938 Naive:97 Naivetime:5.14984e-05 DTcnt:194 Class:TreeChild TimeConsistency:0
```

Print every visited network `-v3`:
```  
> supnet -g "(a,(b,c)); ((a,b),c);" -r1 -R1 --HC -v3 --randseed 13
   i: ((((a)#1,c),#1),b) cost=2
   >: (#1,(((a)#1,c),b)) cost=1
   <: ((((a)#1,c),#1),b) cost=2
   =: ((#1,b),((a)#1,c)) cost=1
   >: (#1,(((a)#1,b),c)) cost=0
   <: (#1,((a)#1,(b,c))) cost=1
   <: (#1,(b,((a)#1,c))) cost=1
   <: ((#1,b),((a)#1,c)) cost=1
   <: (#1,(((a)#1,c),b)) cost=1
   <: (#1,((a)#1,(c,b))) cost=1
   <: (#1,(((a,b))#1,c)) cost=1
   <: (((a)#1,(#1,b)),c) cost=1
   <: ((#1,((a)#1,b)),c) cost=1
   <: (((a)#1,b),(#1,c)) cost=1
   <: ((((a)#1,b),#1),c) cost=1
   <: (#1,(((a)#1,c),b)) cost=1
   <: (#1,((a)#1,(b,c))) cost=1
   <: ((#1,c),((a)#1,b)) cost=1
   <: (#1,(((a,c))#1,b)) cost=2
Cost:0 TopNetworks:1 Steps:18 Climbs:6 HCruns:1 HCTime:5.84126e-05 Naive:19 Naivetime:1.0252e-05 DTcnt:38 Class:TreeChild TimeConsistency:0
```

To show details on DP & BB algorithms use `-v` with 4, 5 and 6. See examples below.


HC with guide trees.
```

> supnet -g "(a,(b,c));((a,(b,g)),(c,(d,(e,(f,h)))));" -R1 --HC -v1 -q2 --guidetree '(a,b);(d,e,(f,g))' -R3 --randseed 13
   i: ((c)#2,((((a,b))#3,(((#3,(((f,g),e),d)))#1,(#1,h))),#2)) cost=12
   >: ((c)#2,((((a,b))#3,(((#3,((f,g),(e,d))))#1,(#1,h))),#2)) cost=11
   >: ((c)#2,(#3,((((a,b))#3,((((f,g),(e,d)))#1,(#1,h))),#2))) cost=10
   >: (((#3,c))#2,((((a,b))#3,((((f,g),(e,d)))#1,(#1,h))),#2)) cost=9
   i: ((((b,a))#3,(((((g,f),(e,d)))#1,((c)#2,(h,#2))),#3)),#1) cost=10
   >: (((((b,a),#2))#3,(((((g,f),(e,d)))#1,((c)#2,h)),#3)),#1) cost=9
Cost:9 TopNetworks:103 Steps:2151 Climbs:10 HCruns:2 HCTime:0.101498 Naive:2153 Naivetime:0.00698018 DTcnt:17224 Class:TreeChild TimeConsistency:0
```

HC with guide clusters.
```
> supnet -g "(a,(b,c));((a,(b,g)),(c,(d,(e,f))));" -R1 --HC -v1 -q2 --guideclusters '(a,b,g);(d,e,f)' -R3 --randseed 13
   i: ((((#3,g))#1,((((b)#3,a),#1),(#2,(e,f)))),((d)#2,c)) cost=4
   >: ((((#3,g))#1,((((b)#3,a),#1),#2)),(((d)#2,(e,f)),c)) cost=3
   >: (((#3,g))#1,(((((b)#3,a),#1),(((d)#2,(e,f)),c)),#2)) cost=2
   >: (((#3,g))#1,(((((b)#3,a),#1),(((d,(e,f)))#2,c)),#2)) cost=1
   >: (((#3,g))#1,(((a,#1),((b)#3,(((d,(e,f)))#2,c))),#2)) cost=0
   i: (#1,((d)#3,((((#2,b),(#3,(f,e))))#1,(((a,g))#2,c)))) cost=5
   >: (#1,((d)#3,((((#2,(b,g)),(#3,(f,e))))#1,((a)#2,c)))) cost=4
   >: (#1,((d)#3,((((#2,g),(#3,(f,e))))#1,(((b,a))#2,c)))) cost=3
   >: (#1,(((d)#3,(f,e)),((((#2,g),#3))#1,(((b,a))#2,c)))) cost=2
   >: (#1,(((d)#3,(f,e)),(((((#2,g),a),#3))#1,((b)#2,c)))) cost=1
   >: (((#2,g),#1),(((d)#3,(f,e)),(((a)#1,((b)#2,c)),#3))) cost=0
Cost:0 TopNetworks:53 Steps:2944 Climbs:16 HCruns:2 HCTime:0.0605209 Naive:2946 Naivetime:0.00981212 DTcnt:23410 Class:TreeChild TimeConsistency:0
```

## Networks/dags comparing and aggregating

 - `--comparedags`: pairwise comparison of all networks or their shapes
 - `--uniquedags`: print unique networks; summary stats printed on stderr
 - `--uniquedagscnts`: as above plus counts of networks
 - `--dagshapes`: two networks are equal if their shapes are isomorphic in --*dags* options (i.e., ignore leaf labels)


Print unique networks from 100 randomly generated. 

```  
> supnet -A2 -R1 -r100 --pnetworks  --randseed 13 | supnet --uniquedags -N-
((a)#1,(#1,b))
(#1,((b)#1,a))
unique=2 all=100
```

Pairwise comparison of random general networks. Here, 
 - `1` - dags are equal,
 - `0` - dags are not equal.

Use `--dagshapes`, to compare shapes.

```  
> supnet -r4 -R1 -A2  --pnetworks --general --comparedags --randseed 13
((a)#1,(#1,b))
(((a)#1,b),#1)
(#1,((b)#1,a))
((b)#1,(#1,a))

((a)#1,(#1,b))	(((a)#1,b),#1)	1
((a)#1,(#1,b))	(#1,((b)#1,a))	0
((a)#1,(#1,b))	((b)#1,(#1,a))	0

(((a)#1,b),#1)	(#1,((b)#1,a))	0
(((a)#1,b),#1)	((b)#1,(#1,a))	0

(#1,((b)#1,a))	((b)#1,(#1,a))	1

```

Print unique random networks with counts
```  
> supnet -r100000 -R1 -A2 --uniquedagscnts --randseed 13
49899	((a)#1,(#1,b))
50101	(#1,((b)#1,a))
unique=2 all=0
```

Print unique random networks with counts using uniform draw
```  
> supnet -r100000 -R1 -A2 --uniquedagscnts --uniformedgesampling --randseed 13
50211	(((b)#1,a),#1)
49789	(((a)#1,b),#1)
unique=2 all=0
```

Print unique random shapes of networks with counts
```  
> supnet -r100000 -R1 -A3 --dagshapes --uniquedagscnts --randseed 13
33530	((((a)#1,c),#1),b)
29997	((((b)#1,c),a),#1)
11760	((((b,c))#1,a),#1)
11528	((a)#1,(#1,(c,b)))
13185	((#1,c),((b)#1,a))
unique=5 all=0
```

Print unique random shapes of networks with counts using uniform draw
```  
> supnet -r100000 -R1 -A3 --dagshapes --uniformedgesampling --uniquedagscnts
28660	((b)#1,(a,(#1,c)))
14394	(((c,b))#1,(#1,a))
14054	(((b)#1,a),(c,#1))
14300	(#1,((b)#1,(a,c)))
28592	(c,(((a)#1,b),#1))
unique=5 all=0
```

## Other options, print, debug, test, etc.

 - `--pspeciesdictionary`: list species dictionary
 - `--pdetailed`: print debug tree and network structures
 - `--reachablenodescnt`: for each v in V(N), print the number of nodes reachable from v
 - `--reachableleafvescnt`: for each v in V(N), print the number of leaves reachable from v
 - `--pstsubtrees`: print species tree subtrees
 - `--maxdisplaytreecachesize`: set limit for the size of the cache of display tree nodes (if DTCACHE is enabled)

Detailed tree/network info (debug)

```
>  supnet -g '(a,((b,a),c))' --pdetailed
Gene trees:
 Nodes(nn)=7 Leaves(lf)=4
 Leaves:0..3 TreeNodes:4..6
 root=6 exactspecies=0
0  Leaf p=6 $0 a subtree=a
1  Leaf p=4 $1 b subtree=b
2  Leaf p=4 $0 a subtree=a
3  Leaf p=5 $2 c subtree=c
4  Tree c=1 c=2 p=5 subtree=(b,a)
5  Tree c=4 c=3 p=6 subtree=((b,a),c)
6  Root c=0 c=5 subtree=(a,((b,a),c))
 parent=  0:6 1:4 2:4 3:5 4:5 5:6 6:32000
 leftchild=  4:1 5:4 6:0
 rightchild=  4:2 5:3 6:5
 lab=  0:0 1:1 2:0 3:2
 depth=uninitialized

```

```
>  supnet -n '((((c)#B,b))#A,(#A,(#B,a)))' --pnetworks --pdetailed
((((c)#B,b))#A,(#A,(#B,a)))
Networks:
 Nodes(nn)=9 Leaves(lf)=3 Reticulations(rt)=2 rtcount()=2
 Leaves:0..2 TreeNodes:3..6 RetNodes:7..8
 rtstartid=7 root=6 exactspecies=1 displaytreemaxid=4
0  Leaf p=7 $0 c subtree=c
1  Leaf p=3 $1 b subtree=b
2  Leaf p=4 $2 a subtree=a
3  Tree c=7 c=1 p=8 subtree=((c)#B,b)
4  Tree c=7 c=2 p=5 subtree=(#B,a)
5  Tree c=8 c=4 p=6 subtree=(#A,(#B,a))
6  Root c=8 c=5 subtree=((((c)#B,b))#A,(#A,(#B,a)))
7  Reti c=0 p=3 p=4 retlabel=#B subtree=(c)#B
8  Reti c=3 p=6 p=5 retlabel=#A subtree=(((c)#B,b))#A
 parent=  0:7 1:3 2:4 3:8 4:5 5:6 6:32000 7:3 8:6
 leftchild=  3:7 4:7 5:8 6:8
 rightchild=  3:1 4:2 5:4 6:5
 lab=  0:0 1:1 2:2
 retchild=  7:0 8:3
 retparent=  7:4 8:5
 spid2retlabel=  7:#B 8:#A

```


### Network iterator

Insert one reticulation in allowed positions.

```
> supnet  -n "((a,b),c)"  --netretiterator  --guideclusters '(a,c)'
((((c)#1,a),b),#1)
((#1,b),((a)#1,c))
```

```
> supnet  -n "((a,b),c)"  --netretiterator  --guidetree '(a,b)'
(#1,(((a,b))#1,c))
(((c)#1,(a,b)),#1)
((c)#1,((a,b),#1))
(((a,b))#1,(#1,c))
```

```
> supnet  -n "((a,b),c)"  --netretiterator   --timeconsistent
((#1,b),((a)#1,c))
((a,#1),((b)#1,c))
```

### Check network classes in the result:
```
supnet -N odt.log --detectclass
```

### Display tree cache 

By default DTCACHE is enabled since version 0.10. Compile with macro NODTCACHE to obtain non-cached version. Option `--maxdisplaytreecachesize SIZE` will limit the number of stored nodes in the cache. The default is `1000000` nodes.

```
supnet  --maxdisplaytreecachesize 5000000 ...
```

### Limit the number of improvements

The number of HC steps in improvements 
```
supnet --hcstopclimb=10 ...
```

# Examples 


### Stats in HC runs.

Do not print detailed stats for individual HC runs.
```
> supnet -g "(a,(b,(c,d))); ((a,b),(c,d)); ((a,c),(b,d)); ((a,d),e)" -r8 -R1 --HC -v0 --randseed 13
Cost:3 TopNetworks:1 Steps:2066 Climbs:43 HCruns:8 HCTime:0.0028038 Naive:2074 Naivetime:0.00117135 DTcnt:4148 Class:TreeChild TimeConsistency:0
```

Print stats after completed HC run for each initial network when only when new optimal network is found.
```
> supnet -g "(a,(b,(c,d))); ((a,b),(c,d)); ((a,c),(b,d)); ((a,d),e)" -r8 -R1 --HC --hcrunstats -v0 --randseed 13
1.  Steps:0 Climbs:0 New optimal cost: 4
5.  Steps:0 Climbs:0 New optimal cost: 3
Cost:3 TopNetworks:1 Steps:2066 Climbs:43 HCruns:8 HCTime:0.0028069 Naive:2074 Naivetime:0.00115442 DTcnt:4148 Class:TreeChild TimeConsistency:0
```

Print stats after completed HC run for each initial network when a new optimal network including equal networks is found.
```
> supnet -g "(a,(b,(c,d))); ((a,b),(c,d)); ((a,c),(b,d)); ((a,d),e)" -r8 -R1 --HC --hcrunstatsext -v0 --randseed 13
1.  Steps:0 Climbs:0 New optimal cost: 4
5.  Steps:0 Climbs:0 New optimal cost: 3
Cost:3 TopNetworks:1 Steps:2066 Climbs:43 HCruns:8 HCTime:0.0029521 Naive:2074 Naivetime:0.0012331 DTcnt:4148 Class:TreeChild TimeConsistency:0
```

Print stats after completed HC run.
```
> supnet -g "(a,(b,(c,d))); ((a,b),(c,d)); ((a,c),(b,d)); ((a,d),e)" -r8 -R1 --HC --hcrunstatsalways -v0 --randseed 13
1.  Steps:0 Climbs:0 New optimal cost: 4
2.  Steps:0 Climbs:0 New optimal networks:0 Total:5
3.  Steps:0 Climbs:0 New optimal networks:0 Total:5
4.  Steps:0 Climbs:0 New optimal networks:0 Total:5
5.  Steps:0 Climbs:0 New optimal cost: 3
6.  Steps:0 Climbs:0 New optimal networks:0 Total:1
7.  Steps:0 Climbs:0 New optimal networks:0 Total:1
8.  Steps:0 Climbs:0 New optimal networks:0 Total:1
Cost:3 TopNetworks:1 Steps:2066 Climbs:43 HCruns:8 HCTime:0.0028584 Naive:2074 Naivetime:0.00119877 DTcnt:4148 Class:TreeChild TimeConsistency:0
```

#### Run naive odt computation threshold `-t THRESHOLD`

Run naive odt computation, when the number of reticulations in a network is `< THRESHOLD`; otherwise run BB & DP algorithms.

Higher `-t NUM`, more naive algorithm computations.

```
> supnet -g "(a,(b,(c,d)))" -n "(#1,((b,((d,c))#1),a))" -t2 --HC -v35 --randseed 13
HC start: hcusenaive=0 runnaiveleqrt=13 tailmove=1
ODT-cost-naive 1 tree(s);rt=1(-t) cost=0
   i: (#1,((b,((d,c))#1),a)) cost=0
ODT-cost-naive 1 tree(s);rt=1(-t) cost=1
   <: (#1,(((d,b),(c)#1),a)) cost=1
ODT-cost-naive 1 tree(s);rt=1(-t) cost=0
   =: (#1,((b,(d,(c)#1)),a)) cost=0
ODT-cost-naive 1 tree(s);rt=1(-t) cost=1
   <: (#1,((d,(b,(c)#1)),a)) cost=1
ODT-cost-naive 1 tree(s);rt=1(-t) cost=2
   <: (#1,((b,(c)#1),(d,a))) cost=2
ODT-cost-naive 1 tree(s);rt=1(-t) cost=1
   <: ((d,#1),((b,(c)#1),a)) cost=1
ODT-cost-naive 1 tree(s);rt=1(-t) cost=3
   <: (#1,(d,((b,(c)#1),a))) cost=3
ODT-cost-naive 1 tree(s);rt=1(-t) cost=1
   <: (#1,(((b,c),(d)#1),a)) cost=1
ODT-cost-naive 1 tree(s);rt=1(-t) cost=0
   =: (#1,((b,((d)#1,c)),a)) cost=0
...
   <: (#1,(b,(a,(d,(c)#1)))) cost=1
ODT-cost-naive 1 tree(s);rt=1(-t) cost=1
   <: ((#1,(d,(c)#1)),(b,a)) cost=1
ODT-cost-naive 1 tree(s);rt=1(-t) cost=1
   <: (#1,((b,a),(d,(c)#1))) cost=1
ODT-cost-naive 1 tree(s);rt=1(-t) cost=2
   <: (#1,(b,((d,a),(c)#1))) cost=2
ODT-cost-naive 1 tree(s);rt=1(-t) cost=2
   <: (#1,(b,(d,((c)#1,a)))) cost=2
ODT-cost-naive 1 tree(s);rt=1(-t) cost=1
   <: (#1,((b,a),(d,(c)#1))) cost=1
ODT-cost-naive 1 tree(s);rt=1(-t) cost=1
   <: (#1,(b,((d,(c)#1),a))) cost=1
ODT-cost-naive 1 tree(s);rt=1(-t) cost=2
   <: (#1,(b,(d,((c,a))#1))) cost=2
HC run completed: (#1,((b,(d,(c)#1)),a)) cost=0
Stats saved to odt.dat
Best networks saved to odt.log
Cost:0 TopNetworks:19 Steps:410 Climbs:2 HCruns:1 HCTime:0.00152683 Naive:411 Naivetime:0.000181913 DTcnt:809 Class:TreeChild TimeConsistency:0
```

Lower `-t`, more BB & DP computations
```
> supnet -g "(a,(b,(c,d)))" -n "(#1,((b,((d,c))#1),a))" -t1 --HC -v35 --randseed 13
HC start: hcusenaive=0 runnaiveleqrt=13 tailmove=1
ODT-cost-naive 1 tree(s);rt=1(-t) cost=0
   i: (#1,((b,((d,c))#1),a)) cost=0
ODT-cost-naive 1 tree(s);rt=1(-t) cost=1
   <: (#1,(((d,b),(c)#1),a)) cost=1
ODT-cost-naive 1 tree(s);rt=1(-t) cost=0
   =: (#1,((b,(d,(c)#1)),a)) cost=0
ODT-cost-naive 1 tree(s);rt=1(-t) cost=1
   <: (#1,((d,(b,(c)#1)),a)) cost=1
ODT-cost-naive 1 tree(s);rt=1(-t) cost=2
   <: (#1,((b,(c)#1),(d,a))) cost=2
ODT-cost-naive 1 tree(s);rt=1(-t) cost=1
   <: ((d,#1),((b,(c)#1),a)) cost=1
ODT-cost-naive 1 tree(s);rt=1(-t) cost=3
   <: (#1,(d,((b,(c)#1),a))) cost=3
ODT-cost-naive 1 tree(s);rt=1(-t) cost=1
   <: (#1,(((b,c),(d)#1),a)) cost=1
ODT-cost-naive 1 tree(s);rt=1(-t) cost=0
   =: (#1,((b,((d)#1,c)),a)) cost=0
...
   <: (#1,(b,(a,(d,(c)#1)))) cost=1
ODT-cost-naive 1 tree(s);rt=1(-t) cost=1
   <: ((#1,(d,(c)#1)),(b,a)) cost=1
ODT-cost-naive 1 tree(s);rt=1(-t) cost=1
   <: (#1,((b,a),(d,(c)#1))) cost=1
ODT-cost-naive 1 tree(s);rt=1(-t) cost=2
   <: (#1,(b,((d,a),(c)#1))) cost=2
ODT-cost-naive 1 tree(s);rt=1(-t) cost=2
   <: (#1,(b,(d,((c)#1,a)))) cost=2
ODT-cost-naive 1 tree(s);rt=1(-t) cost=1
   <: (#1,((b,a),(d,(c)#1))) cost=1
ODT-cost-naive 1 tree(s);rt=1(-t) cost=1
   <: (#1,(b,((d,(c)#1),a))) cost=1
ODT-cost-naive 1 tree(s);rt=1(-t) cost=2
   <: (#1,(b,(d,((c,a))#1))) cost=2
HC run completed: (#1,((b,(d,(c)#1)),a)) cost=0
Stats saved to odt.dat
Best networks saved to odt.log
Cost:0 TopNetworks:19 Steps:410 Climbs:2 HCruns:1 HCTime:0.00163269 Naive:411 Naivetime:0.000174761 DTcnt:809 Class:TreeChild TimeConsistency:0
```

### NNI vs Tail Moves

HC via Tail Moves (default)
```
> supnet -g "(a,(b,(c,d))); ((a,b),(c,d)); ((a,c),(b,d))" -r1 -R1 --HC -v3 --randseed 13
   i: ((#1,d),(((c,a))#1,b)) cost=5
   >: (((c,#1),d),((a)#1,b)) cost=2
   <: ((#1,(c,d)),((a)#1,b)) cost=3
   <: ((c,(#1,d)),((a)#1,b)) cost=3
   <: ((#1,d),(c,((a)#1,b))) cost=5
   <: ((#1,d),((c,(a)#1),b)) cost=5
   <: ((#1,d),((a)#1,(c,b))) cost=6
   <: ((#1,d),(((c,a))#1,b)) cost=5
   <: ((c,(d,#1)),((a)#1,b)) cost=3
   <: (((c,d),#1),((a)#1,b)) cost=3
   <: ((c,d),(((a)#1,b),#1)) cost=3
   <: ((c,d),((a)#1,(b,#1))) cost=3
   <: (d,((c,#1),((a)#1,b))) cost=5
   <: (d,(((c,#1),(a)#1),b)) cost=6
   <: (d,((a)#1,((c,#1),b))) cost=6
   <: (((c,d),#1),((a)#1,b)) cost=3
   <: ((c,(#1,d)),((a)#1,b)) cost=3
   <: ((c,#1),(((a)#1,b),d)) cost=3
   <: ((c,#1),(((a)#1,d),b)) cost=4
...
   <: (((c,a))#1,(d,(b,#1))) cost=4
   <: ((a)#1,(c,((b,d),#1))) cost=4
   <: ((a)#1,(c,(b,(#1,d)))) cost=4
   <: ((a)#1,(c,((b,#1),d))) cost=3
   <: (((a)#1,d),(c,(b,#1))) cost=5
   <: ((a)#1,((c,(b,#1)),d)) cost=4
   <: (((a,d))#1,(c,(b,#1))) cost=6
   <: ((a)#1,(((c,d),b),#1)) cost=3
   <: ((a)#1,(b,((c,d),#1))) cost=3
   <: (((c,d),(a)#1),(b,#1)) cost=3
   <: ((((c,d),a))#1,(b,#1)) cost=4
   <: ((a)#1,((c,(b,#1)),d)) cost=4
   <: ((a)#1,(c,(d,(b,#1)))) cost=3
   <: (((a)#1,(b,#1)),(c,d)) cost=3
   <: ((c,d),(((a)#1,b),#1)) cost=3
   <: ((c,((a)#1,d)),(b,#1)) cost=3
   <: (((a)#1,(c,d)),(b,#1)) cost=3
   <: ((c,d),((a)#1,(b,#1))) cost=3
Cost:2 TopNetworks:3 Steps:96 Climbs:4 HCruns:1 HCTime:0.000226021 Naive:97 Naivetime:5.50747e-05 DTcnt:194 Class:TreeChild TimeConsistency:0
```

HC via NNI moves 

```
> supnet -g "(a,(b,(c,d))); ((a,b),(c,d)); ((a,c),(b,d))" -r1 -R1 --HC --hcnnimove -v3 --randseed 13
   i: ((#1,d),(((c,a))#1,b)) cost=5
   >: ((d,(((c,a))#1,b)),#1) cost=4
   <: (((((c,a))#1,b),#1),d) cost=6
   <: ((#1,d),(((c,a))#1,b)) cost=5
   =: ((((c,a))#1,(b,d)),#1) cost=4
   =: ((b,(d,((c,a))#1)),#1) cost=4
   <: (((d,((c,a))#1),#1),b) cost=5
   <: ((#1,b),(d,((c,a))#1)) cost=5
Cost:4 TopNetworks:3 Steps:7 Climbs:4 HCruns:1 HCTime:5.17368e-05 Naive:8 Naivetime:7.86781e-06 DTcnt:16 Class:TreeChild TimeConsistency:0
```

#### Default tree-child networks and Tail Moves:
```
> supnet -g "(a,(b,(c,d))); ((a,b),(c,d)); ((a,c),(b,d))" -r1 -R3 --HC -v2 --randseed 12
   i: ((#3,((c)#3,((((d)#1,(b,#1)))#2,a))),#2) cost=4
   >: (((c)#3,((((d)#1,((#3,b),#1)))#2,a)),#2) cost=3
   =: (((c)#3,((((d)#1,(#3,(b,#1))))#2,a)),#2) cost=3
   >: (((c)#3,(((((#3,d))#1,(b,#1)))#2,a)),#2) cost=1
   >: ((((c)#3,(((((#3,d))#1,b))#2,a)),#1),#2) cost=0
   =: (#1,(((c)#3,(((((#3,d))#1,b))#2,a)),#2)) cost=0
   =: (((c)#3,((((((#3,d))#1,b))#2,a),#1)),#2) cost=0
   =: (((((((#3,d))#1,b))#2,((c)#3,a)),#1),#2) cost=0
   =: (#1,((((((#3,d))#1,b))#2,((c)#3,a)),#2)) cost=0
Cost:0 TopNetworks:5 Steps:121 Climbs:8 HCruns:1 HCTime:0.000538349 Naive:122 Naivetime:0.000247478 DTcnt:966 Class:TreeChild TimeConsistency:0
```

```
> cat odt.log
((((c)#3,(((((#3,d))#1,b))#2,a)),#1),#2)
(#1,(((c)#3,(((((#3,d))#1,b))#2,a)),#2))
(((c)#3,((((((#3,d))#1,b))#2,a),#1)),#2)
(((((((#3,d))#1,b))#2,((c)#3,a)),#1),#2)
(#1,((((((#3,d))#1,b))#2,((c)#3,a)),#2))
```

#### Relaxed networks. 
```
> supnet -g "(a,(b,(c,d))); ((a,b),(c,d)); ((a,c),(b,d))" -r1 -R3 --HC -v1 --randseed 13 --relaxed
   i: ((#1)#3,((#3,d),(((((a)#2,c),#2))#1,b))) cost=4
   >: ((#1)#3,(((((a)#2,c),#3),d),((#2)#1,b))) cost=1
   >: (((c,#1))#3,((a)#2,((#3,d),((#2)#1,b)))) cost=0
Cost:0 TopNetworks:2 Steps:221 Climbs:6 HCruns:1 HCTime:0.000752449 Naive:222 Naivetime:0.000427246 DTcnt:1770 Class:Relaxed TimeConsistency:0
```

```
> cat odt.log
(((c,#1))#3,((a)#2,((#3,d),((#2)#1,b))))
((((c,#1))#3,((#3,d),((#2)#1,b))),(a)#2)
```

#### General networks run (`--general`)
```
> supnet -g "(a,(b,(c,d))); ((a,b),(c,d)); ((a,c),(b,d))" -r1 -R3 --HC -v2 --randseed 12 --general
   i: ((#2,((#3)#2,a)),((d)#1,(b,((c)#3,#1)))) cost=1
   =: (((#3)#2,a),((d)#1,((#2,b),((c)#3,#1)))) cost=1
   =: (((#3)#2,a),((d)#1,(b,(#2,((c)#3,#1))))) cost=1
   =: (((#3)#2,a),(#2,((d)#1,(b,((c)#3,#1))))) cost=1
   =: (((#3)#2,a),((#2,(d)#1),(b,((c)#3,#1)))) cost=1
   =: (((#3)#2,a),((d)#1,(#2,(b,((c)#3,#1))))) cost=1
   =: (((#3)#2,a),((d)#1,(b,((#2,(c)#3),#1)))) cost=1
   =: (((#3)#2,a),((d)#1,(b,((c)#3,(#2,#1))))) cost=1
   =: (((#3)#2,a),(((#2,d))#1,(b,((c)#3,#1)))) cost=1
   =: (#2,(((d)#1,(b,((c)#3,#1))),((#3)#2,a))) cost=1
   =: (#2,(((d)#1,((#3)#2,a)),(b,((c)#3,#1)))) cost=1
   =: (#2,((d)#1,((b,((c)#3,#1)),((#3)#2,a)))) cost=1
   =: ((#2,((#3)#2,a)),(b,((d)#1,((c)#3,#1)))) cost=1
   =: ((#2,((#3)#2,a)),((b,(d)#1),((c)#3,#1))) cost=1
   =: (((#2,((#3)#2,a)),(d)#1),(b,((c)#3,#1))) cost=1
   =: ((d)#1,((#2,((#3)#2,a)),(b,((c)#3,#1)))) cost=1
   =: ((((d)#1,#2),((#3)#2,a)),(b,((c)#3,#1))) cost=1
   =: ((#2,((d)#1,((#3)#2,a))),(b,((c)#3,#1))) cost=1
   =: ((#2,(((d)#1,(#3)#2),a)),(b,((c)#3,#1))) cost=1
...
   =: ((((((b,#1))#2,a),(((c)#3,d))#1),#3),#2) cost=0
   =: (((((b,((c)#3,#1)))#2,a),((d)#1,#3)),#2) cost=0
   =: (((b,((d)#1,#2)),((#3)#2,a)),((c)#3,#1)) cost=0
   =: (((b,((d)#1,#2)),((#3)#2,a)),((c,#1))#3) cost=0
   =: ((((b,(d)#1),#2),((#3)#2,a)),((c,#1))#3) cost=0
   =: ((((d)#1,(b,#2)),((#3)#2,a)),((c,#1))#3) cost=0
   =: ((((b,((d,#2))#1),((#3)#2,a)),#1),(c)#3) cost=0
   =: ((#2,((c)#3,d)),(a,(((b)#2,#1),(#3)#1))) cost=0
   =: (((#2,((c)#3,d)),(#3)#1),(a,((b)#2,#1))) cost=0
   =: ((#3)#1,((#2,((c)#3,d)),(a,((b)#2,#1)))) cost=0
   =: ((#3)#1,((#2,((c)#3,d)),((a,(b)#2),#1))) cost=0
   =: ((#3)#1,((#2,((c)#3,d)),((b)#2,(a,#1)))) cost=0
   =: ((((b)#2,((c)#3,d)))#1,((#2,a),(#1,#3))) cost=0
   =: ((c)#3,(((b)#2,(a,#1)),(((#3)#1,d),#2))) cost=0
   =: ((c)#3,((a,((b)#2,#1)),(((#3)#1,d),#2))) cost=0
   =: ((c)#3,(((a,(b)#2),#1),(((#3)#1,d),#2))) cost=0
   =: (#3,(((#2,d))#1,((((c,#1))#3,(b)#2),a))) cost=0
   =: ((((#3,d))#1,(((b)#2,#1),((c)#3,a))),#2) cost=0
Cost:0 TopNetworks:680 Steps:52673 Climbs:4 HCruns:1 HCTime:3.65878 Naive:52674 Naivetime:0.10831 DTcnt:420111 Class:General TimeConsistency:0
```

```
> cat odt.log
(#2,(((d)#1,(b,((c,#1))#3)),((#3)#2,a)))
(#2,((b,((d)#1,((c,#1))#3)),((#3)#2,a)))
(#2,(((b,(d)#1),((c,#1))#3),((#3)#2,a)))
(#2,(((c,#1))#3,((((b,(d)#1),#3))#2,a)))
(#2,(((c,#1))#3,((((d)#1,(b,#3)))#2,a)))
(#2,(((c,#1))#3,(((b,((d)#1,#3)))#2,a)))
((#2,(d)#1),(((c,#1))#3,(((b,#3))#2,a)))
((#2,#3),(((c,#1))#3,(((b,(d)#1))#2,a)))
(#2,(((c,#1))#3,((((b,(d)#1))#2,#3),a)))
((((c,#1))#3,#2),((((b,(d)#1))#2,#3),a))
((#2,((((b,(d)#1))#2,#3),a)),((c,#1))#3)
((#2,(d)#1),(((c,#1))#3,(((b)#2,#3),a)))
(((((c,#1))#3,#2),(d)#1),(((b)#2,#3),a))
((#2,(((c,#1))#3,(d)#1)),(((b)#2,#3),a))
((((c,#1))#3,(#2,(d)#1)),(((b)#2,#3),a))
(((#2,(d)#1),(((b)#2,#3),a)),((c,#1))#3)
(((#2,d))#1,(((c,#1))#3,(((b)#2,#3),a)))
(((#2,#3),(d)#1),(((c,#1))#3,((b)#2,a)))
((#2,((d)#1,#3)),(((c,#1))#3,((b)#2,a)))
...
(((((b,#1))#2,a),((((c)#3,d))#1,#3)),#2)
((((((b,#1))#2,a),(((c)#3,d))#1),#3),#2)
(((((b,((c)#3,#1)))#2,a),((d)#1,#3)),#2)
(((b,((d)#1,#2)),((#3)#2,a)),((c)#3,#1))
(((b,((d)#1,#2)),((#3)#2,a)),((c,#1))#3)
((((b,(d)#1),#2),((#3)#2,a)),((c,#1))#3)
((((d)#1,(b,#2)),((#3)#2,a)),((c,#1))#3)
((((b,((d,#2))#1),((#3)#2,a)),#1),(c)#3)
((#2,((c)#3,d)),(a,(((b)#2,#1),(#3)#1)))
(((#2,((c)#3,d)),(#3)#1),(a,((b)#2,#1)))
((#3)#1,((#2,((c)#3,d)),(a,((b)#2,#1))))
((#3)#1,((#2,((c)#3,d)),((a,(b)#2),#1)))
((#3)#1,((#2,((c)#3,d)),((b)#2,(a,#1))))
((((b)#2,((c)#3,d)))#1,((#2,a),(#1,#3)))
((c)#3,(((b)#2,(a,#1)),(((#3)#1,d),#2)))
((c)#3,((a,((b)#2,#1)),(((#3)#1,d),#2)))
((c)#3,(((a,(b)#2),#1),(((#3)#1,d),#2)))
(#3,(((#2,d))#1,((((c,#1))#3,(b)#2),a)))
((((#3,d))#1,(((b)#2,#1),((c)#3,a))),#2)
```

Note that using a broader network class does not guarantee a better cost.


#### General and time consistent networks `--general --timeconsistent`

```
> supnet -g "(a,(b,(c,d))); ((a,b),(c,d)); ((a,c),(b,d))" -r1 -R3 --HC -v2 --randseed 12 --general --timeconsistent
   i: ((c,((#2)#1,((#3)#2,((d)#3,a)))),(b,#1)) cost=4
   >: (((#2)#1,((#3)#2,((d)#3,a))),(b,(c,#1))) cost=3
   >: (((c,(#2)#1),((#3)#2,((d)#3,a))),(b,#1)) cost=2
   =: (((#2)#1,((c,(#3)#2),((d)#3,a))),(b,#1)) cost=2
   =: (((#2)#1,((#3)#2,((c,(d)#3),a))),(b,#1)) cost=2
   =: (((c,(#2)#1),(b,#1)),((#3)#2,((d)#3,a))) cost=2
   =: ((c,(#2)#1),(((#3)#2,((d)#3,a)),(b,#1))) cost=2
   =: (((((#3)#2,c),(#2)#1),((d)#3,a)),(b,#1)) cost=2
   =: (((c,(#2)#1),((d)#3,a)),(((#3)#2,b),#1)) cost=2
   =: (((((d)#3,c),(#2)#1),((#3)#2,a)),(b,#1)) cost=2
   =: (((c,(#2)#1),((#3)#2,a)),(((d)#3,b),#1)) cost=2
   =: (((#2)#1,((#3)#2,a)),(((d)#3,b),(c,#1))) cost=2
   =: (((#2)#1,((c,(#3)#2),a)),(((d)#3,b),#1)) cost=2
   >: ((((c,#2))#1,((#3)#2,a)),(((d)#3,b),#1)) cost=1
   =: (((#2)#1,(((c,#3))#2,a)),(((d)#3,b),#1)) cost=1
   =: ((((c,#2))#1,a),(((#3)#2,((d)#3,b)),#1)) cost=1
   =: (((#2)#1,a),((((c,#3))#2,((d)#3,b)),#1)) cost=1
   =: (((c)#1,a),(((#3)#2,((d)#3,b)),(#1,#2))) cost=1
   =: ((((((d)#3,c),#2))#1,a),(((#3)#2,b),#1)) cost=1
...
   =: (((#3)#2,(((c)#1,(d)#3),b)),((#1,a),#2)) cost=1
   =: (((c)#1,(d)#3),((((#1,#3))#2,b),(a,#2))) cost=1
   =: ((((c)#1,(d)#3),#2),((((#1,#3))#2,b),a)) cost=1
   =: (((((#1,#3))#2,b),((c)#1,(d)#3)),(a,#2)) cost=1
   =: ((((c)#1,(d)#3),(a,#2)),(((#1,#3))#2,b)) cost=1
   =: ((((#1,#3))#2,b),((((c)#1,(d)#3),a),#2)) cost=1
   =: ((((#1,#3))#2,b),(a,(((c)#1,(d)#3),#2))) cost=1
   =: (((((c)#1,(d)#3),((#1,#3))#2),b),(a,#2)) cost=1
   =: ((((#1,#3))#2,(((c)#1,(d)#3),b)),(a,#2)) cost=1
   =: (((c)#1,(d)#3),(((#1)#2,(b,#3)),(a,#2))) cost=1
   =: ((((#1)#2,(b,#3)),((c)#1,(d)#3)),(a,#2)) cost=1
   =: (((((c)#1,(d)#3),(#1)#2),(b,#3)),(a,#2)) cost=1
   =: (((#1)#2,(((c)#1,(d)#3),(b,#3))),(a,#2)) cost=1
   =: (((#1)#2,((c)#1,(d)#3)),((b,#3),(a,#2))) cost=1
   =: ((((a,#1))#2,((c)#1,(d)#3)),((b,#3),#2)) cost=1
   =: (((#1)#2,(d)#3),((b,#3),(((c)#1,a),#2))) cost=1
   =: (((#1)#2,(d)#3),((((c)#1,b),#3),(a,#2))) cost=1
   >: (((#1)#2,(((c)#1,d))#3),((b,#3),(a,#2))) cost=0
Cost:0 TopNetworks:1 Steps:176 Climbs:10 HCruns:1 HCTime:0.00292063 Naive:177 Naivetime:0.000544071 DTcnt:1414 Class:General TimeConsistency:1
```

### Global cache of networks in HC runs

If `--globaldagcache` is set then networks already visited will be omitted in HC.
This feature is not optimized. For large sets of networks, it may be slower that non-cache option.

### odt.dat in labelled format `--odtlabelled`

```
> supnet -g "(a,(b,(c,d))); ((a,b),(c,d))" -r1 -R2 --randseed 13 --HC --odtlabelled -v0
Cost:0 TopNetworks:20 Steps:703 Climbs:8 HCruns:1 HCTime:0.00272775 Naive:704 Naivetime:0.000603914 DTcnt:2795 Class:TreeChild TimeConsistency:0
```

```
> cat odt.dat
optcost=0
time=0.00272799
hctime=0.00272775
mergetime=2.38419e-07
topnets=20
class=0
timeconsistency=0
improvements=8
steps=703
bbruns=0
startnets=1
memoryMB=6
dtcnt=2795
naivetime=0.000603914
naivecnt=704
bbnaivecnt=0
bbnaivetime=0
bbdpcnt=0
bbdptime=0
randseed=13
```

### More examples

Larger instance; tree-child search:
```  
> supnet -g "((i,c),((a,d),(f,(b,((g,(e,j)),h)))));((i,h),((c,f),((a,(d,e)),(j,(g,b)))));(((f,(b,d)),((j,g),(e,i))),(c,(h,a)));(((f,(i,(j,d))),(c,b)),(((h,a),g),e));((((h,e),((c,f),a)),(d,b)),((g,i),j))" -R8 -q1 --HC --hcrunstats --randseed 13
   i: ((((((d)#3,((((#3,((a)#7,h)))#6,c),#8)))#5,e))#2,((((((i)#4,(#7,(#2,g))))#1,(b,(#1,(#4,((j)#8,f))))),#5),#6)) cost=49
   >: ((((((d)#3,((((#3,((a)#7,h)))#6,c),#8)))#5,e))#2,(((((((#2,i))#4,(#7,g)))#1,(b,(#1,(#4,((j)#8,f))))),#5),#6)) cost=48
   >: ((((((d)#3,((((#3,((a)#7,h)))#6,c),#8)))#5,(#7,e)))#2,(((((((#2,i))#4,g))#1,(b,(#1,(#4,((j)#8,f))))),#5),#6)) cost=47
   >: ((((((d)#3,((((#3,((a)#7,h)))#6,c),#8)))#5,(#7,e)))#2,(((#4,(((((#2,i))#4,g))#1,(b,(#1,((j)#8,f))))),#5),#6)) cost=46
   >: ((((((d)#3,((((#3,((a)#7,h)))#6,c),#8)))#5,(#7,e)))#2,(((#4,(((((#2,i))#4,g))#1,((#1,b),((j)#8,f)))),#5),#6)) cost=45
   >: ((((((d)#3,(((((#3,b),((a)#7,h)))#6,c),#8)))#5,(#7,e)))#2,(((#4,(((((#2,i))#4,g))#1,(#1,((j)#8,f)))),#5),#6)) cost=43
   >: ((((((d)#3,(((((#3,b),((a)#7,h)))#6,c),#8)))#5,(#7,e)))#2,((#1,((#4,(((((#2,i))#4,g))#1,((j)#8,f))),#5)),#6)) cost=42
   >: (((((((d,b))#3,((((#3,((a)#7,h)))#6,c),#8)))#5,(#7,e)))#2,((#1,((#4,(((((#2,i))#4,g))#1,((j)#8,f))),#5)),#6)) cost=41
   >: (((((((d,b))#3,((((#3,((a)#7,h)))#6,c),#8)))#5,(#7,e)))#2,((#1,((#4,((((i)#4,(#2,g)))#1,((j)#8,f))),#5)),#6)) cost=39
   >: (((((((d,b))#3,(((#3,((a)#7,h)))#6,c)))#5,(#7,e)))#2,((#1,((#4,((((i)#4,(#2,(g,#8))))#1,((j)#8,f))),#5)),#6)) cost=35
   >: (((((((d,b))#3,((((a)#7,h))#6,c)))#5,(#7,(#3,e))))#2,((#1,((#4,((((i)#4,(#2,(g,#8))))#1,((j)#8,f))),#5)),#6)) cost=33
   >: ((((((((#7,d),b))#3,((((a)#7,h))#6,c)))#5,(#3,e)))#2,((#1,((#4,((((i)#4,(#2,(g,#8))))#1,((j)#8,f))),#5)),#6)) cost=32
   >: ((((((((#7,d),b))#3,((((a)#7,h))#6,c)))#5,(#3,e)))#2,((#1,(#4,(((((i)#4,(#2,(g,#8))))#1,((j)#8,f)),#5))),#6)) cost=31
   >: ((((((((#7,d),b))#3,((((a)#7,h))#6,(#4,c))))#5,(#3,e)))#2,((#1,(((((i)#4,(#2,(g,#8))))#1,((j)#8,f)),#5)),#6)) cost=30
   >: ((((((((#7,d),b))#3,((((a)#7,h))#6,(#4,c))))#5,(#3,(e,#6))))#2,(#1,(((((i)#4,(#2,(g,#8))))#1,((j)#8,f)),#5))) cost=29
   >: (((((((d,b))#3,((((a)#7,h))#6,(#4,c))))#5,(#3,(e,#6))))#2,(#1,((((#7,((i)#4,(#2,(g,#8)))))#1,((j)#8,f)),#5))) cost=28
   >: (((((((d,b))#3,((((a)#7,h))#6,(#4,c))))#5,(e,#6)))#2,(#1,((((#7,(#3,((i)#4,(#2,(g,#8))))))#1,((j)#8,f)),#5))) cost=27
   >: ((((((((d,b))#3,((((a)#7,h))#6,(#4,c))))#5,(e,#6)))#2,((((#7,(#3,((i)#4,(#2,(g,#8))))))#1,((j)#8,f)),#5)),#1) cost=25
   >: (((((((((a)#7,h))#6,(#4,c)))#5,(e,#6)))#2,((((#7,(#3,((i)#4,(#2,(g,#8))))))#1,(((d,b))#3,((j)#8,f))),#5)),#1) cost=24
   >: (((((((((a)#7,h))#6,(#4,c)))#5,(e,#6)))#2,((((#7,(#3,((i)#4,(#2,(g,#8))))))#1,(((((j)#8,d),b))#3,f)),#5)),#1) cost=23
   >: ((((e,#6))#2,((((#7,(#3,((i)#4,(#2,(g,#8))))))#1,(((((((((a)#7,h))#6,(#4,c)))#5,((j)#8,d)),b))#3,f)),#5)),#1) cost=22
   >: (((((#7,((((((a)#7,h))#6,(#4,c)))#5,((j)#8,(((e,#6))#2,d)))),b))#3,((((#3,(#2,((i)#4,(g,#8)))))#1,f),#5)),#1) cost=20
   >: (((((#7,(((#4,c))#5,((j)#8,d))),b))#3,((((((e)#2,a))#7,h))#6,((((#3,((#2,(((i)#4,g),#8)),#6)))#1,f),#5))),#1) cost=19
   >: (((((((e)#2,a))#7,h))#6,(((((i)#4,(#3,((#2,(g,#8)),#6))))#1,((#7,(((j)#8,d),(b)#3)),(((#4,c))#5,f))),#5)),#1) cost=18
   >: (((((((e)#2,a))#7,h))#6,((#4,((((i)#4,(#3,((#2,(g,#8)),#6))))#1,((#7,(((j)#8,d),(b)#3)),((c)#5,f)))),#5)),#1) cost=17
   >: (((#4,(((#7,((((((((a)#7,h))#6,i))#4,e))#2,(d,(b)#3))),(j)#8),(((#3,((#2,(g,#8)),#6)))#1,((c)#5,f)))),#5),#1) cost=16
1.  Steps:0 Climbs:0 New optimal cost: 16
Cost:16 TopNetworks:60 Steps:46020 Climbs:52 HCruns:1 HCTime:13.8784 Naive:46021 Naivetime:8.8291 DTcnt:11781376 Class:TreeChild TimeConsistency:0
```

Recommended with large HC-runs using quasi-consensus rand networks
```  
> supnet -g "(a,(b,(c,d))); ((a,b),(c,d))" -q1000 -R3 --HC -v1 --randseed 13 | tail -20
   >: (((((c)#2,b))#1,((((a,#1))#3,d),#2)),#3) cost=0
   i: ((((#2,((d)#3,((b)#2,c))))#1,(#1,a)),#3) cost=1
   >: ((((#2,((d)#3,c)))#1,(#1,((b)#2,a))),#3) cost=0
   i: ((#1,(#3,((c)#2,((d)#3,((b)#1,a))))),#2) cost=3
   >: ((#1,(((#3,c))#2,((d)#3,((b)#1,a)))),#2) cost=1
   >: ((((#1,(#3,c)))#2,((d)#3,((b)#1,a))),#2) cost=0
   i: (((a)#3,(#1,((((c)#1,b))#2,(#2,d)))),#3) cost=3
   >: ((#2,((a)#3,(#1,((((c)#1,b))#2,d)))),#3) cost=1
   >: ((#2,((((((c)#1,b))#2,a))#3,(#1,d))),#3) cost=0
   i: (#3,((b)#1,((((c)#3,(d,#1)))#2,(a,#2)))) cost=2
   >: (#3,((b)#1,((((d,(c)#3),#1))#2,(a,#2)))) cost=1
   >: (#3,(((b)#1,((((d,(c)#3),#1))#2,a)),#2)) cost=0
   i: (#1,((#3,((d)#2,((((a)#3,c))#1,b))),#2)) cost=2
   >: (#1,(#3,((d)#2,((((a)#3,(c,#2)))#1,b)))) cost=0
   i: ((((((((d)#3,a))#2,(#2,c)))#1,b),#1),#3) cost=4
   >: (((#2,((((((d)#3,a))#2,c))#1,b)),#1),#3) cost=3
   >: ((#2,((((((d)#3,a))#2,(c,#3)))#1,b)),#1) cost=0
   i: ((#1,(#2,((((((c)#2,a))#1,d))#3,b))),#3) cost=3
   >: ((#1,((((((c)#2,a))#1,(#2,d)))#3,b)),#3) cost=0
Cost:0 TopNetworks:254 Steps:470980 Climbs:2623 HCruns:1000 HCTime:7.72257 Naive:471980 Naivetime:0.81058 DTcnt:3687318 Class:TreeChild TimeConsistency:0
```

### Using stopping criterion (--hcstopinit=1000) with quasi-consensus rand networks.

```  
> supnet -g "(a,(b,(c,d))); ((a,b),(c,d))" -q-1 -R3 --HC --hcrunstats --hcstopinit=1000 --randseed 13 | tail -10
   i: ((d)#2,((((b,#2))#3,(((#3,c))#1,a)),#1)) cost=2
   >: ((d)#2,(((b)#3,(((#3,(c,#2)))#1,a)),#1)) cost=0
   i: ((#3,((d)#2,(((((a)#3,b),#2))#1,c))),#1) cost=1
   >: (#3,(#1,((d)#2,(((((a)#3,b),#2))#1,c)))) cost=0
   i: ((d)#2,(#3,(#2,((c)#3,((b)#1,(a,#1)))))) cost=4
   >: ((d)#2,(#3,((#2,((c)#3,((b)#1,a))),#1))) cost=3
   >: ((d)#2,(#3,((((#2,c))#3,((b)#1,a)),#1))) cost=1
   >: ((d)#2,(#3,((((#2,c),#1))#3,((b)#1,a)))) cost=0
   i: ((#1,(((((c)#2,((b)#1,a)))#3,d),#2)),#3) cost=0
Cost:0 TopNetworks:254 Steps:263275 Climbs:1449 HCruns:547 HCTime:4.39665 Naive:263822 Naivetime:0.458773 DTcnt:2060491 Class:TreeChild TimeConsistency:0
```

### Using time stop criterion in HC in seconds `--hcstoptime=TIME`, i.e., stop a HC climb if there is no improvement after TIME seconds.

```
> supnet -g "((i,c),((a,d),(f,(b,((g,(e,j)),h)))));((i,h),((c,f),((a,(d,e)),(j,(g,b)))));(((f,(b,d)),((j,g),(e,i))),(c,(h,a)));(((f,(i,(j,d))),(c,b)),(((h,a),g),e));((((h,e),((c,f),a)),(d,b)),((g,i),j))" -R8 -q1 -v4 --HC --hcrunstats --hcstoptime=0.5 --randseed 13
HC start: hcusenaive=0 runnaiveleqrt=13 tailmove=1
   i: ((((((d)#3,((((#3,((a)#7,h)))#6,c),#8)))#5,e))#2,((((((i)#4,(#7,(#2,g))))#1,(b,(#1,(#4,((j)#8,f))))),#5),#6)) cost=49
   >: ((((((d)#3,((((#3,((a)#7,h)))#6,c),#8)))#5,e))#2,(((((((#2,i))#4,(#7,g)))#1,(b,(#1,(#4,((j)#8,f))))),#5),#6)) cost=48
   >: ((((((d)#3,((((#3,((a)#7,h)))#6,c),#8)))#5,(#7,e)))#2,(((((((#2,i))#4,g))#1,(b,(#1,(#4,((j)#8,f))))),#5),#6)) cost=47
   >: ((((((d)#3,((((#3,((a)#7,h)))#6,c),#8)))#5,(#7,e)))#2,(((#4,(((((#2,i))#4,g))#1,(b,(#1,((j)#8,f))))),#5),#6)) cost=46
   >: ((((((d)#3,((((#3,((a)#7,h)))#6,c),#8)))#5,(#7,e)))#2,(((#4,(((((#2,i))#4,g))#1,((#1,b),((j)#8,f)))),#5),#6)) cost=45
   >: ((((((d)#3,(((((#3,b),((a)#7,h)))#6,c),#8)))#5,(#7,e)))#2,(((#4,(((((#2,i))#4,g))#1,(#1,((j)#8,f)))),#5),#6)) cost=43
   >: ((((((d)#3,(((((#3,b),((a)#7,h)))#6,c),#8)))#5,(#7,e)))#2,((#1,((#4,(((((#2,i))#4,g))#1,((j)#8,f))),#5)),#6)) cost=42
   >: (((((((d,b))#3,((((#3,((a)#7,h)))#6,c),#8)))#5,(#7,e)))#2,((#1,((#4,(((((#2,i))#4,g))#1,((j)#8,f))),#5)),#6)) cost=41
   >: (((((((d,b))#3,((((#3,((a)#7,h)))#6,c),#8)))#5,(#7,e)))#2,((#1,((#4,((((i)#4,(#2,g)))#1,((j)#8,f))),#5)),#6)) cost=39
   >: (((((((d,b))#3,(((#3,((a)#7,h)))#6,c)))#5,(#7,e)))#2,((#1,((#4,((((i)#4,(#2,(g,#8))))#1,((j)#8,f))),#5)),#6)) cost=35
   >: (((((((d,b))#3,((((a)#7,h))#6,c)))#5,(#7,(#3,e))))#2,((#1,((#4,((((i)#4,(#2,(g,#8))))#1,((j)#8,f))),#5)),#6)) cost=33
   >: ((((((((#7,d),b))#3,((((a)#7,h))#6,c)))#5,(#3,e)))#2,((#1,((#4,((((i)#4,(#2,(g,#8))))#1,((j)#8,f))),#5)),#6)) cost=32
   >: ((((((((#7,d),b))#3,((((a)#7,h))#6,c)))#5,(#3,e)))#2,((#1,(#4,(((((i)#4,(#2,(g,#8))))#1,((j)#8,f)),#5))),#6)) cost=31
   >: ((((((((#7,d),b))#3,((((a)#7,h))#6,(#4,c))))#5,(#3,e)))#2,((#1,(((((i)#4,(#2,(g,#8))))#1,((j)#8,f)),#5)),#6)) cost=30
   >: ((((((((#7,d),b))#3,((((a)#7,h))#6,(#4,c))))#5,(#3,(e,#6))))#2,(#1,(((((i)#4,(#2,(g,#8))))#1,((j)#8,f)),#5))) cost=29
   >: (((((((d,b))#3,((((a)#7,h))#6,(#4,c))))#5,(#3,(e,#6))))#2,(#1,((((#7,((i)#4,(#2,(g,#8)))))#1,((j)#8,f)),#5))) cost=28
   >: (((((((d,b))#3,((((a)#7,h))#6,(#4,c))))#5,(e,#6)))#2,(#1,((((#7,(#3,((i)#4,(#2,(g,#8))))))#1,((j)#8,f)),#5))) cost=27
   >: ((((((((d,b))#3,((((a)#7,h))#6,(#4,c))))#5,(e,#6)))#2,((((#7,(#3,((i)#4,(#2,(g,#8))))))#1,((j)#8,f)),#5)),#1) cost=25
   >: (((((((((a)#7,h))#6,(#4,c)))#5,(e,#6)))#2,((((#7,(#3,((i)#4,(#2,(g,#8))))))#1,(((d,b))#3,((j)#8,f))),#5)),#1) cost=24
   >: (((((((((a)#7,h))#6,(#4,c)))#5,(e,#6)))#2,((((#7,(#3,((i)#4,(#2,(g,#8))))))#1,(((((j)#8,d),b))#3,f)),#5)),#1) cost=23
   >: ((((e,#6))#2,((((#7,(#3,((i)#4,(#2,(g,#8))))))#1,(((((((((a)#7,h))#6,(#4,c)))#5,((j)#8,d)),b))#3,f)),#5)),#1) cost=22
   >: (((((#7,((((((a)#7,h))#6,(#4,c)))#5,((j)#8,(((e,#6))#2,d)))),b))#3,((((#3,(#2,((i)#4,(g,#8)))))#1,f),#5)),#1) cost=20
   >: (((((#7,(((#4,c))#5,((j)#8,d))),b))#3,((((((e)#2,a))#7,h))#6,((((#3,((#2,(((i)#4,g),#8)),#6)))#1,f),#5))),#1) cost=19
HC run stopped due to timeout, cost=19
1.  Steps:0 Climbs:0 New optimal cost: 19
Stats saved to odt.dat
Best networks saved to odt.log
Cost:19 TopNetworks:68 Steps:8367 Climbs:46 HCruns:1 HCTime:1.83863 Naive:8368 Naivetime:1.53764 DTcnt:2142208 Class:TreeChild TimeConsistency:0
```

### Save odt files after each improvement

Use `--hcsavewhenimproved`. Note that stats may be not reliable, before the end of computation.

```
> supnet -g "((i,c),((a,d),(f,(b,((g,(e,j)),h)))));((i,h),((c,f),((a,(d,e)),(j,(g,b)))));(((f,(b,d)),((j,g),(e,i))),(c,(h,a)));(((f,(i,(j,d))),(c,b)),(((h,a),g),e));((((h,e),((c,f),a)),(d,b)),((g,i),j))" -R8 -q1 -v4 --HC --hcrunstats --hcstoptime=0.5 --randseed 13 --hcsavewhenimproved
HC start: hcusenaive=0 runnaiveleqrt=13 tailmove=1
   i: ((((((d)#3,((((#3,((a)#7,h)))#6,c),#8)))#5,e))#2,((((((i)#4,(#7,(#2,g))))#1,(b,(#1,(#4,((j)#8,f))))),#5),#6)) cost=49
   >: ((((((d)#3,((((#3,((a)#7,h)))#6,c),#8)))#5,e))#2,(((((((#2,i))#4,(#7,g)))#1,(b,(#1,(#4,((j)#8,f))))),#5),#6)) cost=48
   >: ((((((d)#3,((((#3,((a)#7,h)))#6,c),#8)))#5,(#7,e)))#2,(((((((#2,i))#4,g))#1,(b,(#1,(#4,((j)#8,f))))),#5),#6)) cost=47
   >: ((((((d)#3,((((#3,((a)#7,h)))#6,c),#8)))#5,(#7,e)))#2,(((#4,(((((#2,i))#4,g))#1,(b,(#1,((j)#8,f))))),#5),#6)) cost=46
   >: ((((((d)#3,((((#3,((a)#7,h)))#6,c),#8)))#5,(#7,e)))#2,(((#4,(((((#2,i))#4,g))#1,((#1,b),((j)#8,f)))),#5),#6)) cost=45
   >: ((((((d)#3,(((((#3,b),((a)#7,h)))#6,c),#8)))#5,(#7,e)))#2,(((#4,(((((#2,i))#4,g))#1,(#1,((j)#8,f)))),#5),#6)) cost=43
   >: ((((((d)#3,(((((#3,b),((a)#7,h)))#6,c),#8)))#5,(#7,e)))#2,((#1,((#4,(((((#2,i))#4,g))#1,((j)#8,f))),#5)),#6)) cost=42
   >: (((((((d,b))#3,((((#3,((a)#7,h)))#6,c),#8)))#5,(#7,e)))#2,((#1,((#4,(((((#2,i))#4,g))#1,((j)#8,f))),#5)),#6)) cost=41
   >: (((((((d,b))#3,((((#3,((a)#7,h)))#6,c),#8)))#5,(#7,e)))#2,((#1,((#4,((((i)#4,(#2,g)))#1,((j)#8,f))),#5)),#6)) cost=39
   >: (((((((d,b))#3,(((#3,((a)#7,h)))#6,c)))#5,(#7,e)))#2,((#1,((#4,((((i)#4,(#2,(g,#8))))#1,((j)#8,f))),#5)),#6)) cost=35
   >: (((((((d,b))#3,((((a)#7,h))#6,c)))#5,(#7,(#3,e))))#2,((#1,((#4,((((i)#4,(#2,(g,#8))))#1,((j)#8,f))),#5)),#6)) cost=33
   >: ((((((((#7,d),b))#3,((((a)#7,h))#6,c)))#5,(#3,e)))#2,((#1,((#4,((((i)#4,(#2,(g,#8))))#1,((j)#8,f))),#5)),#6)) cost=32
   >: ((((((((#7,d),b))#3,((((a)#7,h))#6,c)))#5,(#3,e)))#2,((#1,(#4,(((((i)#4,(#2,(g,#8))))#1,((j)#8,f)),#5))),#6)) cost=31
   >: ((((((((#7,d),b))#3,((((a)#7,h))#6,(#4,c))))#5,(#3,e)))#2,((#1,(((((i)#4,(#2,(g,#8))))#1,((j)#8,f)),#5)),#6)) cost=30
   >: ((((((((#7,d),b))#3,((((a)#7,h))#6,(#4,c))))#5,(#3,(e,#6))))#2,(#1,(((((i)#4,(#2,(g,#8))))#1,((j)#8,f)),#5))) cost=29
   >: (((((((d,b))#3,((((a)#7,h))#6,(#4,c))))#5,(#3,(e,#6))))#2,(#1,((((#7,((i)#4,(#2,(g,#8)))))#1,((j)#8,f)),#5))) cost=28
   >: (((((((d,b))#3,((((a)#7,h))#6,(#4,c))))#5,(e,#6)))#2,(#1,((((#7,(#3,((i)#4,(#2,(g,#8))))))#1,((j)#8,f)),#5))) cost=27
   >: ((((((((d,b))#3,((((a)#7,h))#6,(#4,c))))#5,(e,#6)))#2,((((#7,(#3,((i)#4,(#2,(g,#8))))))#1,((j)#8,f)),#5)),#1) cost=25
   >: (((((((((a)#7,h))#6,(#4,c)))#5,(e,#6)))#2,((((#7,(#3,((i)#4,(#2,(g,#8))))))#1,(((d,b))#3,((j)#8,f))),#5)),#1) cost=24
   >: (((((((((a)#7,h))#6,(#4,c)))#5,(e,#6)))#2,((((#7,(#3,((i)#4,(#2,(g,#8))))))#1,(((((j)#8,d),b))#3,f)),#5)),#1) cost=23
   >: ((((e,#6))#2,((((#7,(#3,((i)#4,(#2,(g,#8))))))#1,(((((((((a)#7,h))#6,(#4,c)))#5,((j)#8,d)),b))#3,f)),#5)),#1) cost=22
   >: (((((#7,((((((a)#7,h))#6,(#4,c)))#5,((j)#8,(((e,#6))#2,d)))),b))#3,((((#3,(#2,((i)#4,(g,#8)))))#1,f),#5)),#1) cost=20
HC run stopped due to timeout, cost=20
1.  Steps:0 Climbs:0 New optimal cost: 20
Stats saved to odt.dat
Best networks saved to odt.log
Cost:20 TopNetworks:48 Steps:6056 Climbs:44 HCruns:1 HCTime:1.43186 Naive:6057 Naivetime:1.20076 DTcnt:1550592 Class:TreeChild TimeConsistency:0
```

### Odt files optimizations

Set '--hcignorecostgeq=COST' to omit saving odt files if the cost is too large.
Set '--hcsavefinalodt' to save the final odt/log files, even if `--hcignorecostgeq` is set.


## Display tree sampler in HC

HC algorithm may start with the approximated cost, before switching to exact computations, by sampling display trees in the naive algorithm.

To enable random sampling for display trees, you can utilize the `--displaytreesampling` flag with the following format: `--displaytreesampling="PAR-1 PAR-2 ... PAR-N"`. Here, the parameters are specified as a space-separated string, where each 'PAR' represents a parameter of an exponential distribution.
In such a scenario, the HC algorithm will engage in a sequence of N+1 iterative climbs, first utilizing the provided parameters. The optimal networks obtained from one sampling parameter will serve as the starting point for the subsequent climb with the next parameter. The final (N+1)-th climb is performed using the exact algorithm. Setting `PAR=0.5` will result in approximately half of the display trees being sampled 
(due to the integer sampler implementation the number of sampled trees for parameters >=1 will be smaller than expected). Setting the value `PAR=0.25` will yield roughly a quater of display trees and so on. Parameters should be monotonically increasing.

```
> supnet -g "((i,c),((a,d),(f,(b,((g,(e,j)),h)))));((i,h),((c,f),((a,(d,e)),(j,(g,b)))));(((f,(b,d)),((j,g),(e,i))),(c,(h,a)));(((f,(i,(j,d))),(c,b)),(((h,a),g),e));((((h,e),((c,f),a)),(d,b)),((g,i),j))" -R8 -q1 -v4 --HC --hcrunstats --hcstoptime=0.5 --randseed 13 --hcsavewhenimproved --displaytreesampling="0.125 0.25 0.5"
HC start: hcusenaive=0 runnaiveleqrt=13 tailmove=1
   i: ((((((d)#3,((((#3,((a)#7,h)))#6,c),#8)))#5,e))#2,((((((i)#4,(#7,(#2,g))))#1,(b,(#1,(#4,((j)#8,f))))),#5),#6)) cost=54
   >: (((((((#7,d))#3,((((#3,((a)#7,h)))#6,c),#8)))#5,e))#2,((((((i)#4,(#2,g)))#1,(b,(#1,(#4,((j)#8,f))))),#5),#6)) cost=51
   >: (((((((#7,d))#3,((((#3,((a)#7,(#4,h))))#6,c),#8)))#5,e))#2,((((((i)#4,(#2,g)))#1,(b,(#1,((j)#8,f)))),#5),#6)) cost=50
   >: (((((((#7,d))#3,((((#3,((a)#7,(#4,h))))#6,c),#8)))#5,e))#2,(((((((i)#4,g),#2))#1,(b,(#1,((j)#8,f)))),#5),#6)) cost=49
   >: (((((((#7,d))#3,((((#3,((a)#7,(#4,h))))#6,c),#8)))#5,e))#2,((#1,((((((i)#4,g),#2))#1,(b,((j)#8,f))),#5)),#6)) cost=46
   >: (((((((#7,(b,d)))#3,((((#3,((a)#7,(#4,h))))#6,c),#8)))#5,e))#2,((#1,((((((i)#4,g),#2))#1,((j)#8,f)),#5)),#6)) cost=45
   >: (((((((b,d))#3,((((#3,((a)#7,(#4,h))))#6,c),#8)))#5,e))#2,((#1,((#7,(((((i)#4,g),#2))#1,((j)#8,f))),#5)),#6)) cost=42
   >: (((((((#3,((a)#7,(#4,h))))#6,c))#5,(((b,d))#3,e)))#2,((#1,((#7,(((((i)#4,(g,#8)),#2))#1,((j)#8,f))),#5)),#6)) cost=41
   >: (((((((#3,((a)#7,(#4,h))))#6,c))#5,(((b,d))#3,e)))#2,((#1,((#7,((((((i)#4,g),#8),#2))#1,((j)#8,f))),#5)),#6)) cost=39
   >: ((((((((a)#7,(#4,h)))#6,c))#5,(((b,d))#3,e)))#2,((#1,((#7,(((#3,((((i)#4,g),#8),#2)))#1,((j)#8,f))),#5)),#6)) cost=38
   >: ((((((((a)#7,(#4,h)))#6,c))#5,(((b,d))#3,e)))#2,((#1,((((#3,(#7,((((i)#4,g),#8),#2))))#1,((j)#8,f)),#5)),#6)) cost=37
   >: ((((((((a)#7,(#4,h)))#6,c))#5,(((b,d))#3,e)))#2,((#1,((((#3,(#7,(((i)#4,(g,#8)),#2))))#1,((j)#8,f)),#5)),#6)) cost=35
   >: ((((((((a)#7,(#4,h)))#6,c))#5,(((#7,(b,d)))#3,e)))#2,((#1,((((#3,(((i)#4,(g,#8)),#2)))#1,((j)#8,f)),#5)),#6)) cost=33
   >: (((((#7,(b,d)))#3,e))#2,((#1,((((((a)#7,(#4,h)))#6,c))#5,((((#3,(((j)#8,((i)#4,(g,#8))),#2)))#1,f),#5))),#6)) cost=31
   >: ((((((#7,d),b))#3,e))#2,((#1,((((((a)#7,(#4,h)))#6,c))#5,((((#3,(((j)#8,((i)#4,(g,#8))),#2)))#1,f),#5))),#6)) cost=30
   >: (((((((#7,d),#8),b))#3,e))#2,((#1,((c)#5,((((a)#7,(#4,h)))#6,(((#3,(((j)#8,((i)#4,g)),#2)))#1,(#5,f))))),#6)) cost=29
HC run with DTsampling=0.125 completed: (((((((#7,d),#8),b))#3,e))#2,((#1,((c)#5,((((a)#7,(#4,h)))#6,(((#3,(((j)#8,((i)#4,g)),#2)))#1,(#5,f))))),#6)) cost=29
1.  Steps:0 Climbs:0 New optimal cost: 29
...
   i: ((#4,((((c)#5,((a)#7,(((((#7,(d,#8)))#3,e))#2,h))))#6,((#3,(#1,b)),(#5,f)))),(((i)#4,((((j)#8,g),#2),#6)))#1) cost=19
   >: ((#4,((((c)#5,((a)#7,((e)#2,h))))#6,((#3,(#1,b)),(#5,f)))),((((((#7,(d,#8)))#3,i))#4,((((j)#8,g),#2),#6)))#1) cost=18
   >: ((#4,((((c)#5,((a)#7,h)))#6,((#3,(#1,b)),(#5,f)))),((((((#7,((e)#2,(d,#8))))#3,i))#4,((((j)#8,g),#2),#6)))#1) cost=17
   >: ((#4,((((c)#5,((((e)#2,a))#7,h)))#6,((#3,(#1,b)),(#5,f)))),((((((#7,(d,#8)))#3,i))#4,((((j)#8,g),#2),#6)))#1) cost=16
HC run stopped due to timeout, cost=16
1.  Steps:0 Climbs:0 New optimal cost: 16
   i: ((c)#5,((#4,(((#7,(d,#8)))#3,((((a)#7,((e)#2,h)))#6,((#3,(#1,b)),(#5,f))))),(((i)#4,((((j)#8,g),#2),#6)))#1)) cost=19
   >: ((c)#5,((#4,(((#7,((e)#2,(d,#8))))#3,((((a)#7,h))#6,((#3,(#1,b)),(#5,f))))),(((i)#4,((((j)#8,g),#2),#6)))#1)) cost=18
   >: ((c)#5,((#4,(((#7,(d,#8)))#3,((((((e)#2,a))#7,h))#6,((#3,(#1,b)),(#5,f))))),(((i)#4,((((j)#8,g),#2),#6)))#1)) cost=17
   >: ((c)#5,((#4,((((((e)#2,a))#7,h))#6,((#3,(#1,b)),(#5,f)))),((((((#7,(d,#8)))#3,i))#4,((((j)#8,g),#2),#6)))#1)) cost=16
HC run stopped due to timeout, cost=16
   i: ((c)#5,((#4,((((a)#7,(((((#7,(d,#8)))#3,e))#2,h)))#6,((#3,(#1,b)),(#5,f)))),(((i)#4,((((j)#8,g),#2),#6)))#1)) cost=19
   >: ((c)#5,((#4,((((a)#7,((e)#2,h)))#6,((#3,(#1,b)),(#5,f)))),((((((#7,(d,#8)))#3,i))#4,((((j)#8,g),#2),#6)))#1)) cost=18
   >: ((c)#5,((#4,((((a)#7,h))#6,((#3,(#1,b)),(#5,f)))),((((((#7,((e)#2,(d,#8))))#3,i))#4,((((j)#8,g),#2),#6)))#1)) cost=17
   >: ((c)#5,((#4,((((((e)#2,a))#7,h))#6,((#3,(#1,b)),(#5,f)))),((((((#7,(d,#8)))#3,i))#4,((((j)#8,g),#2),#6)))#1)) cost=16
HC run stopped due to timeout, cost=16
Stats data saved: odt.dat
Best networks saved to odt.log
Cost:16 TopNetworks:59 Steps:19327 Climbs:78 HCruns:6 HCTime:2.94697 Naive:19333 Naivetime:2.54067 DTcnt:3063305 Class:TreeChild TimeConsistency:0
```

Use `hcsamplingmaxnetstonextlevel=MAXNETS` to limit the number of networks passed to the next level sampler. Default is 0 (unlimited).

To test the effectiveness of display tree sampling employ 
the `--testdisplaytreesampling=RETICULATIONCNT` flag, where 'RETICULATIONCNT' indicates the number of reticulations. Here the number of reticulation was set to 10, i.e., to 1024 display trees. 

```
> supnet --HC --displaytreesampling="0.0625 0.125 0.25 0.5 1" --testdisplaytreesampling=10
 Sampler 0.0625 79
 Sampler 0.125 124
 Sampler 0.25 256
 Sampler 0.5 531
 Sampler 1 834
 Sampler 0 1024
```

When sampling if a display tree already yields a better cost that found in HC, the climb can continue immediately if `--cutwhendtimproved` is set. To be tested.

Use `--hcdetailedsummary` to print summary on each sampling level.
Use `--hcdetailedsummarydat` to write additional summary dat file.


## Automatic names of outfiles

Configure the `--autooutfiles` option to generate filenames with the cost value followed by `dat` or `log` In the event of naming conflicts, increment numbers will be added, such as COST.1.dat, COST.2.dat, and so forth. Try with '--outdirectory'.

## More examples

To print the minimum total cost for 10 random gene trees versus random tree-child networks with 5 reticulations over 8 species, and also print the initial network:

```
> supnet -r10 -A8 --pnetworks  | supnet -G- -r1 -A8 -R5 --pnetworks --odtnaivecost
((((((e,(((f)#3,d),#1)))#5,a))#2,((#4,((c)#1,(((#2,g))#4,h))),#3)),(#5,b))
64
```


## Coronavirus dataset processing

Infering tree-child networks without time-consistency, start from 10 networks `-q10`.
Use: `--hcsavewhenimproved` to save odt files after each improvement
and  `--hcstoptime=30` to stop a single HC if there is no improvement after 30 seconds.

```
supnet -G corona.txt -R7 -q10 --HC --hcrunstats --outfiles corona --hcstoptime=30  --hcsavewhenimproved
```
Results in `corona.log` and `corona.dat`.

Inferring relaxed time-consistent networks, start from 10 networks `-q10`.
```
supnet -G corona.txt -R7 -q10 --HC --hcrunstats --timeconsistent --relaxed --outdirectory corona_tc_rel
```

<!-- Inferring relaxed time-consistent networks, start from 10 networks `-q10`, save networks to #first.log.
```
supnet -G corona.txt -R7 -q10 --HC -v1 --timeconsistent --relaxed --outfiles first --hcsavewhenimproved
supnet -G corona.txt -N first.log --hcsavewhenimproved  
supnet -G corona.txt -N first.log --HC -v1 --timeconsistent --relaxed --hcsavewhenimproved
``` -->

### Parallel processing 

TODO: multithreaded supnet

Use parallel with 10 jobs to run 20 supnet computations.
```
parallel --jobs 10 --ungroup supnet -G corona.txt -R7 -q1 --HC --hcrunstats --outdirectory corona --hcstoptime 1 --hcsavewhenimproved --odtlabelled --autooutfiles ::: {1..20}
```

Extracting best networks from log files using supnet:
```
cat corona/*.log | supnet -G corona.txt -N- --bestnetworks -v4
cat odt.log
```

Download `csvmanip` to merge dat files:
```
git clone git@github.com:ppgorecki/csvmanip.git 
```

Use `csvmanip` to merge dat files into csv:
```
> csvmanip/csvmanip.py corona/*.dat
Id,Source,optcost,time,hctime,mergetime,topnets,class,timeconsistency,improvements,steps,bbruns,startnets,memoryMB,dtcnt,naivetime,naivecnt,bbnaivecnt,bbnaivetime,bbdpcnt,bbdptime,randseed
119,corona/119.dat,119,5.2893,5.2893,2.38419e-07,2,0,0,80,10186,0,1,261,1303936,4.52791,10187,0,0,0,0,3894346635
122,corona/122.dat,122,3.55873,3.55873,4.76837e-07,3,0,0,86,8837,0,1,258,1131264,2.75365,8838,0,0,0,0,3894353523
24.1,corona/24.1.dat,24,3.38207,3.38207,4.76837e-07,69,0,0,102,26871,0,2,125,1491739,2.34371,26873,0,0,0,0,3892887752
24.2,corona/24.2.dat,24,7.1368,7.1368,2.38419e-07,64,0,0,128,26560,0,1,261,3399808,6.06366,26561,0,0,0,0,3894353403
24,corona/24.dat,24,6.76618,6.76618,2.38419e-07,2,0,0,110,25438,0,1,262,3256192,5.81489,25439,0,0,0,0,3894346632
25.1,corona/25.1.dat,25,7.52293,7.52293,2.38419e-07,4,0,0,114,24546,0,1,260,3142016,6.63867,24547,0,0,0,0,3892752804
25.2,corona/25.2.dat,25,4.80912,4.80912,2.38419e-07,8,0,0,136,18157,0,1,196,2324224,4.55608,18158,0,0,0,0,3892756406
25.3,corona/25.3.dat,25,5.77002,5.77002,2.38419e-07,78,0,0,94,20104,0,1,235,2573440,5.18477,20105,0,0,0,0,3892757807
25.4,corona/25.4.dat,25,3.02409,3.02409,2.38419e-07,85,0,0,122,22827,0,2,80,1091661,2.17269,22829,0,0,0,0,3892889004
25.5,corona/25.5.dat,25,4.71099,4.71099,9.53674e-07,20,0,0,92,16743,0,1,202,2143232,4.47507,16744,0,0,0,0,3894346634
25.6,corona/25.6.dat,25,4.99702,4.99702,9.53674e-07,17,0,0,76,23172,0,1,218,2966144,4.50657,23173,0,0,0,0,3894346638
25.7,corona/25.7.dat,25,8.88987,8.88987,2.38419e-07,24,0,0,112,29396,0,1,261,3762816,7.75601,29397,0,0,0,0,3894346637
25.8,corona/25.8.dat,25,5.80147,5.80147,2.38419e-07,2,0,0,102,20142,0,1,262,2578304,4.6518,20143,0,0,0,0,3894351422
25,corona/25.dat,25,7.33998,7.33998,2.38419e-07,103,0,0,128,23076,0,1,261,2953856,6.16913,23077,0,0,0,0,3892752802
26.1,corona/26.1.dat,26,6.26758,6.26758,2.38419e-07,100,0,0,96,17760,0,1,261,2273408,5.37885,17761,0,0,0,0,3892757353
26.2,corona/26.2.dat,26,3.72703,3.72703,4.76837e-07,64,0,0,84,21394,0,2,121,1296013,2.84425,21396,0,0,0,0,3892885079
26.3,corona/26.3.dat,26,3.7391,3.7391,2.38419e-07,4,0,0,104,18221,0,2,168,2025251,3.17031,18223,0,0,0,0,3892885074
26.4,corona/26.4.dat,26,2.43439,2.43438,7.15256e-07,80,0,0,108,24853,0,2,82,1280504,1.86578,24855,0,0,0,0,3892889650
...
33.2,corona/33.2.dat,33,6.70923,6.70923,2.38419e-07,35,0,0,114,19246,0,1,262,2463616,5.95569,19247,0,0,0,0,3894346643
33.3,corona/33.3.dat,33,3.27832,3.27832,2.38419e-07,126,0,0,120,22867,0,2,176,1864966,2.74244,22869,0,0,0,0,3892892150
33.4,corona/33.4.dat,33,5.04684,5.04684,2.38419e-07,44,0,0,86,14390,0,1,262,1842048,4.20835,14391,0,0,0,0,3894351965
33,corona/33.dat,33,5.15879,5.15879,2.38419e-07,2,0,0,98,15065,0,1,261,1928448,4.38866,15066,0,0,0,0,3892757202
34.1,corona/34.1.dat,34,12.7193,12.7193,1.19209e-06,98,0,0,194,49735,0,4,263,5835263,10.8442,49739,0,0,0,0,3892891669
34.2,corona/34.2.dat,34,0,0,0,0,0,0,46,0,0,0,197,0,0,0,0,0,0,0,3894341672
34,corona/34.dat,34,1.84298,1.84298,0,85,0,0,90,19013,0,2,72,1109007,1.40869,19015,0,0,0,0,3892892152
35,corona/35.dat,35,2.4352,2.4352,2.38419e-07,174,0,0,154,55996,0,6,124,912761,1.39763,56002,0,0,0,0,3892867240
36,corona/36.dat,36,4.74186,4.74186,4.76837e-07,171,0,0,122,31319,0,3,144,2110112,3.48464,31322,0,0,0,0,3892890750
39,corona/39.dat,39,0,0,0,0,0,0,45,0,0,0,183,0,0,0,0,0,0,0,3894341681
41.1,corona/41.1.dat,41,5.14744,5.14744,2.38419e-07,2,0,0,110,15571,0,1,259,1993216,4.47464,15572,0,0,0,0,3892758128
42,corona/42.dat,42,2.14001,2.14,2.38419e-07,131,0,0,72,22148,0,2,87,1090153,1.52107,22150,0,0,0,0,3892891265
46,corona/46.dat,46,4.49964,4.49964,2.38419e-07,5,0,0,82,13012,0,1,260,1665664,3.71095,13013,0,0,0,0,3892752808
50.1,corona/50.1.dat,50,5.30467,5.30467,4.76837e-07,8,0,0,102,14107,0,1,262,1805824,4.51396,14108,0,0,0,0,3892752801
55,corona/55.dat,55,5.05534,5.05534,2.38419e-07,5,0,0,98,15144,0,1,261,1938560,4.1888,15145,0,0,0,0,3894346636
56,corona/56.dat,56,4.56794,4.56794,4.76837e-07,6,0,0,108,9637,0,1,261,1233664,3.69766,9638,0,0,0,0,3894346639
57,corona/57.dat,57,4.76359,4.76359,9.53674e-07,6,0,0,82,10013,0,1,260,1281792,3.95678,10014,0,0,0,0,3894346640
62,corona/62.dat,62,0,0,0,0,0,0,44,0,0,0,217,0,0,0,0,0,0,0,3894341677
96,corona/96.dat,96,0,0,0,0,0,0,37,0,0,0,172,0,0,0,0,0,0,0,3894341673
```

Print cost + networks with sort:
```
> cat corona/*.log | supnet -N- -G corona.txt --odtcost | sort -k1 -n
24 (#1,(#5,(#6,((((GuangxiPangolinP2V)#7,((GuangdongPangolin12019,#7),(((HuItalyTE48362020,HuWuhan2020),RaTG13),((BatCoVZC45,(BM4831BGR2008)#5),BatCoVZXC21)))))#1,(((HKU312,#4),((((BtCoV2792005)#4,(#2,BtCoV2732005)))#6,((((Rf1)#2,Rs3367))#3,(SARSCoVBJ1824,SARS)))),#3)))))
24 (#1,(#5,(#6,((((GuangxiPangolinP2V)#7,((GuangdongPangolin12019,#7),(((HuItalyTE48362020,HuWuhan2020),RaTG13),((BatCoVZC45,(BM4831BGR2008)#5),BatCoVZXC21)))))#1,(((HKU312,#4),((((BtCoV2792005)#4,(#2,BtCoV2732005)))#6,((Rf1)#2,((Rs3367)#3,(SARSCoVBJ1824,SARS))))),#3)))))
24 (#1,(#5,(#6,((((GuangxiPangolinP2V)#7,((GuangdongPangolin12019,#7),(((HuItalyTE48362020,HuWuhan2020),RaTG13),((BatCoVZC45,(BM4831BGR2008)#5),BatCoVZXC21)))))#1,(((HKU312,#4),((((BtCoV2792005)#4,(#2,BtCoV2732005)))#6,((Rs3367)#3,(((Rf1)#2,SARSCoVBJ1824),SARS)))),#3)))))
24 (#1,(#5,(#6,((((GuangxiPangolinP2V)#7,((GuangdongPangolin12019,#7),(((HuItalyTE48362020,HuWuhan2020),RaTG13),((BatCoVZC45,(BM4831BGR2008)#5),BatCoVZXC21)))))#1,(((HKU312,#4),((((BtCoV2792005)#4,(#2,BtCoV2732005)))#6,((Rs3367)#3,((Rf1)#2,(SARSCoVBJ1824,SARS))))),#3)))))
24 (#1,(#5,(#6,((((GuangxiPangolinP2V)#7,((GuangdongPangolin12019,#7),(((HuItalyTE48362020,HuWuhan2020),RaTG13),((BatCoVZC45,(BM4831BGR2008)#5),BatCoVZXC21)))))#1,(((HKU312,#4),((((BtCoV2792005)#4,(#2,BtCoV2732005)))#6,((Rs3367)#3,(SARSCoVBJ1824,((Rf1)#2,SARS))))),#3)))))
24 (#1,(#5,(#6,((((GuangxiPangolinP2V)#7,((GuangdongPangolin12019,#7),(((HuItalyTE48362020,HuWuhan2020),RaTG13),((BatCoVZC45,(BM4831BGR2008)#5),BatCoVZXC21)))))#1,(((HKU312,#4),((((BtCoV2792005)#4,(#2,((Rf1)#2,BtCoV2732005))))#6,((Rs3367)#3,(SARSCoVBJ1824,SARS)))),#3)))))
24 (#1,(#5,(#6,((((GuangxiPangolinP2V)#7,((GuangdongPangolin12019,#7),(((HuItalyTE48362020,HuWuhan2020),RaTG13),((BatCoVZC45,(BM4831BGR2008)#5),BatCoVZXC21)))))#1,(((HKU312,#4),((((((Rf1)#2,BtCoV2792005))#4,(#2,BtCoV2732005)))#6,((Rs3367)#3,(SARSCoVBJ1824,SARS)))),#3)))))
24 (#1,(#5,(#6,((((GuangxiPangolinP2V)#7,((GuangdongPangolin12019,#7),(((HuItalyTE48362020,HuWuhan2020),RaTG13),((BatCoVZC45,(BM4831BGR2008)#5),BatCoVZXC21)))))#1,(((HKU312,#4),((((Rf1)#2,((BtCoV2792005)#4,(#2,BtCoV2732005))))#6,((Rs3367)#3,(SARSCoVBJ1824,SARS)))),#3)))))
24 (#1,(#5,(#6,((((GuangxiPangolinP2V)#7,((GuangdongPangolin12019,#7),(((HuItalyTE48362020,HuWuhan2020),RaTG13),((BatCoVZC45,(BM4831BGR2008)#5),BatCoVZXC21)))))#1,(((HKU312,#4),((Rf1)#2,((((BtCoV2792005)#4,(#2,BtCoV2732005)))#6,((Rs3367)#3,(SARSCoVBJ1824,SARS))))),#3)))))
24 (#1,(#5,(#6,((((GuangxiPangolinP2V)#7,((GuangdongPangolin12019,#7),(((HuItalyTE48362020,HuWuhan2020),RaTG13),((BatCoVZC45,(BM4831BGR2008)#5),BatCoVZXC21)))))#1,(((((Rf1)#2,HKU312),#4),((((BtCoV2792005)#4,(#2,BtCoV2732005)))#6,((Rs3367)#3,(SARSCoVBJ1824,SARS)))),#3)))))
24 (#1,(#5,(#6,((((GuangxiPangolinP2V)#7,((GuangdongPangolin12019,#7),(((HuItalyTE48362020,HuWuhan2020),RaTG13),((BatCoVZC45,(BM4831BGR2008)#5),BatCoVZXC21)))))#1,((((Rf1)#2,(HKU312,#4)),((((BtCoV2792005)#4,(#2,BtCoV2732005)))#6,((Rs3367)#3,(SARSCoVBJ1824,SARS)))),#3)))))
24 (#1,(#5,(#6,((((GuangxiPangolinP2V)#7,((GuangdongPangolin12019,#7),(((HuItalyTE48362020,HuWuhan2020),RaTG13),((BatCoVZC45,(BM4831BGR2008)#5),BatCoVZXC21)))))#1,(((Rf1)#2,((HKU312,#4),((((BtCoV2792005)#4,(#2,BtCoV2732005)))#6,((Rs3367)#3,(SARSCoVBJ1824,SARS))))),#3)))))
24 (#1,(#5,(#6,((((GuangxiPangolinP2V)#7,((GuangdongPangolin12019,#7),(((HuItalyTE48362020,HuWuhan2020),RaTG13),((BatCoVZC45,(BM4831BGR2008)#5),BatCoVZXC21)))))#1,((Rf1)#2,(((HKU312,#4),((((BtCoV2792005)#4,(#2,BtCoV2732005)))#6,((Rs3367)#3,(SARSCoVBJ1824,SARS)))),#3))))))
24 (#1,(#5,(#6,((((GuangxiPangolinP2V)#7,((GuangdongPangolin12019,#7),(((HuItalyTE48362020,HuWuhan2020),RaTG13),((BatCoVZC45,(BM4831BGR2008)#5),((Rf1)#2,BatCoVZXC21))))))#1,(((HKU312,#4),((((BtCoV2792005)#4,(#2,BtCoV2732005)))#6,((Rs3367)#3,(SARSCoVBJ1824,SARS)))),#3)))))
24 (#1,(#5,(#6,((((GuangxiPangolinP2V)#7,((GuangdongPangolin12019,#7),(((HuItalyTE48362020,HuWuhan2020),RaTG13),((BatCoVZC45,(((Rf1)#2,BM4831BGR2008))#5),BatCoVZXC21)))))#1,(((HKU312,#4),((((BtCoV2792005)#4,(#2,BtCoV2732005)))#6,((Rs3367)#3,(SARSCoVBJ1824,SARS)))),#3)))))
24 (#1,(#5,(#6,((((GuangxiPangolinP2V)#7,((GuangdongPangolin12019,#7),(((HuItalyTE48362020,HuWuhan2020),RaTG13),((((Rf1)#2,BatCoVZC45),(BM4831BGR2008)#5),BatCoVZXC21)))))#1,(((HKU312,#4),((((BtCoV2792005)#4,(#2,BtCoV2732005)))#6,((Rs3367)#3,(SARSCoVBJ1824,SARS)))),#3)))))
24 (#1,(#5,(#6,((((GuangxiPangolinP2V)#7,((GuangdongPangolin12019,#7),(((HuItalyTE48362020,HuWuhan2020),RaTG13),(((Rf1)#2,(BatCoVZC45,(BM4831BGR2008)#5)),BatCoVZXC21)))))#1,(((HKU312,#4),((((BtCoV2792005)#4,(#2,BtCoV2732005)))#6,((Rs3367)#3,(SARSCoVBJ1824,SARS)))),#3)))))
24 (#1,(#5,(#6,((((GuangxiPangolinP2V)#7,((GuangdongPangolin12019,#7),(((HuItalyTE48362020,HuWuhan2020),RaTG13),((Rf1)#2,((BatCoVZC45,(BM4831BGR2008)#5),BatCoVZXC21))))))#1,(((HKU312,#4),((((BtCoV2792005)#4,(#2,BtCoV2732005)))#6,((Rs3367)#3,(SARSCoVBJ1824,SARS)))),#3)))))
24 (#1,(#5,(#6,((((GuangxiPangolinP2V)#7,((GuangdongPangolin12019,#7),(((HuItalyTE48362020,HuWuhan2020),((Rf1)#2,RaTG13)),((BatCoVZC45,(BM4831BGR2008)#5),BatCoVZXC21)))))#1,(((HKU312,#4),((((BtCoV2792005)#4,(#2,BtCoV2732005)))#6,((Rs3367)#3,(SARSCoVBJ1824,SARS)))),#3)))))
...
56 ((((((BatCoVZXC21,BatCoVZC45))#7,BM4831BGR2008),#5),(((((HKU312,(BtCoV2792005,#4)))#1,(Rs3367,((SARSCoVBJ1824)#4,SARS))))#6,((#3,(((((((BtCoV2732005)#3,((#6,HuItalyTE48362020),HuWuhan2020)),RaTG13),#7),GuangdongPangolin12019))#5,((GuangxiPangolinP2V)#2,Rf1))),#1))),#2)
56 ((((((BatCoVZXC21,BatCoVZC45))#7,BM4831BGR2008),#5),(((((HKU312,(BtCoV2792005,#4)))#1,(Rs3367,((SARSCoVBJ1824)#4,SARS))))#6,((#3,((((((((BtCoV2732005)#3,(#6,HuWuhan2020)),HuItalyTE48362020),RaTG13),#7),GuangdongPangolin12019))#5,((GuangxiPangolinP2V)#2,Rf1))),#1))),#2)
56 ((((((BatCoVZXC21,BatCoVZC45))#7,BM4831BGR2008),#5),(((((HKU312,(BtCoV2792005,#4)))#1,(Rs3367,((SARSCoVBJ1824)#4,SARS))))#6,((#3,(((((((BtCoV2732005)#3,(#6,(HuWuhan2020,HuItalyTE48362020))),RaTG13),#7),GuangdongPangolin12019))#5,((GuangxiPangolinP2V)#2,Rf1))),#1))),#2)
56 ((((((BatCoVZXC21,BatCoVZC45))#7,BM4831BGR2008),#5),(((((HKU312,(BtCoV2792005,#4)))#1,(Rs3367,((SARSCoVBJ1824)#4,SARS))))#6,((#3,(((((((BtCoV2732005)#3,((#6,HuWuhan2020),HuItalyTE48362020)),RaTG13),#7),GuangdongPangolin12019))#5,((GuangxiPangolinP2V)#2,Rf1))),#1))),#2)
56 ((((((BatCoVZXC21,BatCoVZC45))#7,BM4831BGR2008),#5),(((((HKU312,(BtCoV2792005,#4)))#1,(Rs3367,((SARSCoVBJ1824)#4,SARS))))#6,((#3,((((((((BtCoV2732005)#3,HuItalyTE48362020),(#6,HuWuhan2020)),RaTG13),#7),GuangdongPangolin12019))#5,((GuangxiPangolinP2V)#2,Rf1))),#1))),#2)
56 ((((((BatCoVZXC21,BatCoVZC45))#7,BM4831BGR2008),#5),((Rs3367,((((HKU312,(BtCoV2792005,#4)))#1,((SARSCoVBJ1824)#4,SARS)))#6),((#3,((((((((BtCoV2732005)#3,(#6,HuWuhan2020)),HuItalyTE48362020),RaTG13),#7),GuangdongPangolin12019))#5,((GuangxiPangolinP2V)#2,Rf1))),#1))),#2)
57 (#6,(((((#5,BtCoV2792005))#1,(((BatCoVZXC21)#5,((BM4831BGR2008)#6,BatCoVZC45)),(HuWuhan2020,HuItalyTE48362020))))#3,((((((((Rf1,BtCoV2732005))#7,(GuangxiPangolinP2V,(#4,GuangdongPangolin12019))))#2,(SARSCoVBJ1824,SARS)),(((RaTG13,#3))#4,Rs3367)),((HKU312,#1),#7)),#2)))
57 (#6,(((((#5,BtCoV2792005))#1,(((BatCoVZXC21,((BM4831BGR2008)#6,BatCoVZC45)))#5,(HuWuhan2020,HuItalyTE48362020))))#3,((((((((Rf1,BtCoV2732005))#7,(GuangxiPangolinP2V,(#4,GuangdongPangolin12019))))#2,(SARSCoVBJ1824,SARS)),(((RaTG13,#3))#4,Rs3367)),((HKU312,#1),#7)),#2)))
57 (#6,(((((#5,BtCoV2792005))#1,((BatCoVZXC21,(((BM4831BGR2008)#6,BatCoVZC45))#5),(HuWuhan2020,HuItalyTE48362020))))#3,((((((((Rf1,BtCoV2732005))#7,(GuangxiPangolinP2V,(#4,GuangdongPangolin12019))))#2,(SARSCoVBJ1824,SARS)),(((RaTG13,#3))#4,Rs3367)),((HKU312,#1),#7)),#2)))
57 (#6,(((((#5,BtCoV2792005))#1,(((BatCoVZXC21,((BM4831BGR2008)#6,BatCoVZC45)))#5,(HuWuhan2020,HuItalyTE48362020))))#3,(((((((((Rf1,BtCoV2732005))#7,(GuangxiPangolinP2V,(#4,GuangdongPangolin12019))))#2,(SARSCoVBJ1824,SARS)),((RaTG13,#3))#4),Rs3367),((HKU312,#1),#7)),#2)))
57 (#6,(((((#5,BtCoV2792005))#1,(((BatCoVZXC21,((BM4831BGR2008)#6,BatCoVZC45)))#5,(HuWuhan2020,HuItalyTE48362020))))#3,((((SARSCoVBJ1824,(((((Rf1,BtCoV2732005))#7,(GuangxiPangolinP2V,(#4,GuangdongPangolin12019))))#2,SARS)),(((RaTG13,#3))#4,Rs3367)),((HKU312,#1),#7)),#2)))
57 (#6,(((((#5,BtCoV2792005))#1,(((BatCoVZXC21,((BM4831BGR2008)#6,BatCoVZC45)))#5,(HuWuhan2020,HuItalyTE48362020))))#3,(((((SARSCoVBJ1824,((((Rf1,BtCoV2732005))#7,(GuangxiPangolinP2V,(#4,GuangdongPangolin12019))))#2),SARS),(((RaTG13,#3))#4,Rs3367)),((HKU312,#1),#7)),#2)))
62 ((((GuangdongPangolin12019)#5,(((((GuangxiPangolinP2V)#6,BM4831BGR2008))#7,(Rs3367,((SARSCoVBJ1824,SARS))#3)),((Rf1,BtCoV2732005),((#3,((BtCoV2792005,#6),(((((#5,RaTG13),(HuWuhan2020,HuItalyTE48362020)),#4))#1,HKU312))))#2))),#2),(#7,(#1,((BatCoVZXC21)#4,BatCoVZC45))))
96 ((((#1,((HKU312)#4,Rf1)))#2,(((((Rs3367)#1,BtCoV2732005),(((SARSCoVBJ1824,SARS))#5,(#6,(#3,GuangdongPangolin12019)))))#7,(BtCoV2792005,#4))),(((((HuWuhan2020,HuItalyTE48362020),RaTG13),(BatCoVZC45,BatCoVZXC21)))#3,((GuangxiPangolinP2V)#6,(#5,(#2,(BM4831BGR2008,#7))))))
119 ((GuangxiPangolinP2V,(((GuangdongPangolin12019)#3,(((((((#7,BM4831BGR2008))#4,HuItalyTE48362020))#5,(SARSCoVBJ1824,SARS)))#1,HuWuhan2020)))#6),((#2,(((BatCoVZC45,BatCoVZXC21),#5),(#4,(((((Rf1,BtCoV2732005))#7,RaTG13),(((HKU312,BtCoV2792005))#2,(#1,Rs3367))),#6)))),#3))
119 ((GuangxiPangolinP2V,(((GuangdongPangolin12019)#3,(((((((#7,BM4831BGR2008))#4,HuItalyTE48362020))#5,(SARSCoVBJ1824,SARS)))#1,HuWuhan2020)))#6),((#2,(((BatCoVZC45,BatCoVZXC21),#5),(#4,(((((Rf1,BtCoV2732005))#7,RaTG13),((HKU312,(BtCoV2792005)#2),(#1,Rs3367))),#6)))),#3))
122 ((#1,(#6,((#5,(((SARSCoVBJ1824,SARS))#3,(((((((BatCoVZC45)#4,BM4831BGR2008))#7,((GuangxiPangolinP2V)#6,((HuItalyTE48362020,(HuWuhan2020)#1),RaTG13))),Rf1),(GuangdongPangolin12019)#2),((((#4,BatCoVZXC21))#5,(HKU312,BtCoV2792005)),#2)))),((Rs3367,#3),BtCoV2732005)))),#7)
122 (#1,(#7,(#6,((#5,(((SARSCoVBJ1824,SARS))#3,(((((((BatCoVZC45)#4,BM4831BGR2008))#7,((GuangxiPangolinP2V)#6,(((HuItalyTE48362020,HuWuhan2020))#1,RaTG13))),Rf1),(GuangdongPangolin12019)#2),((((#4,BatCoVZXC21))#5,(HKU312,BtCoV2792005)),#2)))),((Rs3367,#3),BtCoV2732005)))))
122 (#1,(#7,(#6,((#5,(((SARSCoVBJ1824,SARS))#3,(((((((BatCoVZC45)#4,BM4831BGR2008))#7,((GuangxiPangolinP2V)#6,((HuItalyTE48362020,(HuWuhan2020)#1),RaTG13))),Rf1),(GuangdongPangolin12019)#2),((((#4,BatCoVZXC21))#5,(HKU312,BtCoV2792005)),#2)))),((Rs3367,#3),BtCoV2732005)))))
```


## Wheat dataset procesing

Using sampling and output directory wheat_guided. Single run with 10 initial networks.
```
supnet -G wheat_trees_clean -R6 -q10 -v4 --HC --hcstoptime=4 --hcrunstatsext --relaxed --hcsavewhenimproved --odtlabelled --guidetree "(1,2,3,4);(5,6,7,8);(9,10,11,12);(13,14,15);(16,17,18);(19,20,21);(22,23,24);(25,26);(27,28,29);(30,31,32,33);(34,35);(36,37,38,39);(40,41,42,43)" --outdirectory wheat_guided --displaytreesampling="0.001 0.01 0.03 0.05 0.1" --hcsamplerstats --hcsamplingmaxnetstonextlevel=1  --hcdetailedsummary  --autooutfiles
```   

In parallel with 10 cores and 20 runs:
```
parallel -q --jobs 10 --ungroup supnet -G wheat_trees_clean -R6 -q1 -v4 --HC --hcstoptime=4 --hcrunstatsext --relaxed --hcsavewhenimproved --odtlabelled --guidetree "(1,2,3,4);(5,6,7,8);(9,10,11,12);(13,14,15);(16,17,18);(19,20,21);(22,23,24);(25,26);(27,28,29);(30,31,32,33);(34,35);(36,37,38,39);(40,41,42,43)" --outdirectory wheat_guided --displaytreesampling="0.001 0.01 0.03 0.05 0.1" --hcsamplerstats --hcsamplingmaxnetstonextlevel=1  --hcdetailedsummary  --autooutfiles --outfile={1} ::: {1..20}
```

Results in wheat_guided directory.
```
cat wheat_guided/*.log | supnet -G wheat_trees_clean -N- --bestnetworks -v4
cat odt.log
```


