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

Display tree usage stats after HC run:
```  
> supnet --HC -g '(a,(b,c));(b,(a,c))' -R1 -q1; supnet -N odt.log --pdisplaytrees | sort | uniq -c | sort -k1 -n
   i: (c,((b)#1,(a,#1))) cost=2
   >: ((c,#1),((b)#1,a)) cost=1
   >: (#1,((c,(b)#1),a)) cost=0
Cost:0 Steps:17 Climbs:6 TopNetworks:1 Class:TreeChild TimeConsistency:0 HCruns:1 HCTime:8.9407e-05 Naive:18 Naivetime:9.05991e-06 DTcnt:36
      1 ((c,a),b)
      1 ((c,b),a)
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
((e,(d,c)),(b,a))
((b,a),(d,(e,c)))
```

### Guide clusters: used in quasi consensus 

Guide clusters are defined as a list of multifurcated trees separated by `;`. All clusters from the trees must be present in all quasi consensus tree and networks in hill climbing. Note that random networks are not generated using guide clusters (i.e., with `-r NUM`).

```
> supnet --guideclusters "(a,b,c);(d,e,(f,g))" -q4 --pnetworks
(((f,g),(e,d)),((a,b),c))
(((e,(g,f)),d),((a,c),b))
((d,(e,(f,g))),((a,b),c))
(((g,f),(d,e)),((a,c),b))
```

### Guide tree: forcing tree structure in networks

Similarly to guide clusters guide trees are defined as a multifurcated tree. 
A cluster present in every guide tree must be present as a cluster of a node in a network. Additionally, every such a node is a root of a subtree, i.e., no reticulation is allowed.
Note that random networks are not not generated using guide trees (i.e., with `-r NUM`).

```
> supnet --guidetree "((a,b,c),(d,e,f))" -q4 --pnetworks
(((f,e),d),((a,b),c))
((f,(d,e)),(c,(b,a)))
((c,(b,a)),(e,(f,d)))
((c,(a,b)),((d,e),f))
```

With reticulations; guide trees are separated by `;`.
```
> supnet -A10 -q2  --pnetworks --guidetree '(a,b);(d,e,f)' -R3
((((f,(e,d)),#3),(((c,h))#2,((j)#1,(i,(((a,b))#3,(#1,g)))))),#2)
((#1,(((d,f),e),j)),(#2,(((((h,c))#1,g))#2,((((a,b))#3,i),#3))))
```

### Random trees
    
Print two random networks with no reticulations (i.e., species trees) using species {a,b,c,d,e} (`-A5`). 
```
>  supnet -A5 -r2 --pnetworks
(a,((d,e),(c,b)))
(((a,c),(d,e)),b)
```

Print two quasi-consensus trees.
```
> supnet -g '(a,((b,c),d));(a,(b,d))' -q2 --pnetworks
(d,((c,a),b))
((d,c),(b,a))
```

### Random networks 

Print one random relaxed network with species {a,b,c}  and 5 reticulations.
```
>  supnet -A3 -R5 -r1 --pnetworks --relaxed
(((#4)#2)#5,((#3)#4,((#5,(a,((#2)#1,((c)#3,b)))),#1)))
```

Print one quasi-consensus network with two reticulations.
```
> supnet -G examples/gtrees.txt -q1 -R2 --pnetworks --pnetworkclusters
((c)#1,((d)#2,((#1,(a,#2)),b)))
a 
b 
c 
d 
a d 
a c d 
a b c d 
```

Print one quasi-consensus network with one reticulations and having guide tree clusters, and print all network clusters.
```
> supnet -q1 -R2 --pnetworks --guideclusters '((a,b,c),(d,e,f))' --pnetworkclusters
(((f,((d)#1,e)),#1),(b,(#2,((a)#2,c))))
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
a b c d e f 
```


### Fixed seed for random generator

Use `-z SEED` or `--randseed SEED` to set seed for random generator (srand). 

### Inserting reticulations

Insert 2 reticulations into a network; tree-child network in output (default).
```  
> supnet -R2 -n '(a,((d)#1,(b,(c,#1))))' --pnetworks
(#3,((a)#3,(((c,#1))#2,((d)#1,(b,#2)))))
```

Insert 8 reticulations into a network; general network `--general`
```  
> supnet -R8 -n '(a,((d)#1,(b,(c,#1))))' --pnetworks --general 
(((((#8)#3,(d)#1),(#3,#9)))#7,(((#2)#6,(((#6,#1))#9,(((((c)#2,#4))#8,((a)#5,(b)#4)),#5))),#7))
```

## Cost functions 

Available cost function are {DL,D,L,DC,RF,DCE}; the default is DC. 
Limitations: 
- In DP only DC is available, HC cannot be run using PD via BB.
- Experimental DTCache has only DC cost implementation.

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

 - `--HC`: run hill climbing heuristic for each initial network using cost function and print optimal cost, all optimal networks are written in odt.log file; summary stats are saved to odt.dat;
 - ` --hcnnimove`: use NNI instead of TailMoves
 - ` --hcpstats`: print stats after locating optimal networks after each HC runs
 - ` --hcpstatsext`: print extended stats after each HC run

 - `--hcrunstats`: print improvements stats after each HC runs
 - `--hcrunstatsext`: print improvements stats after locating optimal networks after each HC run
 - `--hcrunstatsalways`: print stats after each HC run
 - `--noodtfiles`: do not save odt.log and odt.dat with optimal networks
 - `--hcusenaive`: use only ODT naive algorithm in HCH (no BB)
 - `--hcstopinit=NUM`: stop when there is no new optimal network after HC with NUM initial networks
 - `--hcstopclimb=NUM`: stop in HC when there is no improvements after NUM steps
 - `--hcmaximprovements=NUM`: stop after NUM improvements in HC climb
 - `-O BASENAMEFILE, --outfiles`: base name of odt outputfiles .log and .dat
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
Cost:0 Steps:703 Climbs:8 TopNetworks:20 Class:TreeChild TimeConsistency:0 HCruns:1 HCTime:0.00435758 Naive:704 Naivetime:0.000680923 DTcnt:2795
```

See more examples in the last section.

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
Cost:2 Steps:96 Climbs:4 TopNetworks:3 Class:TreeChild TimeConsistency:0 HCruns:1 HCTime:0.000315189 Naive:97 Naivetime:6.34193e-05 DTcnt:194
```

Print improvements and equal cost networks `-v2`:
```  
> supnet -g "(a,(b,(c,d))); ((a,b),(c,d)); ((a,c),(b,d))" -r1 -R1 --HC -v2 --randseed 13
   i: ((#1,d),(((c,a))#1,b)) cost=5
   >: (((c,#1),d),((a)#1,b)) cost=2
   =: ((a)#1,(((c,#1),d),b)) cost=2
   =: ((a)#1,((c,d),(b,#1))) cost=2
Cost:2 Steps:96 Climbs:4 TopNetworks:3 Class:TreeChild TimeConsistency:0 HCruns:1 HCTime:0.000326633 Naive:97 Naivetime:6.29425e-05 DTcnt:194
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
Cost:0 Steps:18 Climbs:6 TopNetworks:1 Class:TreeChild TimeConsistency:0 HCruns:1 HCTime:0.000101805 Naive:19 Naivetime:9.77516e-06 DTcnt:38
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
Cost:9 Steps:2151 Climbs:10 TopNetworks:103 Class:TreeChild TimeConsistency:0 HCruns:2 HCTime:0.105543 Naive:2153 Naivetime:0.00626135 DTcnt:17224
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
Cost:0 Steps:2944 Climbs:16 TopNetworks:53 Class:TreeChild TimeConsistency:0 HCruns:2 HCTime:0.0629137 Naive:2946 Naivetime:0.00739336 DTcnt:23410
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
14198	((((b,a))#1,c),#1)
28763	((b)#1,(a,(#1,c)))
14163	(((a)#1,c),(b,#1))
28406	((((a)#1,c),#1),b)
14470	(#1,((c)#1,(a,b)))
unique=5 all=0
```

## Other options, print, debug, test, etc.

 - `--pspeciesdictionary`: list species dictionary
 - `--pdetailed`: print debug tree and network structures
 - `--reachablenodescnt`: for each v in V(N), print the number of nodes reachable from v
 - `--reachableleafvescnt`: for each v in V(N), print the number of leaves reachable from v
 - `--pstsubtrees`: print species tree subtrees
 - `--maxdisplaytreecachesize`: set limit for the size of the cache of display tree nodes (only if DTCACHE is enabled)

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

## Coronavirus dataset processing

Infering tree-child networks without time-consistency, start from 10 networks `-q10`.
Use: `--hcsavewhenimproved` to save odt files after each improvement
and  `--hcstoptime 30` to stop a single HC if there is no improvement after 30 seconds.

```
make supnet_dtcache
supnet_dtcache -G corona.txt -R7 -q10 --HC --hcrunstats --outfiles corona --hcstoptime 30 --randseed 13 --hcsavewhenimproved
```
Results in `corona.log` and `corona.dat`.


Inferring relaxed time-consistent networks, start from 10 networks `-q10`.
```
supnet_dtcache -G corona.txt -R7 -q10 --HC --hcrunstats --timeconsistent --relaxed
```

Inferring relaxed time-consistent networks, start from 10 networks `-q10`, save networks to first.log.
```
supnet_dtcache -G corona.txt -R7 -q10 --HC -v1 --timeconsistent --relaxed --outfiles first --hcsavewhenimproved
supnet_dtcache -G corona.txt -N first.log --hcsavewhenimproved  
supnet_dtcache -G corona.txt -N first.log --HC -v1 --timeconsistent --relaxed --hcsavewhenimproved
```


### Parallel processing 

TODO: multithreded supnet

Use parallel with 10 jobs to run 20 supnet computations.
```
mkdir -p corona
parallel --jobs 10 --ungroup supnet_dtcache -G corona.txt -R7 -q1 --HC --hcrunstats --outfiles corona/{1} --hcstoptime 1  --hcsavewhenimproved --odtlabelled ::: {1..20}
```

Download `csvmanip` to merge dat files:
```
git clone git@github.com:ppgorecki/csvmanip.git 
```

Use `csvmanip` to merge dat files into csv:
```
> csvmanip/csvmanip.py corona/*.dat
Id,Source,optcost,time,hctime,mergetime,topnets,class,timeconsistency,improvements,steps,bbruns,startnets,memoryMB,dtcnt,naivetime,naivecnt,bbnaivecnt,bbnaivetime,bbdpcnt,bbdptime,randseed
1,corona/1.dat,139,1.00023,1.00021,1.90735e-05,5,0,0,48,3298,0,1,72,422272,0.928381,3299,0,0,0,0,1030374630
2,corona/2.dat,113,1.04761,1.04758,2.24113e-05,1,0,0,52,3122,0,1,129,399744,0.993293,3123,0,0,0,0,1030374632
3,corona/3.dat,151,1.00013,1.00012,1.57356e-05,1,0,0,38,2751,0,1,100,352256,0.980743,2752,0,0,0,0,1030374633
4,corona/4.dat,81,1.00017,1.00016,1.19209e-05,1,0,0,74,4714,0,1,77,603520,0.95827,4715,0,0,0,0,1030374634
5,corona/5.dat,102,1.00017,1.00016,1.21593e-05,3,0,0,68,4068,0,1,60,520832,0.959025,4069,0,0,0,0,1030374634
6,corona/6.dat,97,1.00028,1.00027,1.19209e-05,1,0,0,46,4107,0,1,94,525824,0.972525,4108,0,0,0,0,1030374636
7,corona/7.dat,121,1.00034,1.00033,6.91414e-06,17,0,0,42,3289,0,1,95,421120,0.878576,3290,0,0,0,0,1030374637
8,corona/8.dat,129,1.00025,1.00024,1.62125e-05,2,0,0,46,3143,0,1,86,402432,0.964586,3144,0,0,0,0,1030374639
9,corona/9.dat,160,1.00054,1.00052,1.71661e-05,1,0,0,34,2640,0,1,105,338048,0.977057,2641,0,0,0,0,1030374641
10,corona/10.dat,116,1.00017,1.00016,1.23978e-05,1,0,0,58,2871,0,1,101,367616,0.970199,2872,0,0,0,0,1030374642
11,corona/11.dat,156,1.00029,1.00028,1.14441e-05,1,0,0,40,2692,0,1,91,344704,0.978215,2693,0,0,0,0,1030375740
12,corona/12.dat,125,1.00007,1.00007,6.19888e-06,25,0,0,42,3236,0,1,64,414336,0.884311,3237,0,0,0,0,1030375742
13,corona/13.dat,142,1.00032,1.0003,1.85966e-05,1,0,0,32,2688,0,1,94,344192,0.973474,2689,0,0,0,0,1030375752
14,corona/14.dat,205,1.00016,1.00015,1.50204e-05,1,0,0,54,2661,0,1,101,340736,0.975555,2662,0,0,0,0,1030375762
15,corona/15.dat,105,1.00003,1.00002,1.38283e-05,4,0,0,64,3483,0,1,110,445952,0.916026,3484,0,0,0,0,1030375768
16,corona/16.dat,98,1.02999,1.02998,1.35899e-05,3,0,0,38,3626,0,1,128,464256,1.0059,3627,0,0,0,0,1030375776
17,corona/17.dat,95,1.00003,1.00001,1.69277e-05,32,0,0,58,2649,0,1,58,339200,0.750013,2650,0,0,0,0,1030375795
18,corona/18.dat,94,1.00014,1.00012,2.02656e-05,2,0,0,44,3113,0,1,123,398592,0.973034,3114,0,0,0,0,1030375804
19,corona/19.dat,132,1.00026,1.00024,2.09808e-05,1,0,0,36,2902,0,1,102,371584,0.973881,2903,0,0,0,0,1030375813
20,corona/20.dat,134,1.00008,1.00006,1.50204e-05,8,0,0,52,3830,0,1,86,490368,0.951913,3831,0,0,0,0,1030375844
```

Print cost + networks with sort:
```
> cat corona/*.log | supnet -N- -G corona.txt --odtcost | sort -k1 -n
81 ((BM4831BGR2008,((((((HKU312,(BtCoV2792005)#2),((SARS,(#6,Rf1)),#5)),Rs3367))#1,(((((BtCoV2732005)#6,(((BatCoVZXC21,BatCoVZC45))#3,(HuWuhan2020,#7))))#4,GuangxiPangolinP2V),#1)),#2)),((#3,(((((HuItalyTE48362020)#7,SARSCoVBJ1824))#5,RaTG13),#4)),GuangdongPangolin12019))
94 ((#5,(Rs3367,((((#6,(((((HuWuhan2020,HuItalyTE48362020))#2,BatCoVZC45))#5,(((((BtCoV2732005,Rf1))#7,GuangdongPangolin12019))#6,((BatCoVZXC21)#3,(RaTG13,#4))))))#1,(#7,((GuangxiPangolinP2V)#4,(SARSCoVBJ1824,SARS)))),(BtCoV2792005,HKU312)))),(((BM4831BGR2008,#3),#1),#2))
94 ((#5,(Rs3367,((((#6,((((HuWuhan2020,(HuItalyTE48362020)#2),BatCoVZC45))#5,(((((BtCoV2732005,Rf1))#7,GuangdongPangolin12019))#6,((BatCoVZXC21)#3,(RaTG13,#4))))))#1,(#7,((GuangxiPangolinP2V)#4,(SARSCoVBJ1824,SARS)))),(BtCoV2792005,HKU312)))),(((BM4831BGR2008,#3),#1),#2))
...
```


### Check network classes in the result:
```
supnet -N odt.log --detectclass
```



## Development variants

### Display tree cache 

When compiling with macro DTCACHE `make supnet_dtcache` option `--maxdisplaytreecachesize SIZE` will limit the number of stored nodes in the cache. The default is `1000000` nodes.

```
make supnet_dtcache
supnet_dtcache  --maxdisplaytreecachesize 5000000 ...
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
Cost:3 Steps:2066 Climbs:43 TopNetworks:1 Class:TreeChild TimeConsistency:0 HCruns:8 HCTime:0.00680876 Naive:2074 Naivetime:0.00165772 DTcnt:4148
```

Print stats after completed HC run for each initial network when only when new optimal network is found.
```
> supnet -g "(a,(b,(c,d))); ((a,b),(c,d)); ((a,c),(b,d)); ((a,d),e)" -r8 -R1 --HC --hcrunstats -v0 --randseed 13
1. Cost:4 Steps:266 Climbs:3 TopNetworks:5 Class:TreeChild TimeConsistency:0
 New optimal cost: 4
5. Cost:3 Steps:182 Climbs:6 TopNetworks:1 Class:TreeChild TimeConsistency:0
 New optimal cost: 3
Cost:3 Steps:2066 Climbs:43 TopNetworks:1 Class:TreeChild TimeConsistency:0 HCruns:8 HCTime:0.00700283 Naive:2074 Naivetime:0.00165749 DTcnt:4148
```

Print stats after completed HC run for each initial network when a new optimal network including equal networks is found.
```
> supnet -g "(a,(b,(c,d))); ((a,b),(c,d)); ((a,c),(b,d)); ((a,d),e)" -r8 -R1 --HC --hcrunstatsext -v0 --randseed 13
1. Cost:4 Steps:266 Climbs:3 TopNetworks:5 Class:TreeChild TimeConsistency:0
 New optimal cost: 4
5. Cost:3 Steps:182 Climbs:6 TopNetworks:1 Class:TreeChild TimeConsistency:0
 New optimal cost: 3
Cost:3 Steps:2066 Climbs:43 TopNetworks:1 Class:TreeChild TimeConsistency:0 HCruns:8 HCTime:0.00709414 Naive:2074 Naivetime:0.00165248 DTcnt:4148
```

Print stats after completed HC run.
```
> supnet -g "(a,(b,(c,d))); ((a,b),(c,d)); ((a,c),(b,d)); ((a,d),e)" -r8 -R1 --HC --hcrunstatsalways -v0 --randseed 13
1. Cost:4 Steps:266 Climbs:3 TopNetworks:5 Class:TreeChild TimeConsistency:0
 New optimal cost: 4
2. Cost:6 Steps:160 Climbs:4 TopNetworks:2 Class:TreeChild TimeConsistency:0
 New optimal networks:0 Total:5
3. Cost:4 Steps:247 Climbs:1 TopNetworks:5 Class:TreeChild TimeConsistency:0
 New optimal networks:0 Total:5
4. Cost:4 Steps:255 Climbs:3 TopNetworks:5 Class:TreeChild TimeConsistency:0
 New optimal networks:0 Total:5
5. Cost:3 Steps:182 Climbs:6 TopNetworks:1 Class:TreeChild TimeConsistency:0
 New optimal cost: 3
6. Cost:3 Steps:371 Climbs:9 TopNetworks:1 Class:TreeChild TimeConsistency:0
 New optimal networks:0 Total:1
7. Cost:4 Steps:292 Climbs:7 TopNetworks:5 Class:TreeChild TimeConsistency:0
 New optimal networks:0 Total:1
8. Cost:4 Steps:293 Climbs:6 TopNetworks:3 Class:TreeChild TimeConsistency:0
 New optimal networks:0 Total:1
Cost:3 Steps:2066 Climbs:43 TopNetworks:1 Class:TreeChild TimeConsistency:0 HCruns:8 HCTime:0.00736785 Naive:2074 Naivetime:0.00177097 DTcnt:4148
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
Cost:0 Steps:410 Climbs:2 TopNetworks:19 Class:TreeChild TimeConsistency:0 HCruns:1 HCTime:0.0021975 Naive:411 Naivetime:0.000174046 DTcnt:809
Optimal networks saved: odt.log
Stats data saved: odt.dat
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
Cost:0 Steps:410 Climbs:2 TopNetworks:19 Class:TreeChild TimeConsistency:0 HCruns:1 HCTime:0.00220919 Naive:411 Naivetime:0.000168085 DTcnt:809
Optimal networks saved: odt.log
Stats data saved: odt.dat
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
Cost:2 Steps:96 Climbs:4 TopNetworks:3 Class:TreeChild TimeConsistency:0 HCruns:1 HCTime:0.000408649 Naive:97 Naivetime:6.4373e-05 DTcnt:194
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
Cost:4 Steps:7 Climbs:4 TopNetworks:3 Class:TreeChild TimeConsistency:0 HCruns:1 HCTime:7.31945e-05 Naive:8 Naivetime:7.86781e-06 DTcnt:16
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
Cost:0 Steps:121 Climbs:8 TopNetworks:5 Class:TreeChild TimeConsistency:0 HCruns:1 HCTime:0.000751972 Naive:122 Naivetime:0.000237942 DTcnt:966
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
Cost:0 Steps:233 Climbs:6 TopNetworks:2 Class:Relaxed TimeConsistency:0 HCruns:1 HCTime:0.00120163 Naive:234 Naivetime:0.000477314 DTcnt:1866
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
Cost:0 Steps:52673 Climbs:4 TopNetworks:680 Class:General TimeConsistency:0 HCruns:1 HCTime:3.88446 Naive:52674 Naivetime:0.120052 DTcnt:420111
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
Cost:0 Steps:176 Climbs:10 TopNetworks:1 Class:General TimeConsistency:1 HCruns:1 HCTime:0.00196719 Naive:177 Naivetime:0.00036931 DTcnt:1414
```

### Global cache of networks in HC runs

If `--globaldagcache` is set then networks already visited will be omitted in HC.
This feature is not optimized. For large sets of networks, it may be slower that non-cache option.

### odt.dat in labelled format `--odtlabelled`

```
> supnet -g "(a,(b,(c,d))); ((a,b),(c,d))" -r1 -R2 --randseed 13 --HC --odtlabelled -v0
Cost:0 Steps:703 Climbs:8 TopNetworks:20 Class:TreeChild TimeConsistency:0 HCruns:1 HCTime:0.00411057 Naive:704 Naivetime:0.000678778 DTcnt:2795
```

```
> cat odt.dat
optcost=0
time=0.00411081
hctime=0.00411057
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
naivetime=0.000678778
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
1. Cost:16 Steps:46020 Climbs:26 TopNetworks:60 Class:TreeChild TimeConsistency:0
 New optimal cost: 16
Cost:16 Steps:46020 Climbs:52 TopNetworks:60 Class:TreeChild TimeConsistency:0 HCruns:1 HCTime:20.5688 Naive:46021 Naivetime:14.9772 DTcnt:11781376
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
Cost:0 Steps:470980 Climbs:2623 TopNetworks:254 Class:TreeChild TimeConsistency:0 HCruns:1000 HCTime:8.92229 Naive:471980 Naivetime:0.851009 DTcnt:3687318
```

### Using stopping criterion (--hcstopinit=1000) with quasi-consensus rand networks.

```  
> supnet -g "(a,(b,(c,d))); ((a,b),(c,d))" -q-1 -R3 --HC --hcrunstats --hcstopinit=1000 --randseed 13 | tail -10
   >: ((d)#2,(#3,((((a)#1,(c,#2)))#3,(#1,b)))) cost=1
   >: ((d)#2,(#3,(#1,((((a)#1,(c,#2)))#3,b)))) cost=0
   i: (#3,(((#1,((c)#2,((d)#1,a))))#3,(#2,b))) cost=4
   >: (#3,((((c)#2,((d)#1,a)))#3,(#2,(#1,b)))) cost=3
   >: (#3,(((((#1,c))#2,((d)#1,a)))#3,(#2,b))) cost=1
   >: (#3,(#2,(((((#1,c))#2,((d)#1,a)))#3,b))) cost=0
   i: (((((d)#2,((b)#3,c)),#2))#1,((#1,a),#3)) cost=1
   >: ((((d)#2,((b)#3,(c,#2))))#1,((#1,a),#3)) cost=0
   i: (#2,((((((((d)#1,a))#3,c),#1))#2,b),#3)) cost=0
Cost:0 Steps:492723 Climbs:2734 TopNetworks:254 Class:TreeChild TimeConsistency:0 HCruns:1047 HCTime:9.3128 Naive:493770 Naivetime:0.88283 DTcnt:3857600
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
HC run stopped due to timeout, cost=25
1. Cost:25 Steps:1556 Climbs:18 TopNetworks:10 Class:TreeChild TimeConsistency:0
 New optimal cost: 25
Cost:25 Steps:1556 Climbs:36 TopNetworks:10 Class:TreeChild TimeConsistency:0 HCruns:1 HCTime:0.500043 Naive:1557 Naivetime:0.481957 DTcnt:398592
Optimal networks saved: odt.log
Stats data saved: odt.dat
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
HC run stopped due to timeout, cost=25
1. Cost:25 Steps:1488 Climbs:18 TopNetworks:4 Class:TreeChild TimeConsistency:0
 New optimal cost: 25
Cost:25 Steps:1488 Climbs:36 TopNetworks:4 Class:TreeChild TimeConsistency:0 HCruns:1 HCTime:0.500284 Naive:1489 Naivetime:0.468969 DTcnt:381184
Optimal networks saved: odt.log
Stats data saved: odt.dat
```

## More examples

Print min total cost 10 random gene trees vs random tree-child network with 5 reticulations over 8 species; print the initial network.

```
> supnet -r10 -A8 --pnetworks  | supnet -G- -r1 -A8 -R5 --pnetworks --odtnaivecost
(#2,(((#4,((((h)#1,((a)#4,d)))#3,b)))#5,((e,c),(#1,(((f,#5))#2,(#3,g))))))
61
```


