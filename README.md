# SupNet

Phylogenetic networks: visualization, evaluation, and inference

## Compiling

Default version up to 128 species and 64 reticulations
```
make supnet
```

Up to 256 species and 64 reticulations
```
make supnet256
```

Up to 64 species and 64 reticulations
```
make supnet64
```

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
3 ((c,b),a)
```

Display tree usage stats after HC run:
```  
> supnet --HC -g '(a,(b,c));(b,(a,c))' -R1 -q1; supnet -N odt.log --pdisplaytrees | sort | uniq -c | sort -k1 -n
   i: (#1,((b)#1,(c,a))) cost=1
   >: (#1,((c,(b)#1),a)) cost=0
Cost:0 TopNets:1 Steps:17 Climbs:4 CRuns:1 CTime:6.93798e-05 Naive:18 Naivetime:1.0252e-05 DTcnt:36 Class:TreeChild TimeConsistency:0
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
 - `-A SPECIESNUM`: define `SPECIESNUM` species a,b,...z, a1,b1, ..., z1, a2, b2,...
 - `--uniformedgesampling`: draw uniformly a pair of edges to create random reticulation; the default: first draw the source edge, then the uniformly the - destination edge given the source
 - `--preserveroot`: preserve root in quasi-consensus and in HC algorithm
 - `--randseed, -z SEED`: set seed for random generator (srand)

If NUM is -1, the network generator will provide network when requested (e.g., when used with `-K/--hcstop`, see `--HC` option).

By default, the generator produces tree-child networks without time-consistency restriction. 

Print 2 quasi consensus trees with preserved split of the root.
```
> supnet -g '(a,((b,c),d));(a,(b,e))' -s'((a,b),(c,(d,e)))' -q2 --preserveroot --pnetworks
(((e,d),c),(a,b))
((b,a),(c,(e,d)))
```

### Guide clusters: used in quasi consensus 

Guide clusters are defined as a list of multifurcated trees separated by `;`. All clusters from the trees must be present in all quasi consensus tree and networks in hill climbing. 

```
> supnet --guideclusters "(a,b,c);(d,e,(f,g))" -q4 --pnetworks
(((d,e),(f,g)),(b,(c,a)))
((d,(e,(g,f))),(c,(a,b)))
((d,(e,(g,f))),((b,a),c))
(((c,b),a),((d,e),(f,g)))
```

### Guide tree: forcing tree structure in networks

Similarly to guide clusters guide trees are defined as a multifurcated tree. 
A cluster present in every guide tree must be present as a cluster of a node in a network. Additionally, every such a node is a root of a subtree, i.e., no reticulation is allowed.
Note that random networks are not not generated using guide trees (i.e., with `-r NUM`).

```
> supnet --guidetree "((a,b,c),(d,e,f))" -q4 --pnetworks
((b,(c,a)),(f,(e,d)))
((c,(b,a)),((f,e),d))
(((f,e),d),((c,a),b))
((b,(c,a)),(e,(d,f)))
```

With reticulations; guide trees are separated by `;`.
```
> supnet -A10 -q2  --pnetworks --guidetree '(a,b);(d,e,f)' -R3
(#1,((((i)#2,(h,((#2,(((e,(f,d)))#3,((b,a),#3))),g))))#1,(c,j)))
(#1,(((g,((h)#2,(((#2,j))#1,((f,d),e)))))#3,((#3,c),((a,b),i))))
```

### Random trees
    
Print two random networks with no reticulations (i.e., species trees) using species {a,b,c,d,e} (`-A5`). 
```
>  supnet -A5 -r2 --pnetworks
(c,(e,((a,d),b)))
((a,b),((c,d),e))
```

Print two quasi-consensus trees.
```
> supnet -g '(a,((b,c),d));(a,(b,d))' -q2 --pnetworks
((d,(c,b)),a)
(a,(d,(c,b)))
```

### Random networks 

Print one random relaxed network with species {a,b,c}  and 5 reticulations.
```
>  supnet -A3 -R5 -r1 --pnetworks --relaxed
(#5,(#2,((((((#4)#5,((a)#4,b)))#1)#2)#3,((#3,c),#1))))
```

Print one quasi-consensus network with two reticulations.
```
> supnet -G examples/gtrees.txt -q1 -R2 --pnetworks --pnetworkclusters
((((d)#1,b))#2,((#2,c),(#1,a)))
1 a 
1 b 
1 c 
2 d 
1 a d 
2 b d 
1 b c d 
2 a b c d 
```

Print one quasi-consensus network with one reticulations and having guide tree clusters, and print all network clusters.
```
> supnet -q1 -R2 --pnetworks --guideclusters '((a,b,c),(d,e,f))' --pnetworkclusters
((e)#1,((d)#2,((#2,(#1,f)),((c,a),b))))
2 a 
2 b 
2 c 
4 d 
5 e 
2 f 
1 a c 
1 e f 
1 a b c 
1 d e f 
3 a b c d e f 
```


### Fixed seed for random generator

Use `-z SEED` or `--randseed SEED` to set seed for random generator (srand). 

### Inserting reticulations

Insert 2 reticulations into a network; tree-child network in output (default).
```  
> supnet -R2 -n '(a,((d)#1,(b,(c,#1))))' --pnetworks
((a)#2,(#2,((b)#3,((d)#1,(#3,(c,#1))))))
```

Insert 8 reticulations into a network; general network `--general`
```  
> supnet -R8 -n '(a,((d)#1,(b,(c,#1))))' --pnetworks --general 
(#6,(((((b)#7,(((((((#7)#9,c),#8))#4,#2))#5,a)))#6,((#9)#2,(d)#1)),(((#1)#8,((#4)#3,#5)),#3)))
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
Cost:0 TopNets:20 Steps:703 Climbs:8 CRuns:1 CTime:0.00284195 Naive:704 Naivetime:0.000617743 DTcnt:2795 Class:TreeChild TimeConsistency:0
```

See more examples in the last section.

## Saving extended newick and inserting prefix in reticulation ids

```
> supnet -g "(a,(b,(c,d))); ((a,b),(c,d))" -r1 -R2 --HC --randseed 13 --saveextnewick --retidprefix=H
   i: ((#H1,d),(((((a)#H2,c),#H2))#H1,b)) cost=4
   >: (((#H1,d),#H2),((((a)#H2,c))#H1,b)) cost=2
   >: (#H2,((#H1,d),((((a)#H2,c))#H1,b))) cost=1
   >: (#H2,((#H1,d),((a)#H2,((c)#H1,b)))) cost=0
Cost:0 TopNets:20 Steps:703 Climbs:8 CRuns:1 CTime:0.00298977 Naive:704 Naivetime:0.000622511 DTcnt:2795 Class:TreeChild TimeConsistency:0
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
Cost:2 TopNets:3 Steps:96 Climbs:4 CRuns:1 CTime:0.00013423 Naive:97 Naivetime:5.48363e-05 DTcnt:194 Class:TreeChild TimeConsistency:0
```

Print improvements and equal cost networks `-v2`:
```  
> supnet -g "(a,(b,(c,d))); ((a,b),(c,d)); ((a,c),(b,d))" -r1 -R1 --HC -v2 --randseed 13
   i: ((#1,d),(((c,a))#1,b)) cost=5
   >: (((c,#1),d),((a)#1,b)) cost=2
   =: ((a)#1,(((c,#1),d),b)) cost=2
   =: ((a)#1,((c,d),(b,#1))) cost=2
Cost:2 TopNets:3 Steps:96 Climbs:4 CRuns:1 CTime:0.000133038 Naive:97 Naivetime:5.36442e-05 DTcnt:194 Class:TreeChild TimeConsistency:0
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
Cost:0 TopNets:1 Steps:18 Climbs:6 CRuns:1 CTime:6.17504e-05 Naive:19 Naivetime:1.00136e-05 DTcnt:38 Class:TreeChild TimeConsistency:0
```

To show details on DP & BB algorithms use `-v` with 4, 5 and 6. See examples below.


HC with guide trees.
```

> supnet -g "(a,(b,c));((a,(b,g)),(c,(d,(e,(f,h)))));" -R1 --HC -v1 -q2 --guidetree '(a,b);(d,e,(f,g))' -R3 --randseed 13
   i: ((c)#2,(((h)#3,((((a,b),#2))#1,((e,(d,(f,g))),#1))),#3)) cost=12
   >: ((c)#2,(((h)#3,((((a,b),#2))#1,(((d,e),(f,g)),#1))),#3)) cost=11
   >: ((c)#2,((((h)#3,((((a,b),#2))#1,((d,e),(f,g)))),#1),#3)) cost=9
   i: (#3,(#2,((((((d,(f,g)),e))#2,(#1,c)))#3,(((b,a))#1,h)))) cost=11
   >: (#3,(#2,((((((f,g),(d,e)))#2,(#1,c)))#3,(((b,a))#1,h)))) cost=10
   >: (#3,((((((f,g),(d,e)))#2,(#1,c)))#3,(((b,a))#1,(#2,h)))) cost=9
Cost:9 TopNets:105 Steps:2238 Climbs:9 CRuns:2 CTime:0.100263 Naive:2240 Naivetime:0.00695682 DTcnt:17920 Class:TreeChild TimeConsistency:0
```

HC with guide clusters.
```
> supnet -g "(a,(b,c));((a,(b,g)),(c,(d,(e,f))));" -R1 --HC -v1 -q2 --guideclusters '(a,b,g);(d,e,f)' -R3 --randseed 13
   i: ((a)#2,((((d,(#3,((e)#3,f))))#1,(#1,c)),(#2,(g,b)))) cost=2
   >: (((a)#2,((((d,(#3,((e)#3,f))))#1,c),(#2,(g,b)))),#1) cost=1
   i: (((((g,#2))#1,c),((b)#2,(a,#1))),(((d)#3,(e,f)),#3)) cost=2
   >: ((((g,#2))#1,(((b)#2,c),(a,#1))),(((d)#3,(e,f)),#3)) cost=1
   >: (((g,#2))#1,(((d,(e,f)))#3,(((b)#2,(#3,c)),(a,#1)))) cost=0
Cost:0 TopNets:6 Steps:2078 Climbs:8 CRuns:2 CTime:0.0344622 Naive:2080 Naivetime:0.00645876 DTcnt:16640 Class:TreeChild TimeConsistency:0
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
28642	((((a)#1,b),c),#1)
14420	(((a)#1,(c,b)),#1)
14300	(((a)#1,b),(c,#1))
28289	(((b)#1,(#1,c)),a)
14349	((((c,b))#1,a),#1)
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
Cost:3 TopNets:1 Steps:2066 Climbs:43 CRuns:8 CTime:0.00282931 Naive:2074 Naivetime:0.00113845 DTcnt:4148 Class:TreeChild TimeConsistency:0
```

Print stats after completed HC run for each initial network when only when new optimal network is found.
```
> supnet -g "(a,(b,(c,d))); ((a,b),(c,d)); ((a,c),(b,d)); ((a,d),e)" -r8 -R1 --HC --hcrunstats -v0 --randseed 13
1.  Steps:0 Climbs:0 New optimal cost: 4
5.  Steps:0 Climbs:0 New optimal cost: 3
Cost:3 TopNets:1 Steps:2066 Climbs:43 CRuns:8 CTime:0.00295901 Naive:2074 Naivetime:0.00128436 DTcnt:4148 Class:TreeChild TimeConsistency:0
```

Print stats after completed HC run for each initial network when a new optimal network including equal networks is found.
```
> supnet -g "(a,(b,(c,d))); ((a,b),(c,d)); ((a,c),(b,d)); ((a,d),e)" -r8 -R1 --HC --hcrunstatsext -v0 --randseed 13
1.  Steps:0 Climbs:0 New optimal cost: 4
5.  Steps:0 Climbs:0 New optimal cost: 3
Cost:3 TopNets:1 Steps:2066 Climbs:43 CRuns:8 CTime:0.00298548 Naive:2074 Naivetime:0.00119138 DTcnt:4148 Class:TreeChild TimeConsistency:0
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
Cost:3 TopNets:1 Steps:2066 Climbs:43 CRuns:8 CTime:0.00304008 Naive:2074 Naivetime:0.00122499 DTcnt:4148 Class:TreeChild TimeConsistency:0
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
Cost:0 TopNets:19 Steps:410 Climbs:2 CRuns:1 CTime:0.00153136 Naive:411 Naivetime:0.000172853 DTcnt:809 Class:TreeChild TimeConsistency:0
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
Cost:0 TopNets:19 Steps:410 Climbs:2 CRuns:1 CTime:0.00155401 Naive:411 Naivetime:0.000174999 DTcnt:809 Class:TreeChild TimeConsistency:0
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
Cost:2 TopNets:3 Steps:96 Climbs:4 CRuns:1 CTime:0.000230551 Naive:97 Naivetime:5.65052e-05 DTcnt:194 Class:TreeChild TimeConsistency:0
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
Cost:4 TopNets:3 Steps:7 Climbs:4 CRuns:1 CTime:5.22137e-05 Naive:8 Naivetime:8.58307e-06 DTcnt:16 Class:TreeChild TimeConsistency:0
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
Cost:0 TopNets:5 Steps:121 Climbs:8 CRuns:1 CTime:0.000532389 Naive:122 Naivetime:0.000237465 DTcnt:966 Class:TreeChild TimeConsistency:0
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
Cost:0 TopNets:2 Steps:221 Climbs:6 CRuns:1 CTime:0.000851154 Naive:222 Naivetime:0.000449181 DTcnt:1770 Class:Relaxed TimeConsistency:0
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
Cost:0 TopNets:680 Steps:52673 Climbs:4 CRuns:1 CTime:3.83927 Naive:52674 Naivetime:0.10928 DTcnt:420111 Class:General TimeConsistency:0
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
Cost:0 TopNets:1 Steps:176 Climbs:10 CRuns:1 CTime:0.00164437 Naive:177 Naivetime:0.000346184 DTcnt:1414 Class:General TimeConsistency:1
```

### Global cache of networks in HC runs

If `--globaldagcache` is set then networks already visited will be omitted in HC.
This feature is not optimized. For large sets of networks, it may be slower that non-cache option.

### odt.dat in labelled format `--odtlabelled`

```
> supnet -g "(a,(b,(c,d))); ((a,b),(c,d))" -r1 -R2 --randseed 13 --HC --odtlabelled -v0
Cost:0 TopNets:20 Steps:703 Climbs:8 CRuns:1 CTime:0.00279284 Naive:704 Naivetime:0.000623226 DTcnt:2795 Class:TreeChild TimeConsistency:0
```

```
> cat odt.dat
optcost=0
time=0.00279307
climbtime=0.00279284
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
naivetime=0.000623226
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
   i: ((((c)#5,e))#2,(#7,(#3,(((((((g)#4,(#4,((((h)#3,((((f,#5))#6,d),b)))#8,i))))#1,(#6,(#2,(j,#1)))),#8))#7,a)))) cost=46
   >: ((((c)#5,e))#2,(#7,((#3,(((((((g)#4,(#4,((((h)#3,(((f)#6,d),b)))#8,i))))#1,(#6,(#2,(j,#1)))),#8))#7,a)),#5))) cost=45
   >: ((((c)#5,e))#2,(#7,((#4,(#3,(((((((g)#4,((((h)#3,(((f)#6,d),b)))#8,i)))#1,(#6,(#2,(j,#1)))),#8))#7,a))),#5))) cost=42
   >: ((((c)#5,e))#2,(#7,((#6,(#4,(#3,(((((((g)#4,((((h)#3,(((f)#6,d),b)))#8,i)))#1,(#2,(j,#1))),#8))#7,a)))),#5))) cost=37
   >: ((((c)#5,e))#2,(#7,(#6,(#4,(#3,(((((((g)#4,((((h)#3,((((f,#5))#6,d),b)))#8,i)))#1,(#2,(j,#1))),#8))#7,a)))))) cost=36
   >: ((((c)#5,e))#2,(#7,(#6,(#4,(#3,(((((((g)#4,((((h)#3,b))#8,i)))#1,((((f,#5))#6,d),(#2,(j,#1)))),#8))#7,a)))))) cost=35
   >: ((((c)#5,e))#2,(#7,(#6,(#4,(#3,(((((((g)#4,((((h)#3,b))#8,i)))#1,(((f,#5))#6,(d,(#2,(j,#1))))),#8))#7,a)))))) cost=34
   >: (((((c)#5,e))#2,(#6,(#4,(#3,(((((((g)#4,((((h)#3,b))#8,i)))#1,(((f,#5))#6,(d,(#2,(j,#1))))),#8))#7,a))))),#7) cost=31
   >: (((((c)#5,e))#2,(#6,(#4,((((g)#4,((((h)#3,b))#8,i)))#1,(#3,((((((f,#5))#6,(d,(#2,(j,#1)))),#8))#7,a)))))),#7) cost=30
   >: (((((c)#5,e))#2,(#6,((((g)#4,((((h)#3,b))#8,i)))#1,(#3,(#4,((((((f,#5))#6,(d,(#2,(j,#1)))),#8))#7,a)))))),#7) cost=29
   >: (((((c)#5,e))#2,(#6,((((g)#4,((((h)#3,b))#8,i)))#1,(#3,(#4,(((((f,#5))#6,((d,#8),(#2,(j,#1)))))#7,a)))))),#7) cost=28
   >: (((((c)#5,e))#2,((h)#3,(#6,((((g)#4,((b)#8,i)))#1,(#3,(#4,(((((f,#5))#6,((d,#8),(#2,(j,#1)))))#7,a))))))),#7) cost=27
   >: (((e)#2,((h)#3,(#6,((((g)#4,((b)#8,i)))#1,(#3,(#4,((c)#5,(((((f,#5))#6,((d,#8),(#2,(j,#1)))))#7,a)))))))),#7) cost=26
   >: (((e)#2,((h)#3,(#6,((((g)#4,((b)#8,i)))#1,(#4,((c)#5,(#3,(((((f,#5))#6,((d,#8),(#2,(j,#1)))))#7,a)))))))),#7) cost=25
   >: (((e)#2,((h)#3,(#6,((((b)#8,i))#1,(#4,((c)#5,(#3,(((((f,#5))#6,((d,#8),(#2,(((g)#4,j),#1)))))#7,a)))))))),#7) cost=24
   >: (((e)#2,((h)#3,((((b)#8,i))#1,(#4,((c)#5,(#6,(#3,(((((f,#5))#6,((d,#8),(#2,(((g)#4,j),#1)))))#7,a)))))))),#7) cost=23
   >: (((e)#2,((h)#3,((((b)#8,i))#1,(#4,((c)#5,(#6,(((((f,#5))#6,((d,#8),(#3,(#2,(((g)#4,j),#1))))))#7,a))))))),#7) cost=22
   >: (((e)#2,(#4,((h)#3,((((b)#8,i))#1,((c)#5,(#6,(((((f,#5))#6,((d,#8),(#3,(#2,(((g)#4,j),#1))))))#7,a))))))),#7) cost=21
   >: (((e)#2,(#4,((((b)#8,i))#1,((c)#5,((h)#3,(#6,(((((f,#5))#6,((d,#8),(#3,(#2,(((g)#4,j),#1))))))#7,a))))))),#7) cost=20
   >: (((e)#2,(#4,((i)#1,((c)#5,((h)#3,(#6,(((((((b)#8,f),#5))#6,((d,#8),(#3,(#2,(((g)#4,j),#1))))))#7,a))))))),#7) cost=19
   >: (((e)#2,(#4,((i)#1,((c)#5,((h)#3,(#6,(((((f,#5))#6,((d,#8),(#3,(#2,(((b)#8,((g)#4,j)),#1))))))#7,a))))))),#7) cost=18
   >: (((e)#2,(#4,((i)#1,((c)#5,((h)#3,(#6,(((((f,#5))#6,((d,#8),(#3,(#2,(((((b)#8,g))#4,j),#1))))))#7,a))))))),#7) cost=17
1.  Steps:0 Climbs:0 New optimal cost: 17
Cost:17 TopNets:180 Steps:83473 Climbs:44 CRuns:1 CTime:40.3006 Naive:83474 Naivetime:18.7321 DTcnt:21369344 Class:TreeChild TimeConsistency:0
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
Cost:0 TopNets:254 Steps:476101 Climbs:2448 CRuns:1000 CTime:8.33668 Naive:477101 Naivetime:0.853947 DTcnt:3727948 Class:TreeChild TimeConsistency:0
```

### Using stopping criterion (--hcstopinit=1000) with quasi-consensus rand networks.

```  
> supnet -g "(a,(b,(c,d))); ((a,b),(c,d))" -q-1 -R3 --HC --hcrunstats --hcstopinit=1000 --randseed 13 | tail -10
   >: ((((((d)#3,b))#1,(#3,((a)#2,c))),#2),#1) cost=1
   >: (((((d)#3,b))#1,(#3,(((a,#1))#2,c))),#2) cost=0
   i: (#3,(((#2,((((d)#1,(#1,b)))#2,c)))#3,a)) cost=3
   >: (#3,(((#2,(#1,((((d)#1,b))#2,c))))#3,a)) cost=1
   >: (#3,((((d)#1,b))#2,(((#2,(#1,c)))#3,a))) cost=0
   i: (#2,(((((d)#3,a))#1,((b)#2,(c,#3))),#1)) cost=1
   >: (#2,((((d)#3,a))#1,(((b,#1))#2,(c,#3)))) cost=0
   i: ((d)#2,((((b,#2))#3,(((#3,c))#1,a)),#1)) cost=2
   >: ((d)#2,(((b)#3,(((#3,(c,#2)))#1,a)),#1)) cost=0
Cost:0 TopNets:254 Steps:267241 Climbs:1267 CRuns:544 CTime:4.55414 Naive:267785 Naivetime:0.453857 DTcnt:2092084 Class:TreeChild TimeConsistency:0
```

### Using time stop criterion in HC in seconds `--hcstoptime=TIME`, i.e., stop a HC climb if there is no improvement after TIME seconds.

```
> supnet -g "((i,c),((a,d),(f,(b,((g,(e,j)),h)))));((i,h),((c,f),((a,(d,e)),(j,(g,b)))));(((f,(b,d)),((j,g),(e,i))),(c,(h,a)));(((f,(i,(j,d))),(c,b)),(((h,a),g),e));((((h,e),((c,f),a)),(d,b)),((g,i),j))" -R8 -q1 -v4 --HC --hcrunstats --hcstoptime=0.5 --randseed 13
HC start: hcusenaive=0 runnaiveleqrt=13 tailmove=1
   i: ((((c)#5,e))#2,(#7,(#3,(((((((g)#4,(#4,((((h)#3,((((f,#5))#6,d),b)))#8,i))))#1,(#6,(#2,(j,#1)))),#8))#7,a)))) cost=46
   >: ((((c)#5,e))#2,(#7,((#3,(((((((g)#4,(#4,((((h)#3,(((f)#6,d),b)))#8,i))))#1,(#6,(#2,(j,#1)))),#8))#7,a)),#5))) cost=45
   >: ((((c)#5,e))#2,(#7,((#4,(#3,(((((((g)#4,((((h)#3,(((f)#6,d),b)))#8,i)))#1,(#6,(#2,(j,#1)))),#8))#7,a))),#5))) cost=42
   >: ((((c)#5,e))#2,(#7,((#6,(#4,(#3,(((((((g)#4,((((h)#3,(((f)#6,d),b)))#8,i)))#1,(#2,(j,#1))),#8))#7,a)))),#5))) cost=37
   >: ((((c)#5,e))#2,(#7,(#6,(#4,(#3,(((((((g)#4,((((h)#3,((((f,#5))#6,d),b)))#8,i)))#1,(#2,(j,#1))),#8))#7,a)))))) cost=36
   >: ((((c)#5,e))#2,(#7,(#6,(#4,(#3,(((((((g)#4,((((h)#3,b))#8,i)))#1,((((f,#5))#6,d),(#2,(j,#1)))),#8))#7,a)))))) cost=35
   >: ((((c)#5,e))#2,(#7,(#6,(#4,(#3,(((((((g)#4,((((h)#3,b))#8,i)))#1,(((f,#5))#6,(d,(#2,(j,#1))))),#8))#7,a)))))) cost=34
   >: (((((c)#5,e))#2,(#6,(#4,(#3,(((((((g)#4,((((h)#3,b))#8,i)))#1,(((f,#5))#6,(d,(#2,(j,#1))))),#8))#7,a))))),#7) cost=31
   >: (((((c)#5,e))#2,(#6,(#4,((((g)#4,((((h)#3,b))#8,i)))#1,(#3,((((((f,#5))#6,(d,(#2,(j,#1)))),#8))#7,a)))))),#7) cost=30
   >: (((((c)#5,e))#2,(#6,((((g)#4,((((h)#3,b))#8,i)))#1,(#3,(#4,((((((f,#5))#6,(d,(#2,(j,#1)))),#8))#7,a)))))),#7) cost=29
   >: (((((c)#5,e))#2,(#6,((((g)#4,((((h)#3,b))#8,i)))#1,(#3,(#4,(((((f,#5))#6,((d,#8),(#2,(j,#1)))))#7,a)))))),#7) cost=28
   >: (((((c)#5,e))#2,((h)#3,(#6,((((g)#4,((b)#8,i)))#1,(#3,(#4,(((((f,#5))#6,((d,#8),(#2,(j,#1)))))#7,a))))))),#7) cost=27
   >: (((e)#2,((h)#3,(#6,((((g)#4,((b)#8,i)))#1,(#3,(#4,((c)#5,(((((f,#5))#6,((d,#8),(#2,(j,#1)))))#7,a)))))))),#7) cost=26
   >: (((e)#2,((h)#3,(#6,((((g)#4,((b)#8,i)))#1,(#4,((c)#5,(#3,(((((f,#5))#6,((d,#8),(#2,(j,#1)))))#7,a)))))))),#7) cost=25
   >: (((e)#2,((h)#3,(#6,((((b)#8,i))#1,(#4,((c)#5,(#3,(((((f,#5))#6,((d,#8),(#2,(((g)#4,j),#1)))))#7,a)))))))),#7) cost=24
   >: (((e)#2,((h)#3,((((b)#8,i))#1,(#4,((c)#5,(#6,(#3,(((((f,#5))#6,((d,#8),(#2,(((g)#4,j),#1)))))#7,a)))))))),#7) cost=23
   >: (((e)#2,((h)#3,((((b)#8,i))#1,(#4,((c)#5,(#6,(((((f,#5))#6,((d,#8),(#3,(#2,(((g)#4,j),#1))))))#7,a))))))),#7) cost=22
   >: (((e)#2,(#4,((h)#3,((((b)#8,i))#1,((c)#5,(#6,(((((f,#5))#6,((d,#8),(#3,(#2,(((g)#4,j),#1))))))#7,a))))))),#7) cost=21
   >: (((e)#2,(#4,((((b)#8,i))#1,((c)#5,((h)#3,(#6,(((((f,#5))#6,((d,#8),(#3,(#2,(((g)#4,j),#1))))))#7,a))))))),#7) cost=20
   >: (((e)#2,(#4,((i)#1,((c)#5,((h)#3,(#6,(((((((b)#8,f),#5))#6,((d,#8),(#3,(#2,(((g)#4,j),#1))))))#7,a))))))),#7) cost=19
   >: (((e)#2,(#4,((i)#1,((c)#5,((h)#3,(#6,(((((f,#5))#6,((d,#8),(#3,(#2,(((b)#8,((g)#4,j)),#1))))))#7,a))))))),#7) cost=18
   >: (((e)#2,(#4,((i)#1,((c)#5,((h)#3,(#6,(((((f,#5))#6,((d,#8),(#3,(#2,(((((b)#8,g))#4,j),#1))))))#7,a))))))),#7) cost=17
HC run stopped due to timeout, cost=17
1.  Steps:0 Climbs:0 New optimal cost: 17
Stats saved to odt.dat
Best networks saved to odt.log
Cost:17 TopNets:33 Steps:6461 Climbs:44 CRuns:1 CTime:1.53098 Naive:6462 Naivetime:1.39497 DTcnt:1654272 Class:TreeChild TimeConsistency:0
```

### Save odt files after each improvement

Use `--savewhenimproved`. Note that stats may be not reliable, before the end of computation.

```
> supnet -g "((i,c),((a,d),(f,(b,((g,(e,j)),h)))));((i,h),((c,f),((a,(d,e)),(j,(g,b)))));(((f,(b,d)),((j,g),(e,i))),(c,(h,a)));(((f,(i,(j,d))),(c,b)),(((h,a),g),e));((((h,e),((c,f),a)),(d,b)),((g,i),j))" -R8 -q1 -v4 --HC --hcrunstats --hcstoptime=0.5 --randseed 13 --savewhenimproved
HC start: hcusenaive=0 runnaiveleqrt=13 tailmove=1
   i: ((((c)#5,e))#2,(#7,(#3,(((((((g)#4,(#4,((((h)#3,((((f,#5))#6,d),b)))#8,i))))#1,(#6,(#2,(j,#1)))),#8))#7,a)))) cost=46
   >: ((((c)#5,e))#2,(#7,((#3,(((((((g)#4,(#4,((((h)#3,(((f)#6,d),b)))#8,i))))#1,(#6,(#2,(j,#1)))),#8))#7,a)),#5))) cost=45
   >: ((((c)#5,e))#2,(#7,((#4,(#3,(((((((g)#4,((((h)#3,(((f)#6,d),b)))#8,i)))#1,(#6,(#2,(j,#1)))),#8))#7,a))),#5))) cost=42
   >: ((((c)#5,e))#2,(#7,((#6,(#4,(#3,(((((((g)#4,((((h)#3,(((f)#6,d),b)))#8,i)))#1,(#2,(j,#1))),#8))#7,a)))),#5))) cost=37
   >: ((((c)#5,e))#2,(#7,(#6,(#4,(#3,(((((((g)#4,((((h)#3,((((f,#5))#6,d),b)))#8,i)))#1,(#2,(j,#1))),#8))#7,a)))))) cost=36
   >: ((((c)#5,e))#2,(#7,(#6,(#4,(#3,(((((((g)#4,((((h)#3,b))#8,i)))#1,((((f,#5))#6,d),(#2,(j,#1)))),#8))#7,a)))))) cost=35
   >: ((((c)#5,e))#2,(#7,(#6,(#4,(#3,(((((((g)#4,((((h)#3,b))#8,i)))#1,(((f,#5))#6,(d,(#2,(j,#1))))),#8))#7,a)))))) cost=34
   >: (((((c)#5,e))#2,(#6,(#4,(#3,(((((((g)#4,((((h)#3,b))#8,i)))#1,(((f,#5))#6,(d,(#2,(j,#1))))),#8))#7,a))))),#7) cost=31
   >: (((((c)#5,e))#2,(#6,(#4,((((g)#4,((((h)#3,b))#8,i)))#1,(#3,((((((f,#5))#6,(d,(#2,(j,#1)))),#8))#7,a)))))),#7) cost=30
   >: (((((c)#5,e))#2,(#6,((((g)#4,((((h)#3,b))#8,i)))#1,(#3,(#4,((((((f,#5))#6,(d,(#2,(j,#1)))),#8))#7,a)))))),#7) cost=29
   >: (((((c)#5,e))#2,(#6,((((g)#4,((((h)#3,b))#8,i)))#1,(#3,(#4,(((((f,#5))#6,((d,#8),(#2,(j,#1)))))#7,a)))))),#7) cost=28
   >: (((((c)#5,e))#2,((h)#3,(#6,((((g)#4,((b)#8,i)))#1,(#3,(#4,(((((f,#5))#6,((d,#8),(#2,(j,#1)))))#7,a))))))),#7) cost=27
   >: (((e)#2,((h)#3,(#6,((((g)#4,((b)#8,i)))#1,(#3,(#4,((c)#5,(((((f,#5))#6,((d,#8),(#2,(j,#1)))))#7,a)))))))),#7) cost=26
   >: (((e)#2,((h)#3,(#6,((((g)#4,((b)#8,i)))#1,(#4,((c)#5,(#3,(((((f,#5))#6,((d,#8),(#2,(j,#1)))))#7,a)))))))),#7) cost=25
   >: (((e)#2,((h)#3,(#6,((((b)#8,i))#1,(#4,((c)#5,(#3,(((((f,#5))#6,((d,#8),(#2,(((g)#4,j),#1)))))#7,a)))))))),#7) cost=24
   >: (((e)#2,((h)#3,((((b)#8,i))#1,(#4,((c)#5,(#6,(#3,(((((f,#5))#6,((d,#8),(#2,(((g)#4,j),#1)))))#7,a)))))))),#7) cost=23
   >: (((e)#2,((h)#3,((((b)#8,i))#1,(#4,((c)#5,(#6,(((((f,#5))#6,((d,#8),(#3,(#2,(((g)#4,j),#1))))))#7,a))))))),#7) cost=22
   >: (((e)#2,(#4,((h)#3,((((b)#8,i))#1,((c)#5,(#6,(((((f,#5))#6,((d,#8),(#3,(#2,(((g)#4,j),#1))))))#7,a))))))),#7) cost=21
   >: (((e)#2,(#4,((((b)#8,i))#1,((c)#5,((h)#3,(#6,(((((f,#5))#6,((d,#8),(#3,(#2,(((g)#4,j),#1))))))#7,a))))))),#7) cost=20
   >: (((e)#2,(#4,((i)#1,((c)#5,((h)#3,(#6,(((((((b)#8,f),#5))#6,((d,#8),(#3,(#2,(((g)#4,j),#1))))))#7,a))))))),#7) cost=19
   >: (((e)#2,(#4,((i)#1,((c)#5,((h)#3,(#6,(((((f,#5))#6,((d,#8),(#3,(#2,(((b)#8,((g)#4,j)),#1))))))#7,a))))))),#7) cost=18
   >: (((e)#2,(#4,((i)#1,((c)#5,((h)#3,(#6,(((((f,#5))#6,((d,#8),(#3,(#2,(((((b)#8,g))#4,j),#1))))))#7,a))))))),#7) cost=17
HC run stopped due to timeout, cost=17
1.  Steps:0 Climbs:0 New optimal cost: 17
Stats saved to odt.dat
Best networks saved to odt.log
Cost:17 TopNets:32 Steps:6438 Climbs:44 CRuns:1 CTime:1.63616 Naive:6439 Naivetime:1.47265 DTcnt:1648384 Class:TreeChild TimeConsistency:0
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
> supnet -g "((i,c),((a,d),(f,(b,((g,(e,j)),h)))));((i,h),((c,f),((a,(d,e)),(j,(g,b)))));(((f,(b,d)),((j,g),(e,i))),(c,(h,a)));(((f,(i,(j,d))),(c,b)),(((h,a),g),e));((((h,e),((c,f),a)),(d,b)),((g,i),j))" -R8 -q1 -v4 --HC --hcrunstats --hcstoptime=0.5 --randseed 13 --savewhenimproved --displaytreesampling="0.125 0.25 0.5"
HC start: hcusenaive=0 runnaiveleqrt=13 tailmove=1
   i: ((((c)#5,e))#2,(#7,(#3,(((((((g)#4,(#4,((((h)#3,((((f,#5))#6,d),b)))#8,i))))#1,(#6,(#2,(j,#1)))),#8))#7,a)))) cost=48
   >: ((((c)#5,e))#2,(#7,((#3,(((((((g)#4,(#4,((((h)#3,(((f)#6,d),b)))#8,i))))#1,(#6,(#2,(j,#1)))),#8))#7,a)),#5))) cost=47
   >: ((((c)#5,e))#2,(#7,(#4,((#3,(((((((g)#4,((((h)#3,(((f)#6,d),b)))#8,i)))#1,(#6,(#2,(j,#1)))),#8))#7,a)),#5)))) cost=46
   >: ((((c)#5,e))#2,(#7,(((#3,(((((((g)#4,((((h)#3,(((f)#6,d),b)))#8,i)))#1,(#6,(#2,(j,#1)))),#8))#7,a)),#4),#5))) cost=45
   >: ((((c)#5,e))#2,(#7,(((#3,(((((((g)#4,((((h)#3,(((f)#6,d),b)))#8,i)))#1,(#6,((#2,j),#1))),#8))#7,a)),#4),#5))) cost=43
   >: ((((c)#5,e))#2,(#7,(((#6,(#3,(((((((g)#4,((((h)#3,(((f)#6,d),b)))#8,i)))#1,((#2,j),#1)),#8))#7,a))),#4),#5))) cost=42
   >: ((((c)#5,e))#2,(#7,(#4,((#6,(#3,(((((((g)#4,((((h)#3,(((f)#6,d),b)))#8,i)))#1,((#2,j),#1)),#8))#7,a))),#5)))) cost=39
   >: ((((c)#5,e))#2,(#7,(((#6,(#3,(((((((g)#4,((((h)#3,(((f)#6,d),b)))#8,i)))#1,((#2,j),#1)),#8))#7,a))),#4),#5))) cost=38
   >: (#7,(((#6,((((c)#5,e))#2,(#3,(((((((g)#4,((((h)#3,(((f)#6,d),b)))#8,i)))#1,((#2,j),#1)),#8))#7,a)))),#4),#5)) cost=37
   >: (#7,(((#6,((((c)#5,e))#2,(#3,((((((((g)#4,((((h)#3,b))#8,i)))#1,((#2,j),#1)),((f)#6,d)),#8))#7,a)))),#4),#5)) cost=36
   >: (#7,(((#6,((e)#2,(#3,((((((((g)#4,((((h)#3,b))#8,i)))#1,((#2,j),#1)),((((c)#5,f))#6,d)),#8))#7,a)))),#4),#5)) cost=35
   >: (#7,(((#6,((e)#2,(#3,((((((((g)#4,((((h)#3,b))#8,i)))#1,((#2,j),#1)),((((c)#5,f))#6,d)),#8))#7,a)))),#5),#4)) cost=34
   >: ((#7,((#2,j),#1)),((#6,((e)#2,((#3,(((((((g)#4,((((h)#3,b))#8,i)))#1,((((c)#5,f))#6,d)),#8))#7,a)),#5))),#4)) cost=32
HC run with DTsampling=0.125 completed: ((#7,((#2,j),#1)),((#6,((e)#2,((#3,(((((((g)#4,((((h)#3,b))#8,i)))#1,((((c)#5,f))#6,d)),#8))#7,a)),#5))),#4)) cost=32
1.  Steps:0 Climbs:0 New optimal cost: 32
   i: ((#7,((#2,j),#1)),((#6,((e)#2,((#3,(((((((g)#4,((((h)#3,b))#8,i)))#1,((((c)#5,f))#6,d)),#8))#7,a)),#5))),#4)) cost=36
   >: ((#7,((#2,j),#1)),(#6,((e)#2,(((#3,(((((((g)#4,((((h)#3,b))#8,i)))#1,((((c)#5,f))#6,d)),#8))#7,a)),#4),#5)))) cost=28
   >: (((((c)#5,f))#6,(#7,((#2,((h)#3,j)),#1))),(#6,((e)#2,(((#3,(((((((g)#4,((b)#8,i)))#1,d),#8))#7,a)),#4),#5)))) cost=27
...
   >: (((((c)#5,f))#6,((((h)#3,((((b)#8,i))#1,((g)#4,(#2,j)))),#7),#1)),(((e)#2,(#4,(#3,(#6,(((d,#8))#7,a))))),#5)) cost=21
   >: (((((c)#5,f))#6,((((h)#3,((((b)#8,i))#1,(((g)#4,j),#2))),#7),#1)),(((e)#2,(#4,(#3,(#6,(((d,#8))#7,a))))),#5)) cost=20
HC run with DTsampling=0.5 completed: (((((c)#5,f))#6,((((h)#3,((((b)#8,i))#1,(((g)#4,j),#2))),#7),#1)),((e)#2,((#4,(#3,(#6,(((d,#8))#7,a)))),#5))) cost=20
1.  Steps:0 Climbs:0 New optimal cost: 20
   i: (((((c)#5,f))#6,((((h)#3,((((b)#8,i))#1,(((g)#4,j),#2))),#7),#1)),((e)#2,((#4,(#3,(#6,(((d,#8))#7,a)))),#5))) cost=20
   >: (((((c)#5,f))#6,((((h)#3,((((b)#8,i))#1,(((g)#4,j),#2))),#7),#1)),(((d,#8))#7,((e)#2,((#4,(#3,(#6,a))),#5)))) cost=19
   >: (((((c)#5,f))#6,(((h)#3,((b)#8,((i)#1,(((g)#4,(j,#7)),#2)))),#1)),(((d,#8))#7,((e)#2,((#4,(#3,(#6,a))),#5)))) cost=18
   >: (((((c)#5,f))#6,((b)#8,((((h)#3,i))#1,((g)#4,((j,#2),#7))))),((((d,#8))#7,((e)#2,((#4,(#3,(#6,a))),#5))),#1)) cost=17
   >: ((((((c)#5,f))#6,((b)#8,((((h)#3,i))#1,((g)#4,((j,#2),#7))))),(((d,#8))#7,((e)#2,((#4,(#3,(#6,a))),#5)))),#1) cost=16
HC run stopped due to timeout, cost=16
1.  Steps:0 Climbs:0 New optimal cost: 16
   i: (((((c)#5,f))#6,((((h)#3,((((b)#8,i))#1,(((g)#4,j),#2))),#7),#1)),(((e)#2,(#4,(#3,(#6,(((d,#8))#7,a))))),#5)) cost=20
   >: (((((c)#5,f))#6,((((h)#3,((((b)#8,i))#1,(((g)#4,j),#2))),#7),#1)),((((d,#8))#7,((e)#2,(#4,(#3,(#6,a))))),#5)) cost=19
   >: (((((c)#5,f))#6,(((h)#3,((b)#8,(((i)#1,(((g)#4,j),#2)),#7))),#1)),((((d,#8))#7,((e)#2,(#4,(#3,(#6,a))))),#5)) cost=18
   >: (((((c)#5,f))#6,(((h)#3,((b)#8,((((g)#4,j),#2),#7))),#1)),((((d,#8))#7,((((i)#1,e))#2,(#4,(#3,(#6,a))))),#5)) cost=17
HC run stopped due to timeout, cost=17
Stats data saved: odt.dat
Best networks saved to odt.log
Cost:16 TopNets:18 Steps:14278 Climbs:64 CRuns:5 CTime:2.34564 Naive:14283 Naivetime:2.01334 DTcnt:2457082 Class:TreeChild TimeConsistency:0
```

Use `hcsamplingmaxnetstonextlevel=MAXNETS` to limit the number of networks passed to the next level sampler. Default is 0 (unlimited).

To test the effectiveness of display tree sampling employ 
the `--testdisplaytreesampling=RETICULATIONCNT` flag, where 'RETICULATIONCNT' indicates the number of reticulations. Here the number of reticulation was set to 10, i.e., to 1024 display trees. 

```
> supnet --HC --displaytreesampling="0.0625 0.125 0.25 0.5 1" --testdisplaytreesampling=10
 Sampler 0.0625 79
 Sampler 0.125 145
 Sampler 0.25 286
 Sampler 0.5 569
 Sampler 1 790
 Sampler 0 1024
```

When sampling if a display tree already yields a better cost that found in HC, the climb can continue immediately if `--cutwhendtimproved` is set. To be tested.

Use `--hcdetailedsummary` to print summary on each sampling level.
Use `--hcdetailedsummarydat` to write additional summary dat file.


## Output files

Configure the `--autooutfiles` option to generate filenames with the cost value followed by `dat` or `log` In the event of naming conflicts, increment numbers will be added, such as COST.1.dat, COST.2.dat, and so forth. Try with '--outdirectory'.

Output: sp/odt.dat and sp/odt.log 
```
> rm -f sp/* && supnet -g '(a,(b,(c,d)));((a,b),(c,a));((b,c),(d,a))' -q1 -R3 --HC --outdirectory sp && ls sp
   i: (#3,((((d)#1,a))#2,(#2,(((c,#1))#3,b)))) cost=1
   >: (#3,((((a)#2,d))#1,(#2,(((c,#1))#3,b)))) cost=0
Cost:0 TopNets:3 Steps:114 Climbs:4 CRuns:1 CTime:0.000672817 Naive:115 Naivetime:0.000227451 DTcnt:920 Class:TreeChild TimeConsistency:0
odt.dat
odt.log
```

Output: test.dat and test.log
```
> supnet -g '(a,(b,(c,d)));((a,b),(c,a));((b,c),(d,a))' -q1 -R3 --HC --outdirectory sp --outfiles=test && ls sp
   i: ((((d)#1,(((a)#2,b),#2)))#3,(#3,(c,#1))) cost=4
   >: ((((d)#1,((a)#2,b)))#3,((#3,(c,#1)),#2)) cost=3
   >: ((((d)#1,((a)#2,b)))#3,(((c,#3),#1),#2)) cost=2
   >: ((((c,#3),#1),(((d)#1,((a)#2,b)))#3),#2) cost=1
   >: ((((c,#3),#1),(((((d)#1,a))#2,b))#3),#2) cost=0
Cost:0 TopNets:5 Steps:179 Climbs:10 CRuns:1 CTime:0.000797033 Naive:180 Naivetime:0.000341654 DTcnt:1440 Class:TreeChild TimeConsistency:0
odt.dat
odt.log
test.dat
test.log
```

Output: 0.dat and 0.log, where 0 is the resulting optimal cost
```
> supnet -g '(a,(b,(c,d)));((a,b),(c,a));((b,c),(d,a))' -q1 -R3 --HC --outdirectory sp --autooutfiles && ls sp
   i: (#3,((c)#2,(((((#2,d))#3,a))#1,(#1,b)))) cost=2
   >: (#3,(#1,((c)#2,(((((#2,d))#3,a))#1,b)))) cost=1
   >: ((#3,((c)#2,(((((#2,d))#3,a))#1,b))),#1) cost=0
Cost:0 TopNets:2 Steps:67 Climbs:6 CRuns:1 CTime:0.000241756 Naive:68 Naivetime:0.000129461 DTcnt:544 Class:TreeChild TimeConsistency:0
0.dat
0.log
odt.dat
odt.log
test.dat
test.log
```

Output: 0.1.dat and 0.1.log
```
> supnet -g '(a,(b,(c,d)));((a,b),(c,a));((b,c),(d,a))' -q1 -R3 --HC --outdirectory sp --autooutfiles && ls sp
   i: (((#3,b))#1,((c)#2,((a)#3,((d,#2),#1)))) cost=2
   >: (((#3,(b,#2)))#1,((c)#2,((a)#3,(d,#1)))) cost=1
   >: (((#3,(b,#2)))#1,((a)#3,(((c)#2,d),#1))) cost=0
Cost:0 TopNets:3 Steps:106 Climbs:6 CRuns:1 CTime:0.000394344 Naive:107 Naivetime:0.000205517 DTcnt:856 Class:TreeChild TimeConsistency:0
0.1.dat
0.1.log
0.dat
0.log
odt.dat
odt.log
test.dat
test.log
```

Output: 0.2.dat and 0.2.log
```
> supnet -g '(a,(b,(c,d)));((a,b),(c,a));((b,c),(d,a))' -q1 -R3 --HC --outdirectory sp --autooutfiles && ls sp
   i: ((c)#3,((#1,(((((a)#1,d),#3))#2,b)),#2)) cost=2
   >: ((c)#3,((#1,(((((a)#1,d))#2,b),#3)),#2)) cost=1
   >: ((((((a)#1,d))#2,c))#3,((#1,(b,#3)),#2)) cost=0
Cost:0 TopNets:3 Steps:104 Climbs:6 CRuns:1 CTime:0.000464916 Naive:105 Naivetime:0.000216246 DTcnt:840 Class:TreeChild TimeConsistency:0
0.1.dat
0.1.log
0.2.dat
0.2.log
0.dat
0.log
odt.dat
odt.log
test.dat
test.log
```





## More examples

To print the minimum total cost for 10 random gene trees versus random tree-child networks with 5 reticulations over 8 species, and also print the initial network:

```
> supnet -r10 -A8 --pnetworks  | supnet -G- -r1 -A8 -R5 --pnetworks --odtnaivecost
(((#4,g),#5),(#1,(((((((((e)#5,(#3,((h)#3,f))))#2,c),(b,#2)))#1,a))#4,d)))
86
```


## Coronavirus dataset processing

Infering tree-child networks without time-consistency, start from 10 networks `-q10`.
Use: `--savewhenimproved` to save odt files after each improvement
and  `--hcstoptime=30` to stop a single HC if there is no improvement after 30 seconds.

```
supnet -G corona.txt -R7 -q10 --HC --hcrunstats --outfiles corona --hcstoptime=30  --savewhenimproved
```
Results in `corona.log` and `corona.dat`.

Inferring relaxed time-consistent networks, start from 10 networks `-q10`.
```
supnet -G corona.txt -R7 -q10 --HC --hcrunstats --timeconsistent --relaxed --outdirectory corona_tc_rel
```

<!-- Inferring relaxed time-consistent networks, start from 10 networks `-q10`, save networks to #first.log.
```
supnet -G corona.txt -R7 -q10 --HC -v1 --timeconsistent --relaxed --outfiles first --savewhenimproved
supnet -G corona.txt -N first.log --savewhenimproved  
supnet -G corona.txt -N first.log --HC -v1 --timeconsistent --relaxed --savewhenimproved
``` -->

### Parallel processing 

TODO: multithreaded supnet

Use parallel with 10 jobs to run 20 supnet computations.
```
parallel --jobs 10 --ungroup supnet -G corona.txt -R7 -q1 --HC --hcrunstats --outdirectory corona --hcstoptime 1 --savewhenimproved --odtlabelled --autooutfiles ::: {1..20}
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
102,corona/102.dat,102,4.7727,4.7727,2.38419e-07,2,0,0,106,15472,0,1,260,1980544,4.22166,15473,0,0,0,0,251948443
119,corona/119.dat,119,5.2893,5.2893,2.38419e-07,2,0,0,80,10186,0,1,261,1303936,4.52791,10187,0,0,0,0,3894346635
122,corona/122.dat,122,3.55873,3.55873,4.76837e-07,3,0,0,86,8837,0,1,258,1131264,2.75365,8838,0,0,0,0,3894353523
22,corona/22.dat,22,7.05259,7.05259,0,14,0,0,112,32983,0,1,261,4221952,6.17192,32984,0,0,0,0,251951075
23.1,corona/23.1.dat,23,5.14809,5.14809,2.38419e-07,6,0,0,98,17436,0,1,259,2231936,4.42186,17437,0,0,0,0,251942764
23.2,corona/23.2.dat,23,8.18323,8.18323,2.38419e-07,64,0,0,148,24496,0,1,261,3135616,7.07571,24497,0,0,0,0,251942761
23,corona/23.dat,23,3.22308,3.22308,2.38419e-07,14,0,0,92,19389,0,1,130,2481920,3.07349,19390,0,0,0,0,251942765
24.1,corona/24.1.dat,24,3.38207,3.38207,4.76837e-07,69,0,0,102,26871,0,2,125,1491739,2.34371,26873,0,0,0,0,3892887752
24.2,corona/24.2.dat,24,7.1368,7.1368,2.38419e-07,64,0,0,128,26560,0,1,261,3399808,6.06366,26561,0,0,0,0,3894353403
24.3,corona/24.3.dat,24,5.25572,5.25572,2.38419e-07,36,0,0,106,18876,0,1,264,2416256,4.59196,18877,0,0,0,0,251947305
24,corona/24.dat,24,6.76618,6.76618,2.38419e-07,2,0,0,110,25438,0,1,262,3256192,5.81489,25439,0,0,0,0,3894346632
25.1,corona/25.1.dat,25,7.52293,7.52293,2.38419e-07,4,0,0,114,24546,0,1,260,3142016,6.63867,24547,0,0,0,0,3892752804
25.2,corona/25.2.dat,25,4.80912,4.80912,2.38419e-07,8,0,0,136,18157,0,1,196,2324224,4.55608,18158,0,0,0,0,3892756406
25.3,corona/25.3.dat,25,5.77002,5.77002,2.38419e-07,78,0,0,94,20104,0,1,235,2573440,5.18477,20105,0,0,0,0,3892757807
25.4,corona/25.4.dat,25,3.02409,3.02409,2.38419e-07,85,0,0,122,22827,0,2,80,1091661,2.17269,22829,0,0,0,0,3892889004
25.5,corona/25.5.dat,25,4.71099,4.71099,9.53674e-07,20,0,0,92,16743,0,1,202,2143232,4.47507,16744,0,0,0,0,3894346634
25.6,corona/25.6.dat,25,4.99702,4.99702,9.53674e-07,17,0,0,76,23172,0,1,218,2966144,4.50657,23173,0,0,0,0,3894346638
25.7,corona/25.7.dat,25,8.88987,8.88987,2.38419e-07,24,0,0,112,29396,0,1,261,3762816,7.75601,29397,0,0,0,0,3894346637
...
35,corona/35.dat,35,2.4352,2.4352,2.38419e-07,174,0,0,154,55996,0,6,124,912761,1.39763,56002,0,0,0,0,3892867240
36.2,corona/36.2.dat,36,5.88282,5.88282,4.76837e-07,2,0,0,106,13710,0,1,261,1755008,5.09027,13711,0,0,0,0,251942768
36,corona/36.dat,36,4.74186,4.74186,4.76837e-07,171,0,0,122,31319,0,3,144,2110112,3.48464,31322,0,0,0,0,3892890750
39.1,corona/39.1.dat,39,6.09231,6.09231,2.38419e-07,5,0,0,94,15385,0,1,262,1969408,5.24247,15386,0,0,0,0,251942762
39,corona/39.dat,39,0,0,0,0,0,0,45,0,0,0,183,0,0,0,0,0,0,0,3894341681
41.1,corona/41.1.dat,41,5.14744,5.14744,2.38419e-07,2,0,0,110,15571,0,1,259,1993216,4.47464,15572,0,0,0,0,3892758128
42,corona/42.dat,42,2.14001,2.14,2.38419e-07,131,0,0,72,22148,0,2,87,1090153,1.52107,22150,0,0,0,0,3892891265
46,corona/46.dat,46,4.49964,4.49964,2.38419e-07,5,0,0,82,13012,0,1,260,1665664,3.71095,13013,0,0,0,0,3892752808
50.1,corona/50.1.dat,50,5.30467,5.30467,4.76837e-07,8,0,0,102,14107,0,1,262,1805824,4.51396,14108,0,0,0,0,3892752801
51,corona/51.dat,51,5.65186,5.65186,2.38419e-07,2,0,0,94,15355,0,1,263,1965568,4.98923,15356,0,0,0,0,251942766
54,corona/54.dat,54,5.19397,5.19397,2.38419e-07,2,0,0,70,11313,0,1,262,1448192,4.48229,11314,0,0,0,0,251942760
55,corona/55.dat,55,5.05534,5.05534,2.38419e-07,5,0,0,98,15144,0,1,261,1938560,4.1888,15145,0,0,0,0,3894346636
56,corona/56.dat,56,4.56794,4.56794,4.76837e-07,6,0,0,108,9637,0,1,261,1233664,3.69766,9638,0,0,0,0,3894346639
57,corona/57.dat,57,4.76359,4.76359,9.53674e-07,6,0,0,82,10013,0,1,260,1281792,3.95678,10014,0,0,0,0,3894346640
62,corona/62.dat,62,0,0,0,0,0,0,44,0,0,0,217,0,0,0,0,0,0,0,3894341677
63.1,corona/63.1.dat,63,3.36043,3.36043,2.38419e-07,5,0,0,78,8156,0,1,261,1044096,2.77138,8157,0,0,0,0,251942759
65,corona/65.dat,65,4.52661,4.52661,0,3,0,0,78,11950,0,1,262,1529728,3.7903,11951,0,0,0,0,251942755
81,corona/81.dat,81,4.83534,4.83534,2.38419e-07,10,0,0,100,11665,0,1,260,1493248,4.18813,11666,0,0,0,0,251946209
96,corona/96.dat,96,0,0,0,0,0,0,37,0,0,0,172,0,0,0,0,0,0,0,3894341673
```

Print cost + networks with sort:
```
> cat corona/*.log | supnet -N- -G corona.txt --odtcost | sort -k1 -n
22 ((((((((BM4831BGR2008)#7,BatCoVZC45),BatCoVZXC21),((HuWuhan2020,HuItalyTE48362020))#6),#3))#4,(GuangxiPangolinP2V,((((RaTG13,#6))#3,GuangdongPangolin12019),#4))),(#7,(((Rs3367)#1,(((BtCoV2792005)#5,(BtCoV2732005,Rf1)),(#1,(((HKU312,#5))#2,(SARSCoVBJ1824,SARS))))),#2)))
22 ((((((((BM4831BGR2008)#7,BatCoVZC45),BatCoVZXC21),((HuWuhan2020,HuItalyTE48362020))#6),#3))#4,(GuangxiPangolinP2V,((((RaTG13,#6))#3,GuangdongPangolin12019),#4))),(#7,((Rs3367)#1,((((BtCoV2792005)#5,(BtCoV2732005,Rf1)),(#1,(((HKU312,#5))#2,(SARSCoVBJ1824,SARS)))),#2))))
22 ((((((((BM4831BGR2008)#7,BatCoVZC45),BatCoVZXC21),((HuWuhan2020,HuItalyTE48362020))#6),#3))#4,(GuangxiPangolinP2V,((((RaTG13,#6))#3,GuangdongPangolin12019),#4))),(#7,(((Rs3367)#1,(((BtCoV2792005)#5,(BtCoV2732005,Rf1)),(((HKU312,#5))#2,(#1,(SARSCoVBJ1824,SARS))))),#2)))
22 ((((((((BM4831BGR2008)#7,BatCoVZC45),BatCoVZXC21),((HuWuhan2020,HuItalyTE48362020))#6),#3))#4,(GuangxiPangolinP2V,((((RaTG13,#6))#3,GuangdongPangolin12019),#4))),(#7,((Rs3367)#1,((((BtCoV2792005)#5,(BtCoV2732005,Rf1)),(((HKU312,#5))#2,(#1,(SARSCoVBJ1824,SARS)))),#2))))
22 ((((((((BM4831BGR2008)#7,BatCoVZC45),BatCoVZXC21),((HuWuhan2020,HuItalyTE48362020))#6),#3))#4,(GuangxiPangolinP2V,((((RaTG13,#6))#3,GuangdongPangolin12019),#4))),(#7,((Rs3367)#1,((((BtCoV2792005)#5,(((HKU312,#5))#2,(BtCoV2732005,Rf1))),(#1,(SARSCoVBJ1824,SARS))),#2))))
22 ((((((((BM4831BGR2008)#7,BatCoVZC45),BatCoVZXC21),((HuWuhan2020,HuItalyTE48362020))#6),#3))#4,(GuangxiPangolinP2V,((((RaTG13,#6))#3,GuangdongPangolin12019),#4))),(#7,((Rs3367)#1,(((((HKU312,#5))#2,((BtCoV2792005)#5,(BtCoV2732005,Rf1))),(#1,(SARSCoVBJ1824,SARS))),#2))))
22 ((((((HuWuhan2020,HuItalyTE48362020))#6,((((BM4831BGR2008)#7,BatCoVZC45),BatCoVZXC21),#3)))#4,(GuangxiPangolinP2V,((((RaTG13,#6))#3,GuangdongPangolin12019),#4))),(#7,(((Rs3367)#1,((BtCoV2732005,Rf1),((BtCoV2792005)#5,(((HKU312,#5))#2,(#1,(SARSCoVBJ1824,SARS)))))),#2)))
22 ((((((HuWuhan2020,HuItalyTE48362020))#6,((((BM4831BGR2008)#7,BatCoVZC45),BatCoVZXC21),#3)))#4,(GuangxiPangolinP2V,((((RaTG13,#6))#3,GuangdongPangolin12019),#4))),(#7,((Rs3367)#1,(((BtCoV2732005,Rf1),((BtCoV2792005)#5,(((HKU312,#5))#2,(#1,(SARSCoVBJ1824,SARS))))),#2))))
22 ((((((HuWuhan2020,HuItalyTE48362020))#6,((((BM4831BGR2008)#7,BatCoVZC45),BatCoVZXC21),#3)))#4,(GuangxiPangolinP2V,((((RaTG13,#6))#3,GuangdongPangolin12019),#4))),(#7,(((Rs3367)#1,(((BtCoV2792005)#5,(BtCoV2732005,Rf1)),(#1,(((HKU312,#5))#2,(SARSCoVBJ1824,SARS))))),#2)))
22 ((((((HuWuhan2020,HuItalyTE48362020))#6,((((BM4831BGR2008)#7,BatCoVZC45),BatCoVZXC21),#3)))#4,(GuangxiPangolinP2V,((((RaTG13,#6))#3,GuangdongPangolin12019),#4))),(#7,((Rs3367)#1,((((BtCoV2792005)#5,(BtCoV2732005,Rf1)),(#1,(((HKU312,#5))#2,(SARSCoVBJ1824,SARS)))),#2))))
22 ((((((HuWuhan2020,HuItalyTE48362020))#6,((((BM4831BGR2008)#7,BatCoVZC45),BatCoVZXC21),#3)))#4,(GuangxiPangolinP2V,((((RaTG13,#6))#3,GuangdongPangolin12019),#4))),(#7,(((Rs3367)#1,(((BtCoV2792005)#5,(BtCoV2732005,Rf1)),(((HKU312,#5))#2,(#1,(SARSCoVBJ1824,SARS))))),#2)))
22 ((((((HuWuhan2020,HuItalyTE48362020))#6,((((BM4831BGR2008)#7,BatCoVZC45),BatCoVZXC21),#3)))#4,(GuangxiPangolinP2V,((((RaTG13,#6))#3,GuangdongPangolin12019),#4))),(#7,((Rs3367)#1,((((BtCoV2792005)#5,(BtCoV2732005,Rf1)),(((HKU312,#5))#2,(#1,(SARSCoVBJ1824,SARS)))),#2))))
22 ((((((HuWuhan2020,HuItalyTE48362020))#6,((((BM4831BGR2008)#7,BatCoVZC45),BatCoVZXC21),#3)))#4,(GuangxiPangolinP2V,((((RaTG13,#6))#3,GuangdongPangolin12019),#4))),(#7,((Rs3367)#1,((((BtCoV2792005)#5,(((HKU312,#5))#2,(BtCoV2732005,Rf1))),(#1,(SARSCoVBJ1824,SARS))),#2))))
22 ((((((HuWuhan2020,HuItalyTE48362020))#6,((((BM4831BGR2008)#7,BatCoVZC45),BatCoVZXC21),#3)))#4,(GuangxiPangolinP2V,((((RaTG13,#6))#3,GuangdongPangolin12019),#4))),(#7,((Rs3367)#1,(((((HKU312,#5))#2,((BtCoV2792005)#5,(BtCoV2732005,Rf1))),(#1,(SARSCoVBJ1824,SARS))),#2))))
23 ((((#4,(((#3,(((SARSCoVBJ1824,SARS))#5,BtCoV2792005)))#6,(Rs3367,#5))),((BtCoV2732005,Rf1))#3),#7),(#2,(((((RaTG13,((HuItalyTE48362020,HuWuhan2020))#1),(((((((#6,HKU312))#4,BM4831BGR2008))#7,BatCoVZC45),BatCoVZXC21))#2),#1),GuangdongPangolin12019),GuangxiPangolinP2V)))
23 ((((#4,(((((SARSCoVBJ1824,SARS))#5,(#3,BtCoV2792005)))#6,(Rs3367,#5))),((BtCoV2732005,Rf1))#3),#7),(#2,(((((RaTG13,((HuItalyTE48362020,HuWuhan2020))#1),(((((((#6,HKU312))#4,BM4831BGR2008))#7,BatCoVZC45),BatCoVZXC21))#2),#1),GuangdongPangolin12019),GuangxiPangolinP2V)))
23 (((#5,((#1,((((HuWuhan2020,HuItalyTE48362020))#1,RaTG13),(BatCoVZXC21,(#3,BatCoVZC45)))),((GuangxiPangolinP2V)#5,GuangdongPangolin12019))))#7,((BM4831BGR2008)#3,(((((#2,Rs3367),(((SARS,SARSCoVBJ1824))#2,(((BtCoV2732005,Rf1))#6,(#4,BtCoV2792005)))),#6),#7),(HKU312)#4)))
23 (((#5,((#1,((((HuWuhan2020,HuItalyTE48362020))#1,RaTG13),(BatCoVZXC21,(#3,BatCoVZC45)))),((GuangxiPangolinP2V)#5,GuangdongPangolin12019))))#7,((BM4831BGR2008)#3,(#6,((#7,((#2,Rs3367),(((SARS,SARSCoVBJ1824))#2,(((BtCoV2732005,Rf1))#6,(#4,BtCoV2792005))))),(HKU312)#4))))
23 (((#5,((#1,((((HuWuhan2020,HuItalyTE48362020))#1,RaTG13),(BatCoVZXC21,(#3,BatCoVZC45)))),((GuangxiPangolinP2V)#5,GuangdongPangolin12019))))#7,((BM4831BGR2008)#3,(#7,((((#2,Rs3367),#6),(((SARS,SARSCoVBJ1824))#2,(((BtCoV2732005,Rf1))#6,(#4,BtCoV2792005)))),(HKU312)#4))))
...
65 (#5,(#6,((((((((HuWuhan2020)#4,HKU312),(BtCoV2792005,#1)))#3,((((((GuangxiPangolinP2V)#5,(#2,BatCoVZXC21)),(BM4831BGR2008)#7),(GuangdongPangolin12019,(RaTG13,(#4,HuItalyTE48362020)))))#1,(Rf1,((BatCoVZC45)#2,BtCoV2732005)))))#6,(((SARS,SARSCoVBJ1824),Rs3367),#3)),#7)))
81 (#3,(((#2,(#5,(BM4831BGR2008,(#1,((GuangxiPangolinP2V)#6,((Rf1)#4,HuWuhan2020)))))),#6),((((BatCoVZC45,BatCoVZXC21))#5,(((HuItalyTE48362020)#1,((SARS,SARSCoVBJ1824),#7)),((#4,BtCoV2732005),(((((RaTG13)#3,GuangdongPangolin12019))#2,HKU312),BtCoV2792005)))),(Rs3367)#7)))
81 (#3,(((#2,(#5,(BM4831BGR2008,((GuangxiPangolinP2V)#6,(#1,((Rf1)#4,HuWuhan2020)))))),#6),((((BatCoVZC45,BatCoVZXC21))#5,(((HuItalyTE48362020)#1,((SARS,SARSCoVBJ1824),#7)),((#4,BtCoV2732005),(((((RaTG13)#3,GuangdongPangolin12019))#2,HKU312),BtCoV2792005)))),(Rs3367)#7)))
81 (#3,(((#2,(#5,(BM4831BGR2008,((GuangxiPangolinP2V)#6,((Rf1)#4,(#1,HuWuhan2020)))))),#6),((((BatCoVZC45,BatCoVZXC21))#5,(((HuItalyTE48362020)#1,((SARS,SARSCoVBJ1824),#7)),((#4,BtCoV2732005),(((((RaTG13)#3,GuangdongPangolin12019))#2,BtCoV2792005),HKU312)))),(Rs3367)#7)))
81 (#3,(((#2,(#5,(BM4831BGR2008,((GuangxiPangolinP2V)#6,((Rf1)#4,(#1,HuWuhan2020)))))),#6),((((BatCoVZC45,BatCoVZXC21))#5,(((HuItalyTE48362020)#1,((SARS,SARSCoVBJ1824),#7)),((#4,BtCoV2732005),(((((RaTG13)#3,GuangdongPangolin12019))#2,HKU312),BtCoV2792005)))),(Rs3367)#7)))
81 (#3,(((#2,(#5,(BM4831BGR2008,((GuangxiPangolinP2V)#6,((Rf1)#4,(#1,HuWuhan2020)))))),#6),((((BatCoVZC45,BatCoVZXC21))#5,(((HuItalyTE48362020)#1,((SARS,SARSCoVBJ1824),#7)),((#4,BtCoV2732005),((((RaTG13)#3,GuangdongPangolin12019))#2,(HKU312,BtCoV2792005))))),(Rs3367)#7)))
81 (#3,(((#2,(#5,(BM4831BGR2008,((GuangxiPangolinP2V)#6,((Rf1)#4,(#1,HuWuhan2020)))))),#6),(((BatCoVZC45,BatCoVZXC21))#5,((((HuItalyTE48362020)#1,((SARS,SARSCoVBJ1824),#7)),((#4,BtCoV2732005),(((((RaTG13)#3,GuangdongPangolin12019))#2,HKU312),BtCoV2792005))),(Rs3367)#7))))
81 (#3,(((#2,(#5,(BM4831BGR2008,((GuangxiPangolinP2V)#6,((Rf1)#4,(#1,HuWuhan2020)))))),#6),((((BatCoVZC45,BatCoVZXC21))#5,((((SARS,SARSCoVBJ1824),(HuItalyTE48362020)#1),#7),((#4,BtCoV2732005),(((((RaTG13)#3,GuangdongPangolin12019))#2,HKU312),BtCoV2792005)))),(Rs3367)#7)))
81 (#3,(((#2,(#5,(BM4831BGR2008,((GuangxiPangolinP2V)#6,(((Rf1,BtCoV2732005))#4,(#1,HuWuhan2020)))))),#6),((((BatCoVZC45,BatCoVZXC21))#5,(((HuItalyTE48362020)#1,((SARS,SARSCoVBJ1824),#7)),(#4,(((((RaTG13)#3,GuangdongPangolin12019))#2,HKU312),BtCoV2792005)))),(Rs3367)#7)))
81 (#3,(#6,((#2,(#5,(BM4831BGR2008,((GuangxiPangolinP2V)#6,((Rf1)#4,(#1,HuWuhan2020)))))),((((BatCoVZC45,BatCoVZXC21))#5,(((HuItalyTE48362020)#1,((SARS,SARSCoVBJ1824),#7)),((#4,BtCoV2732005),(((((RaTG13)#3,GuangdongPangolin12019))#2,HKU312),BtCoV2792005)))),(Rs3367)#7))))
81 (#3,(#6,((((#2,(#5,(BM4831BGR2008,((GuangxiPangolinP2V)#6,((Rf1)#4,(#1,HuWuhan2020)))))),((BatCoVZC45,BatCoVZXC21))#5),(((HuItalyTE48362020)#1,((SARS,SARSCoVBJ1824),#7)),((#4,BtCoV2732005),(((((RaTG13)#3,GuangdongPangolin12019))#2,HKU312),BtCoV2792005)))),(Rs3367)#7)))
96 ((((#1,((HKU312)#4,Rf1)))#2,(((((Rs3367)#1,BtCoV2732005),(((SARSCoVBJ1824,SARS))#5,(#6,(#3,GuangdongPangolin12019)))))#7,(BtCoV2792005,#4))),(((((HuWuhan2020,HuItalyTE48362020),RaTG13),(BatCoVZC45,BatCoVZXC21)))#3,((GuangxiPangolinP2V)#6,(#5,(#2,(BM4831BGR2008,#7))))))
102 ((#6,(((((#5,(((((HuItalyTE48362020,HuWuhan2020))#3,(((Rs3367,#3))#4,(BtCoV2732005,Rf1))))#2,(((((GuangdongPangolin12019,(BatCoVZXC21,BatCoVZC45)))#1,BM4831BGR2008))#6,GuangxiPangolinP2V))),#4),(((HKU312,BtCoV2792005),(SARS,SARSCoVBJ1824)))#5))#7,(#1,(#2,RaTG13)))),#7)
102 ((#6,(((((#5,(((((HuItalyTE48362020,HuWuhan2020))#3,((((Rs3367,(SARS,SARSCoVBJ1824)),#3))#4,(BtCoV2732005,Rf1))))#2,(((((GuangdongPangolin12019,(BatCoVZXC21,BatCoVZC45)))#1,BM4831BGR2008))#6,GuangxiPangolinP2V))),#4),((HKU312,BtCoV2792005))#5))#7,(#1,(#2,RaTG13)))),#7)
119 ((GuangxiPangolinP2V,(((GuangdongPangolin12019)#3,(((((((#7,BM4831BGR2008))#4,HuItalyTE48362020))#5,(SARSCoVBJ1824,SARS)))#1,HuWuhan2020)))#6),((#2,(((BatCoVZC45,BatCoVZXC21),#5),(#4,(((((Rf1,BtCoV2732005))#7,RaTG13),(((HKU312,BtCoV2792005))#2,(#1,Rs3367))),#6)))),#3))
119 ((GuangxiPangolinP2V,(((GuangdongPangolin12019)#3,(((((((#7,BM4831BGR2008))#4,HuItalyTE48362020))#5,(SARSCoVBJ1824,SARS)))#1,HuWuhan2020)))#6),((#2,(((BatCoVZC45,BatCoVZXC21),#5),(#4,(((((Rf1,BtCoV2732005))#7,RaTG13),((HKU312,(BtCoV2792005)#2),(#1,Rs3367))),#6)))),#3))
122 ((#1,(#6,((#5,(((SARSCoVBJ1824,SARS))#3,(((((((BatCoVZC45)#4,BM4831BGR2008))#7,((GuangxiPangolinP2V)#6,((HuItalyTE48362020,(HuWuhan2020)#1),RaTG13))),Rf1),(GuangdongPangolin12019)#2),((((#4,BatCoVZXC21))#5,(HKU312,BtCoV2792005)),#2)))),((Rs3367,#3),BtCoV2732005)))),#7)
122 (#1,(#7,(#6,((#5,(((SARSCoVBJ1824,SARS))#3,(((((((BatCoVZC45)#4,BM4831BGR2008))#7,((GuangxiPangolinP2V)#6,(((HuItalyTE48362020,HuWuhan2020))#1,RaTG13))),Rf1),(GuangdongPangolin12019)#2),((((#4,BatCoVZXC21))#5,(HKU312,BtCoV2792005)),#2)))),((Rs3367,#3),BtCoV2732005)))))
122 (#1,(#7,(#6,((#5,(((SARSCoVBJ1824,SARS))#3,(((((((BatCoVZC45)#4,BM4831BGR2008))#7,((GuangxiPangolinP2V)#6,((HuItalyTE48362020,(HuWuhan2020)#1),RaTG13))),Rf1),(GuangdongPangolin12019)#2),((((#4,BatCoVZXC21))#5,(HKU312,BtCoV2792005)),#2)))),((Rs3367,#3),BtCoV2732005)))))
```


## Wheat dataset procesing

Using sampling and output directory wheat_guided. Single run with 10 initial networks.
```
supnet -G wheat_trees_clean -R6 -q10 -v4 --HC --hcstoptime=4 --hcrunstatsext --relaxed --savewhenimproved --odtlabelled --guidetree "(1,2,3,4);(5,6,7,8);(9,10,11,12);(13,14,15);(16,17,18);(19,20,21);(22,23,24);(25,26);(27,28,29);(30,31,32,33);(34,35);(36,37,38,39);(40,41,42,43)" --outdirectory wheat_guided --displaytreesampling="0.001 0.01 0.03 0.05 0.1" --hcsamplerstats --hcsamplingmaxnetstonextlevel=1  --hcdetailedsummary  --autooutfiles
```   

In parallel with 10 cores and 20 runs:
```
parallel -q --jobs 10 --ungroup supnet -G wheat_trees_clean -R6 -q1 -v4 --HC --hcstoptime=4 --hcrunstatsext --relaxed --savewhenimproved --odtlabelled --guidetree "(1,2,3,4);(5,6,7,8);(9,10,11,12);(13,14,15);(16,17,18);(19,20,21);(22,23,24);(25,26);(27,28,29);(30,31,32,33);(34,35);(36,37,38,39);(40,41,42,43)" --outdirectory wheat_guided --displaytreesampling="0.001 0.01 0.03 0.05 0.1" --hcsamplerstats --hcsamplingmaxnetstonextlevel=1  --hcdetailedsummary  --autooutfiles --outfile={1} ::: {1..20}
```

Results in wheat_guided directory.
```
cat wheat_guided/*.log | supnet -G wheat_trees_clean -N- --bestnetworks -v4
cat odt.log
```


### Parallel processing from mulitple dir files


