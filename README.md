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
> supnet --HC -g '(a,(b,c));(b,(a,c))'; supnet -N odt.log --pdisplaytrees | sort | uniq -c | sort -k1 -n
Cost:0 Steps:0 Climbs:0 TopNetworks:0 Class:TreeChild TimeConsistency:0 HCruns:0 HCTime:0
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
> supnet -r10 -R2 -A5 --general --detectclass  -z121
0 0 0 (#1,(((c)#2,((b,d))#1),(#2,(e,a))))
1 0 1 ((((e)#1)#2,a),(c,(#2,(b,(#1,d)))))
0 1 1 (#1,((((c,e))#1,(((d)#2,b),a)),#2))
1 0 1 (((#2,((a)#1,e)),(c,b)),((#1)#2,d))
0 1 1 ((#2,b),((a)#1,(((c,#1))#2,(e,d))))
0 0 1 (((a)#2,((#2)#1,(e,((c,d),b)))),#1)
1 0 0 (((c)#2,(b)#1),((e,(#1,#2)),(a,d)))
0 0 0 ((((((d)#2,a))#1,c),(#1,#2)),(e,b))
0 1 1 (((#1,d))#2,(c,((b)#1,(a,(#2,e)))))
0 1 1 ((c,((((b)#2,a))#1,d)),(#2,(#1,e)))
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

Print 2 quasi consensus treeS with preserved split of the root.
```
> supnet -g '(a,((b,c),d));(a,(b,e))' -s'((a,b),(c,(d,e)))' -q2 --preserveroot --pnetworks
(((e,d),c),(a,b))
((a,b),((e,c),d))
```


### Random trees
    
Print two random networks with no reticulations (i.e., species trees) using species {a,b,c,d,e} (`-A5`). 
```
>  supnet -A5 -r2 --pnetworks
((a,e),((d,b),c))
((c,(b,e)),(d,a))
```

Print two quasi-consensus trees.
```
> supnet -g '(a,((b,c),d));(a,(b,d))' -q2 --pnetworks
((a,d),(b,c))
((d,(c,b)),a)
```

### Random networks 

Print one random relaxed network with species {a,b,c}  and 5 reticulations.
```
>  supnet -A3 -R5 -r1 --pnetworks --relaxed
((#3,((#2)#5,(((#1)#4,((b)#1,c)),#5))),(((#4)#2)#3,a))
```

Print one quasi-consensus network with two reticulations.
```
> supnet -G examples/gtrees.txt -q1 -R2 --pnetwork
((((b)#1,(d,a)))#2,((#1,c),#2))
```

### Fixed seed for random generator

Use `-z SEED` or `--randseed SEED` to set seed for random generator (srand). 

### Inserting reticulations

Insert 2 reticulations into a network; tree-child network in output (default).
```  
> supnet -R2 -n '(a,((d)#1,(b,(c,#1))))' --pnetworks
(#3,((d)#1,((((b)#2,a))#3,(#2,(c,#1)))))
```

Insert 8 reticulations into a network; general network `--general`
```  
> supnet -R8 -n '(a,((d)#1,(b,(c,#1))))' --pnetworks --general 
((((#3)#1,#4))#2,((#5,(((b,((#7,(((((#6)#9,d))#3)#5,#8)),#1)))#4,a)),((#9)#8,(((c)#6)#7,#2))))
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
 - ` --noodtfiles`: do not save odt.log and odt.dat with optimal networks
 - ` --hcusenaive`: use only ODT naive algorithm in HCH (no BB)
  - `--hcstopinit=NUM`: stop when there is no new optimal network after HC with NUM initial networks
 - ` --hcstopclimb=NUM`: stop in HC when there is no improvements after NUM steps
 - ` --hcmaximprovements=NUM`: stop after NUM improvements in HC climb
 - ` -O ODTFILE`: the name of odt.log 
 - ` -D DATFILE`: the name of odt.dat; see also --odtlabelled and --noodtfiles 
 - ` --noodtfiles`: do not generate odt and dat files
 - ` --odtlabelled`: odt and dat file in labelled format

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
         cost := cost'
         N := N'
         continue
     if cost' = cost: 
         store N'

ODTcost(GSet, N): 
    if `--hcusenaive` is set or N has less than '-t THRESHOLD` reticulations:
       compute the cost by ODT naive algorithm, i.e., enumeration of all possible display trees)
    otherwise:        
       compute the cost by BB algorithm (using multiple calls of DP and/or ODT naive)
```

By default each HC step starts from the set of predefined networks (see `-n`, `-N` or random networks); otherwise by using quasi-random networks.


Minimalistic hill climbing (HC) run using Tail Moves (default) with a single random network (and one HC run) `-r1` with two reticulations `-R2`. Print cost and save result to odt.log. 

```
> supnet -g "(a,(b,(c,d))); ((a,b),(c,d))" -r1 -R2 --HC
   i: ((#1,((((c)#1,d))#2,(b,a))),#2) cost=1
   >: ((#1,((b,(((c)#1,d))#2),a)),#2) cost=0
Cost:0 Steps:46 Climbs:2 TopNetworks:6 Class:TreeChild TimeConsistency:0 HCruns:1 HCTime:0.000201702 Naive:47 Naivetime:4.41074e-05 DTcnt:186
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
> supnet -g "(a,(b,(c,d))); ((a,b),(c,d)); ((a,c),(b,d))" -r1 -R1 --HC -v1 -z 1001
   i: ((d)#1,(c,((#1,b),a))) cost=5
   >: ((d)#1,((#1,c),(b,a))) cost=3
   >: ((b,(d)#1),((#1,c),a)) cost=2
Cost:2 Steps:51 Climbs:3 TopNetworks:3 Class:TreeChild TimeConsistency:0 HCruns:1 HCTime:0.000171661 Naive:52 Naivetime:3.67165e-05 DTcnt:104
```

Print improvements and equal cost networks `-v2`:
```  
> supnet -g "(a,(b,(c,d))); ((a,b),(c,d)); ((a,c),(b,d))" -r1 -R1 --HC -v2 -z 1001
   i: ((d)#1,(c,((#1,b),a))) cost=5
   >: ((d)#1,((#1,c),(b,a))) cost=3
   =: (((d)#1,c),(#1,(b,a))) cost=3
   =: (((d,c))#1,(#1,(b,a))) cost=3
   =: ((d)#1,((b,(#1,c)),a)) cost=3
   >: ((b,(d)#1),((#1,c),a)) cost=2
   =: (((#1,c),(b,(d)#1)),a) cost=2
   =: (((b,(d)#1),a),(#1,c)) cost=2
   =: (((b,(d)#1),(#1,c)),a) cost=2
   =: ((#1,c),((b,(d)#1),a)) cost=2
Cost:2 Steps:51 Climbs:3 TopNetworks:3 Class:TreeChild TimeConsistency:0 HCruns:1 HCTime:0.000185013 Naive:52 Naivetime:3.55244e-05 DTcnt:104
```

Print every visited network `-v3`:
```  
> supnet -g "(a,(b,c)); ((a,b),c);" -r1 -R1 --HC -v3 -z 1001
   i: (((c)#1,(#1,a)),b) cost=2
   >: (((c)#1,a),(#1,b)) cost=1
   <: ((#1,((c)#1,a)),b) cost=2
   <: (((c)#1,(#1,a)),b) cost=2
   =: ((((c)#1,a),b),#1) cost=1
   >: ((((c)#1,b),a),#1) cost=0
   <: (((c)#1,(b,a)),#1) cost=1
   <: ((b,((c)#1,a)),#1) cost=1
   <: ((((c)#1,a),b),#1) cost=1
   <: (((c)#1,a),(#1,b)) cost=1
   <: (((c)#1,(a,b)),#1) cost=1
   <: ((((c,b))#1,a),#1) cost=1
   <: (((c)#1,(b,#1)),a) cost=1
   <: ((((c)#1,b),#1),a) cost=1
   <: (((c)#1,b),(a,#1)) cost=1
   <: (a,(((c)#1,b),#1)) cost=1
   <: ((((c)#1,a),b),#1) cost=1
   <: (((c)#1,(b,a)),#1) cost=1
   <: (((c)#1,b),(#1,a)) cost=1
   <: ((((c,a))#1,b),#1) cost=2
Cost:0 Steps:19 Climbs:3 TopNetworks:1 Class:TreeChild TimeConsistency:0 HCruns:1 HCTime:9.65595e-05 Naive:20 Naivetime:1.0252e-05 DTcnt:40
```

To show details on DP & BB algorithms use `-v` with 4, 5 and 6. See examples below.





## Networks/dags comparing and aggregating

 - `--comparedags`: pairwise comparison of all networks or their shapes
 - `--uniquedags`: print unique networks; summary stats printed on stderr
 - `--uniquedagscnts`: as above plus counts of networks
 - `--dagshapes`: two networks are equal if their shapes are isomorphic in --*dags* options (i.e., ignore leaf labels)


Print unique networks from 100 randomly generated. 

```  
> supnet -A2 -R1 -r100 --pnetworks | supnet --uniquedags -N-
((b)#1,(#1,a))
(#1,((a)#1,b))
unique=2 all=100
```

Pairwise comparison of random general networks. Here, 
 - `1` - dags are equal,
 - `0` - dags are not equal.

Use `--dagshapes`, to compare shapes.

```  
> supnet -r4 -R1 -A2  --pnetworks --general --comparedags
(#1,((a)#1,b))
(#1,((a)#1,b))
(#1,((a)#1,b))
(#1,((a)#1,b))

(#1,((a)#1,b))	(#1,((a)#1,b))	1
(#1,((a)#1,b))	(#1,((a)#1,b))	1
(#1,((a)#1,b))	(#1,((a)#1,b))	1

(#1,((a)#1,b))	(#1,((a)#1,b))	1
(#1,((a)#1,b))	(#1,((a)#1,b))	1

(#1,((a)#1,b))	(#1,((a)#1,b))	1

```

Print unique random networks with counts
```  
> supnet -r100000 -R1 -A2 --uniquedagscnts
49983	(#1,((b)#1,a))
50017	((a)#1,(#1,b))
unique=2 all=0
```

Print unique random networks with counts using uniform draw
```  
> supnet -r100000 -R1 -A2 --uniquedagscnts --uniformedgesampling
50027	(((a)#1,b),#1)
49973	((b)#1,(#1,a))
unique=2 all=0
```

Print unique random shapes of networks with counts
```  
> supnet -r100000 -R1 -A3 --dagshapes --uniquedagscnts
11613	((a)#1,((c,b),#1))
29950	((b,((c)#1,a)),#1)
33670	((#1,((c)#1,a)),b)
11481	(((a,b))#1,(#1,c))
13286	((b,#1),((c)#1,a))
unique=5 all=0
```

Print unique random shapes of networks with counts using uniform draw
```  
> supnet -r100000 -R1 -A3 --dagshapes --uniformedgesampling --uniquedagscnts
14506	((b)#1,((c,a),#1))
28614	(#1,(((a)#1,b),c))
14240	(#1,(((b,c))#1,a))
14094	(((c)#1,a),(#1,b))
28546	(((c)#1,(b,#1)),a)
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

```
make supnet_dtcache
supnet_dtcache -G corona.txt -R7 -q10 --HC --hcrunstats 
```

Infering relaxed time-consistent networks, start from 10 networks `-q10`.
```
supnet_dtcache -G corona.txt -R7 -q10 --HC --hcrunstats --timeconsistent --relaxed
```

Check network classes in the result:
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
> supnet -g "(a,(b,(c,d))); ((a,b),(c,d)); ((a,c),(b,d)); ((a,d),e)" -r8 -R1 --HC -v0 -z 1003
Cost:3 Steps:814 Climbs:35 TopNetworks:1 Class:TreeChild TimeConsistency:0 HCruns:8 HCTime:0.00227785 Naive:822 Naivetime:0.000627518 DTcnt:1644
```

Print stats after completed HC run for each initial network when only when new optimal network is found.
```
> supnet -g "(a,(b,(c,d))); ((a,b),(c,d)); ((a,c),(b,d)); ((a,d),e)" -r8 -R1 --HC --hcrunstats -v0 -z 1003
1. Cost:5 Steps:67 Climbs:4 TopNetworks:9 Class:TreeChild TimeConsistency:0
 New optimal cost: 5
3. Cost:4 Steps:144 Climbs:6 TopNetworks:6 Class:TreeChild TimeConsistency:0
 New optimal cost: 4
7. Cost:3 Steps:129 Climbs:4 TopNetworks:1 Class:TreeChild TimeConsistency:0
 New optimal cost: 3
Cost:3 Steps:814 Climbs:35 TopNetworks:1 Class:TreeChild TimeConsistency:0 HCruns:8 HCTime:0.00239778 Naive:822 Naivetime:0.000628948 DTcnt:1644
```

Print stats after completed HC run for each initial network when a new optimal network including equal networks is found.
```
> supnet -g "(a,(b,(c,d))); ((a,b),(c,d)); ((a,c),(b,d)); ((a,d),e)" -r8 -R1 --HC --hcrunstatsext -v0 -z 1003
1. Cost:5 Steps:67 Climbs:4 TopNetworks:9 Class:TreeChild TimeConsistency:0
 New optimal cost: 5
3. Cost:4 Steps:144 Climbs:6 TopNetworks:6 Class:TreeChild TimeConsistency:0
 New optimal cost: 4
5. Cost:4 Steps:88 Climbs:5 TopNetworks:3 Class:TreeChild TimeConsistency:0
 New optimal networks:3 Total:9
7. Cost:3 Steps:129 Climbs:4 TopNetworks:1 Class:TreeChild TimeConsistency:0
 New optimal cost: 3
Cost:3 Steps:814 Climbs:35 TopNetworks:1 Class:TreeChild TimeConsistency:0 HCruns:8 HCTime:0.00229216 Naive:822 Naivetime:0.00066638 DTcnt:1644
```

Print stats after completed HC run.
```
> supnet -g "(a,(b,(c,d))); ((a,b),(c,d)); ((a,c),(b,d)); ((a,d),e)" -r8 -R1 --HC --hcrunstatsalways -v0 -z 1003
1. Cost:5 Steps:67 Climbs:4 TopNetworks:9 Class:TreeChild TimeConsistency:0
 New optimal cost: 5
2. Cost:9 Steps:55 Climbs:2 TopNetworks:12 Class:TreeChild TimeConsistency:0
3. Cost:4 Steps:144 Climbs:6 TopNetworks:6 Class:TreeChild TimeConsistency:0
 New optimal cost: 4
4. Cost:4 Steps:183 Climbs:8 TopNetworks:6 Class:TreeChild TimeConsistency:0
 New optimal networks:0 Total:6
5. Cost:4 Steps:88 Climbs:5 TopNetworks:3 Class:TreeChild TimeConsistency:0
 New optimal networks:3 Total:9
6. Cost:4 Steps:66 Climbs:2 TopNetworks:3 Class:TreeChild TimeConsistency:0
 New optimal networks:0 Total:9
7. Cost:3 Steps:129 Climbs:4 TopNetworks:1 Class:TreeChild TimeConsistency:0
 New optimal cost: 3
8. Cost:6 Steps:82 Climbs:4 TopNetworks:2 Class:TreeChild TimeConsistency:0
Cost:3 Steps:814 Climbs:35 TopNetworks:1 Class:TreeChild TimeConsistency:0 HCruns:8 HCTime:0.00243497 Naive:822 Naivetime:0.000667572 DTcnt:1644
```

#### Run naive odt computation threshold `-t THRESHOLD`

Run naive odt computation, when the number of reticulations in a network is `< THRESHOLD`; otherwise run BB & DP algorithms.

Higher `-t NUM`, more naive algorithm computations.

```
> supnet -g "(a,(b,(c,d)))" -n "(#1,((b,((d,c))#1),a))" -t2 --HC -v35
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
ODT-cost-naive 1 tree(s);rt=1(-t) cost=1
   <: (#1,(((b,(d)#1),c),a)) cost=1
ODT-cost-naive 1 tree(s);rt=1(-t) cost=2
   <: (#1,((b,(d)#1),(a,c))) cost=2
ODT-cost-naive 1 tree(s);rt=1(-t) cost=1
   <: ((#1,c),((b,(d)#1),a)) cost=1
ODT-cost-naive 1 tree(s);rt=1(-t) cost=3
   <: (#1,(((b,(d)#1),a),c)) cost=3
ODT-cost-naive 1 tree(s);rt=1(-t) cost=1
   <: (#1,((((b,d),c))#1,a)) cost=1
ODT-cost-naive 1 tree(s);rt=1(-t) cost=1
   <: (#1,(((d,(b,c)))#1,a)) cost=1
ODT-cost-naive 1 tree(s);rt=1(-t) cost=1
   <: (#1,(((d,c))#1,(b,a))) cost=1
ODT-cost-naive 1 tree(s);rt=1(-t) cost=0
   =: ((b,#1),(((d,c))#1,a)) cost=0
ODT-cost-naive 1 tree(s);rt=1(-t) cost=1
   <: (#1,(b,(((d,c))#1,a))) cost=1
ODT-cost-naive 1 tree(s);rt=1(-t) cost=0
   =: (#1,(((b,(d,c)))#1,a)) cost=0
ODT-cost-naive 1 tree(s);rt=1(-t) cost=1
   <: (#1,(b,(a,((d,c))#1))) cost=1
ODT-cost-naive 1 tree(s);rt=1(-t) cost=1
   <: (#1,((b,a),((d,c))#1)) cost=1
ODT-cost-naive 1 tree(s);rt=1(-t) cost=0
   =: (((b,((d,c))#1),#1),a) cost=0
ODT-cost-naive 1 tree(s);rt=1(-t) cost=2
   <: (#1,(b,(((d,a),c))#1)) cost=2
ODT-cost-naive 1 tree(s);rt=1(-t) cost=2
   <: (#1,(b,((d,(c,a)))#1)) cost=2
ODT-cost-naive 1 tree(s);rt=1(-t) cost=1
   <: (#1,((b,a),((d,c))#1)) cost=1
ODT-cost-naive 1 tree(s);rt=1(-t) cost=1
   <: (#1,(b,(((d,c))#1,a))) cost=1
ODT-cost-naive 1 tree(s);rt=1(-t) cost=0
   =: ((#1,a),(b,((d,c))#1)) cost=0
ODT-cost-naive 1 tree(s);rt=1(-t) cost=1
   <: (#1,(b,(((d,c),a))#1)) cost=1
ODT-cost-naive 1 tree(s);rt=1(-t) cost=0
   =: (((#1,b),((d,c))#1),a) cost=0
ODT-cost-naive 1 tree(s);rt=1(-t) cost=0
   =: ((#1,(b,((d,c))#1)),a) cost=0
ODT-cost-naive 1 tree(s);rt=1(-t) cost=0
   =: ((b,((d,c))#1),(#1,a)) cost=0
HC run completed: (#1,((b,((d,c))#1),a)) cost=0
Cost:0 Steps:30 Climbs:1 TopNetworks:6 Class:TreeChild TimeConsistency:0 HCruns:1 HCTime:0.000174999 Naive:31 Naivetime:4.52995e-05 DTcnt:57
Optimal networks saved: odt.log
Stats data save to: odt.dat
```

Lower `-t`, more BB & DP computations
```
> supnet -g "(a,(b,(c,d)))" -n "(#1,((b,((d,c))#1),a))" -t1 --HC -v35
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
ODT-cost-naive 1 tree(s);rt=1(-t) cost=1
   <: (#1,(((b,(d)#1),c),a)) cost=1
ODT-cost-naive 1 tree(s);rt=1(-t) cost=2
   <: (#1,((b,(d)#1),(a,c))) cost=2
ODT-cost-naive 1 tree(s);rt=1(-t) cost=1
   <: ((#1,c),((b,(d)#1),a)) cost=1
ODT-cost-naive 1 tree(s);rt=1(-t) cost=3
   <: (#1,(((b,(d)#1),a),c)) cost=3
ODT-cost-naive 1 tree(s);rt=1(-t) cost=1
   <: (#1,((((b,d),c))#1,a)) cost=1
ODT-cost-naive 1 tree(s);rt=1(-t) cost=1
   <: (#1,(((d,(b,c)))#1,a)) cost=1
ODT-cost-naive 1 tree(s);rt=1(-t) cost=1
   <: (#1,(((d,c))#1,(b,a))) cost=1
ODT-cost-naive 1 tree(s);rt=1(-t) cost=0
   =: ((b,#1),(((d,c))#1,a)) cost=0
ODT-cost-naive 1 tree(s);rt=1(-t) cost=1
   <: (#1,(b,(((d,c))#1,a))) cost=1
ODT-cost-naive 1 tree(s);rt=1(-t) cost=0
   =: (#1,(((b,(d,c)))#1,a)) cost=0
ODT-cost-naive 1 tree(s);rt=1(-t) cost=1
   <: (#1,(b,(a,((d,c))#1))) cost=1
ODT-cost-naive 1 tree(s);rt=1(-t) cost=1
   <: (#1,((b,a),((d,c))#1)) cost=1
ODT-cost-naive 1 tree(s);rt=1(-t) cost=0
   =: (((b,((d,c))#1),#1),a) cost=0
ODT-cost-naive 1 tree(s);rt=1(-t) cost=2
   <: (#1,(b,(((d,a),c))#1)) cost=2
ODT-cost-naive 1 tree(s);rt=1(-t) cost=2
   <: (#1,(b,((d,(c,a)))#1)) cost=2
ODT-cost-naive 1 tree(s);rt=1(-t) cost=1
   <: (#1,((b,a),((d,c))#1)) cost=1
ODT-cost-naive 1 tree(s);rt=1(-t) cost=1
   <: (#1,(b,(((d,c))#1,a))) cost=1
ODT-cost-naive 1 tree(s);rt=1(-t) cost=0
   =: ((#1,a),(b,((d,c))#1)) cost=0
ODT-cost-naive 1 tree(s);rt=1(-t) cost=1
   <: (#1,(b,(((d,c),a))#1)) cost=1
ODT-cost-naive 1 tree(s);rt=1(-t) cost=0
   =: (((#1,b),((d,c))#1),a) cost=0
ODT-cost-naive 1 tree(s);rt=1(-t) cost=0
   =: ((#1,(b,((d,c))#1)),a) cost=0
ODT-cost-naive 1 tree(s);rt=1(-t) cost=0
   =: ((b,((d,c))#1),(#1,a)) cost=0
HC run completed: (#1,((b,((d,c))#1),a)) cost=0
Cost:0 Steps:30 Climbs:1 TopNetworks:6 Class:TreeChild TimeConsistency:0 HCruns:1 HCTime:0.000142813 Naive:31 Naivetime:1.45435e-05 DTcnt:57
Optimal networks saved: odt.log
Stats data save to: odt.dat
```

### NNI vs Tail Moves

HC via Tail Moves (default)
```
> supnet -g "(a,(b,(c,d))); ((a,b),(c,d)); ((a,c),(b,d))" -r1 -R1 --HC -v3 -z 1003
   i: ((#1,d),(((c)#1,b),a)) cost=3
   <: (d,((#1,((c)#1,b)),a)) cost=7
   <: (d,(((c)#1,b),(#1,a))) cost=6
   <: (d,(#1,(((c)#1,b),a))) cost=6
   <: (d,(((c)#1,(#1,b)),a)) cost=7
   <: (#1,((((c)#1,b),d),a)) cost=4
   <: (#1,(((c)#1,b),(a,d))) cost=6
   <: (#1,((((c)#1,b),a),d)) cost=6
   >: (#1,((((c)#1,d),b),a)) cost=2
   <: (#1,(((c)#1,(d,b)),a)) cost=4
   <: (#1,((d,b),((c)#1,a))) cost=4
   <: (#1,((c)#1,((d,b),a))) cost=6
   <: (#1,((d,((c)#1,b)),a)) cost=4
   <: (#1,((((c)#1,b),d),a)) cost=4
   <: (#1,(((c)#1,b),(a,d))) cost=6
   <: ((#1,d),(((c)#1,b),a)) cost=3
   <: (#1,((((c)#1,b),a),d)) cost=6
   <: (#1,(((c)#1,(b,d)),a)) cost=4
   =: (#1,((((c,d))#1,b),a)) cost=2
   =: (((((c)#1,d),b),#1),a) cost=2
   <: (#1,((((c)#1,a),d),b)) cost=5
   <: (#1,(((c)#1,(d,a)),b)) cost=6
   >: ((#1,a),(((c)#1,d),b)) cost=1
   <: ((#1,((c)#1,a)),(d,b)) cost=4
   <: (((c)#1,(#1,a)),(d,b)) cost=4
   <: ((#1,a),((c)#1,(d,b))) cost=3
   <: ((#1,a),(d,((c)#1,b))) cost=3
   <: (((#1,d),a),((c)#1,b)) cost=4
   <: ((#1,(a,d)),((c)#1,b)) cost=6
   <: (((#1,a),d),((c)#1,b)) cost=5
   <: ((#1,a),(((c)#1,b),d)) cost=3
   <: ((#1,a),((c)#1,(b,d))) cost=3
   <: ((#1,a),(((c,d))#1,b)) cost=3
   <: (a,(((c)#1,(#1,d)),b)) cost=3
   <: (a,(#1,(((c)#1,d),b))) cost=2
   <: (a,((#1,((c)#1,d)),b)) cost=3
   <: (a,(((c)#1,d),(#1,b))) cost=3
   <: (#1,((((c)#1,a),d),b)) cost=5
   <: (#1,(((c)#1,(d,a)),b)) cost=6
   <: (#1,((((c)#1,d),b),a)) cost=2
   <: (#1,((((c)#1,d),a),b)) cost=4
   <: (#1,(((c)#1,d),(b,a))) cost=3
   <: (#1,((((c,a))#1,d),b)) cost=4
   <: ((((#1,a),(c)#1),d),b) cost=5
   <: (((c)#1,((#1,a),d)),b) cost=5
   <: (((#1,a),((c)#1,d)),b) cost=3
   <: (((c)#1,d),((#1,a),b)) cost=2
   <: ((#1,(((c)#1,d),b)),a) cost=2
   <: (#1,(a,(((c)#1,d),b))) cost=2
   <: (((((c)#1,d),#1),a),b) cost=4
   <: ((#1,(((c)#1,d),a)),b) cost=4
   <: ((((c)#1,d),(#1,a)),b) cost=3
   <: ((#1,a),(((c)#1,b),d)) cost=3
   <: ((#1,a),((c)#1,(d,b))) cost=3
   <: (((#1,b),a),((c)#1,d)) cost=3
   <: ((#1,(a,b)),((c)#1,d)) cost=3
   <: (((#1,a),b),((c)#1,d)) cost=2
   <: ((#1,a),(((c,b))#1,d)) cost=5
Cost:1 Steps:57 Climbs:3 TopNetworks:1 Class:TreeChild TimeConsistency:0 HCruns:1 HCTime:0.000224352 Naive:58 Naivetime:3.74317e-05 DTcnt:116
```

HC via NNI moves 

```
> supnet -g "(a,(b,(c,d))); ((a,b),(c,d)); ((a,c),(b,d))" -r1 -R1 --HC --hcnnimove -v3 -z 1003
   i: ((#1,d),(((c)#1,b),a)) cost=3
   <: ((d,(((c)#1,b),a)),#1) cost=6
   <: (((((c)#1,b),a),#1),d) cost=6
   <: (((c)#1,b),(a,(#1,d))) cost=4
   =: (a,((#1,d),((c)#1,b))) cost=3
   =: ((#1,d),((b,a),(c)#1)) cost=3
   >: ((#1,d),((a,(c)#1),b)) cost=2
   <: ((d,((a,(c)#1),b)),#1) cost=5
   <: ((((a,(c)#1),b),#1),d) cost=5
   >: ((a,(c)#1),(b,(#1,d))) cost=1
   <: ((a,(c)#1),(#1,(d,b))) cost=3
   <: ((a,(c)#1),(d,(b,#1))) cost=3
   <: (b,((#1,d),(a,(c)#1))) cost=3
   <: ((#1,d),((a,(c)#1),b)) cost=2
   <: (((c)#1,(b,(#1,d))),a) cost=2
   <: (((b,(#1,d)),a),(c)#1) cost=2
Cost:1 Steps:15 Climbs:3 TopNetworks:1 Class:TreeChild TimeConsistency:0 HCruns:1 HCTime:8.55923e-05 Naive:16 Naivetime:1.26362e-05 DTcnt:32
```

#### Default tree-child networks and Tail Moves:
```
> supnet -g "(a,(b,(c,d))); ((a,b),(c,d)); ((a,c),(b,d))" -r1 -R3 --HC -v2 -z 12
   i: ((#3,((c)#3,((((d)#1,(b,#1)))#2,a))),#2) cost=4
   >: (((c)#3,((((d)#1,((#3,b),#1)))#2,a)),#2) cost=3
   =: (((c)#3,((((d)#1,(#3,(b,#1))))#2,a)),#2) cost=3
   >: (((c)#3,(((((#3,d))#1,(b,#1)))#2,a)),#2) cost=1
   >: ((((c)#3,(((((#3,d))#1,b))#2,a)),#1),#2) cost=0
   =: (#1,(((c)#3,(((((#3,d))#1,b))#2,a)),#2)) cost=0
   =: (((c)#3,((((((#3,d))#1,b))#2,a),#1)),#2) cost=0
   =: ((((c)#3,(((((#3,d))#1,b))#2,a)),#2),#1) cost=0
   =: (((((((#3,d))#1,b))#2,((c)#3,a)),#1),#2) cost=0
   =: (((((c)#3,a),((((#3,d))#1,b))#2),#1),#2) cost=0
   =: (((c)#3,((((((#3,d))#1,b))#2,a),#1)),#2) cost=0
   =: (((((((#3,d))#1,b))#2,((c)#3,a)),#1),#2) cost=0
   =: (((c)#3,(#1,(((((#3,d))#1,b))#2,a))),#2) cost=0
Cost:0 Steps:41 Climbs:4 TopNetworks:4 Class:TreeChild TimeConsistency:0 HCruns:1 HCTime:0.000254154 Naive:42 Naivetime:7.55787e-05 DTcnt:318
```

```
> cat odt.log
((((c)#3,(((((#3,d))#1,b))#2,a)),#1),#2)
(#1,(((c)#3,(((((#3,d))#1,b))#2,a)),#2))
(((c)#3,((((((#3,d))#1,b))#2,a),#1)),#2)
(((((((#3,d))#1,b))#2,((c)#3,a)),#1),#2)
```

#### Relaxed networks. 
```
> supnet -g "(a,(b,(c,d))); ((a,b),(c,d)); ((a,c),(b,d))" -r1 -R3 --HC -v2 -z 12 --relaxed
   i: ((#3,((c)#3,(((d)#1)#2,a))),(#2,(b,#1))) cost=4
   >: (((c)#3,(((d)#1)#2,a)),(#2,((#3,b),#1))) cost=3
   =: (((c)#3,(((d)#1)#2,a)),(#3,(#2,(b,#1)))) cost=3
   =: (((c)#3,(((d)#1)#2,a)),(#2,(#3,(b,#1)))) cost=3
   >: (((c)#3,((((#3,d))#1)#2,a)),(#2,(b,#1))) cost=1
   =: (((c)#3,(((#3)#1)#2,a)),(#2,(b,(#1,d)))) cost=1
   =: (((c)#3,(((#3)#1)#2,a)),((#2,d),(b,#1))) cost=1
   =: (((c)#3,((((#3)#1,d))#2,a)),(#2,(b,#1))) cost=1
   =: (((c)#3,((((#3,d))#1)#2,a)),((b,#2),#1)) cost=1
   =: ((((c)#3,((((#3,d))#1)#2,a)),#1),(#2,b)) cost=1
   =: (((c)#3,((((#3,d))#1)#2,a)),((#2,b),#1)) cost=1
   =: (((c)#3,((((#3,d))#1)#2,(a,#1))),(#2,b)) cost=1
   =: (((c)#3,(((((#3,d))#1)#2,a),#1)),(#2,b)) cost=1
   =: ((((c)#3,((((#3,d))#1)#2,a)),#2),(b,#1)) cost=1
   >: (#2,(((c)#3,((((#3,d))#1)#2,a)),(b,#1))) cost=0
   =: (#2,(((c)#3,(((b,((#3,d))#1))#2,a)),#1)) cost=0
   =: (#2,(((c)#3,a),(((((#3,d))#1)#2,b),#1))) cost=0
   =: (#2,(((((#3,d))#1)#2,((c)#3,a)),(b,#1))) cost=0
   =: (#2,(((c)#3,a),((((#3,d))#1)#2,(b,#1)))) cost=0
   =: (#2,((((#3,d))#1)#2,(((c)#3,a),(b,#1)))) cost=0
   =: (#2,((((c)#3,a),(((#3,d))#1)#2),(b,#1))) cost=0
   =: (#2,(((((#3,d))#1)#2,((c)#3,a)),(b,#1))) cost=0
Cost:0 Steps:98 Climbs:4 TopNetworks:6 Class:Relaxed TimeConsistency:0 HCruns:1 HCTime:0.000496626 Naive:99 Naivetime:0.000190973 DTcnt:769
```

```
> cat odt.log
(#2,(((c)#3,((((#3,d))#1)#2,a)),(b,#1)))
(#2,(((c)#3,(((b,((#3,d))#1))#2,a)),#1))
(#2,(((c)#3,a),(((((#3,d))#1)#2,b),#1)))
(#2,(((((#3,d))#1)#2,((c)#3,a)),(b,#1)))
(#2,(((c)#3,a),((((#3,d))#1)#2,(b,#1))))
(#2,((((#3,d))#1)#2,(((c)#3,a),(b,#1))))
```

#### General networks run (`--general`)
```
> supnet -g "(a,(b,(c,d))); ((a,b),(c,d)); ((a,c),(b,d))" -r1 -R3 --HC -v2 -z 12 --general
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
   =: ((#2,((#3)#2,a)),(((d)#1,b),((c)#3,#1))) cost=1
   =: ((#2,((#3)#2,a)),(((d)#1,((c)#3,#1)),b)) cost=1
   =: (((#2,((#3)#2,a)),(d)#1),(b,((c)#3,#1))) cost=1
   =: ((d)#1,((#2,((#3)#2,a)),(b,((c)#3,#1)))) cost=1
   =: ((#2,((d)#1,(b,((c)#3,#1)))),((#3)#2,a)) cost=1
   =: (#2,(((#3)#2,a),((d)#1,(b,((c)#3,#1))))) cost=1
   =: ((((d)#1,#2),((#3)#2,a)),(b,((c)#3,#1))) cost=1
   =: ((#2,((d)#1,((#3)#2,a))),(b,((c)#3,#1))) cost=1
   =: ((#2,((#3)#2,a)),(((d)#1,b),((c)#3,#1))) cost=1
   =: ((#2,((#3)#2,a)),(b,((d)#1,((c)#3,#1)))) cost=1
   =: (((d)#1,(#2,((#3)#2,a))),(b,((c)#3,#1))) cost=1
   =: ((#2,(((d)#1,(#3)#2),a)),(b,((c)#3,#1))) cost=1
   =: ((#2,((#3)#2,((d)#1,a))),(b,((c)#3,#1))) cost=1
   =: ((#2,((#3)#2,a)),(b,(((d)#1,(c)#3),#1))) cost=1
   =: ((#2,((((d)#1,#3))#2,a)),(b,((c)#3,#1))) cost=1
   =: ((#2,((#3)#2,a)),(b,((((d)#1,c))#3,#1))) cost=1
   =: (((#2,(b,((c)#3,#1))),((#3)#2,a)),(d)#1) cost=1
   =: ((#2,(((#3)#2,a),(b,((c)#3,#1)))),(d)#1) cost=1
   =: (((#2,((#3)#2,a)),(b,((c)#3,#1))),(d)#1) cost=1
   =: ((#2,a),((d)#1,(((#3)#2,b),((c)#3,#1)))) cost=1
   =: ((#2,a),((d)#1,(b,((#3)#2,((c)#3,#1))))) cost=1
   =: ((#2,a),((#3)#2,((d)#1,(b,((c)#3,#1))))) cost=1
   =: ((#2,a),(((#3)#2,(d)#1),(b,((c)#3,#1)))) cost=1
   =: ((#2,a),((d)#1,((#3)#2,(b,((c)#3,#1))))) cost=1
   =: ((#2,a),((d)#1,(b,(((#3)#2,(c)#3),#1)))) cost=1
   =: ((#2,a),((d)#1,(b,((c)#3,((#3)#2,#1))))) cost=1
   =: ((#2,a),((((#3)#2,d))#1,(b,((c)#3,#1)))) cost=1
   =: ((#2,((#3)#2,a)),(((c)#3,(d)#1),(b,#1))) cost=1
   =: ((#2,((#3)#2,a)),((((c)#3,d))#1,(b,#1))) cost=1
   =: ((#2,((#3)#2,a)),((d)#1,(b,((c,#1))#3))) cost=1
Cost:1 Steps:92 Climbs:1 TopNetworks:26 Class:General TimeConsistency:0 HCruns:1 HCTime:0.000586033 Naive:93 Naivetime:0.000185251 DTcnt:744
```

```
> cat odt.log
((#2,((#3)#2,a)),((d)#1,(b,((c)#3,#1))))
(((#3)#2,a),((d)#1,((#2,b),((c)#3,#1))))
(((#3)#2,a),((d)#1,(b,(#2,((c)#3,#1)))))
(((#3)#2,a),(#2,((d)#1,(b,((c)#3,#1)))))
(((#3)#2,a),((#2,(d)#1),(b,((c)#3,#1))))
(((#3)#2,a),((d)#1,(#2,(b,((c)#3,#1)))))
(((#3)#2,a),((d)#1,(b,((#2,(c)#3),#1))))
(((#3)#2,a),((d)#1,(b,((c)#3,(#2,#1)))))
(((#3)#2,a),(((#2,d))#1,(b,((c)#3,#1))))
(#2,(((d)#1,(b,((c)#3,#1))),((#3)#2,a)))
(#2,(((d)#1,((#3)#2,a)),(b,((c)#3,#1))))
(#2,((d)#1,((b,((c)#3,#1)),((#3)#2,a))))
((#2,((#3)#2,a)),(b,((d)#1,((c)#3,#1))))
((#2,((#3)#2,a)),((b,(d)#1),((c)#3,#1)))
(((#2,((#3)#2,a)),(d)#1),(b,((c)#3,#1)))
((d)#1,((#2,((#3)#2,a)),(b,((c)#3,#1))))
((((d)#1,#2),((#3)#2,a)),(b,((c)#3,#1)))
((#2,((d)#1,((#3)#2,a))),(b,((c)#3,#1)))
((#2,(((d)#1,(#3)#2),a)),(b,((c)#3,#1)))
((#2,((#3)#2,((d)#1,a))),(b,((c)#3,#1)))
((#2,((((d)#1,#3))#2,a)),(b,((c)#3,#1)))
((#2,((#3)#2,a)),(b,((((d)#1,c))#3,#1)))
(((#2,(b,((c)#3,#1))),((#3)#2,a)),(d)#1)
((#2,(((#3)#2,a),(b,((c)#3,#1)))),(d)#1)
((#2,((#3)#2,a)),((((c)#3,d))#1,(b,#1)))
((#2,((#3)#2,a)),((d)#1,(b,((c,#1))#3)))
```

Note that using a broader network class does not guarantee a better cost.


#### General and time consistent networks `--general --timeconsistent`

```
> supnet -g "(a,(b,(c,d))); ((a,b),(c,d)); ((a,c),(b,d))" -r1 -R3 --HC -v2 -z 12 --general --timeconsistent
   i: ((c,((#2)#1,((#3)#2,((d)#3,a)))),(b,#1)) cost=4
   >: (((#2)#1,((#3)#2,((d)#3,a))),(b,(c,#1))) cost=3
   >: (((c,(#2)#1),((#3)#2,((d)#3,a))),(b,#1)) cost=2
   =: (((#2)#1,((c,(#3)#2),((d)#3,a))),(b,#1)) cost=2
   =: (((#2)#1,((#3)#2,((c,(d)#3),a))),(b,#1)) cost=2
   =: (((c,(#2)#1),(b,#1)),((#3)#2,((d)#3,a))) cost=2
   =: ((c,(#2)#1),(((#3)#2,((d)#3,a)),(b,#1))) cost=2
   =: (((#3)#2,((d)#3,a)),((c,(#2)#1),(b,#1))) cost=2
   =: ((c,(#2)#1),((b,#1),((#3)#2,((d)#3,a)))) cost=2
   =: (((((#3)#2,c),(#2)#1),((d)#3,a)),(b,#1)) cost=2
   =: (((c,(#2)#1),((d)#3,a)),(((#3)#2,b),#1)) cost=2
   =: (((((d)#3,c),(#2)#1),((#3)#2,a)),(b,#1)) cost=2
   =: (((c,(#2)#1),((#3)#2,a)),(((d)#3,b),#1)) cost=2
Cost:2 Steps:33 Climbs:3 TopNetworks:9 Class:General TimeConsistency:1 HCruns:1 HCTime:0.000252724 Naive:34 Naivetime:6.69956e-05 DTcnt:272
```


### Detecting network class 

TODO

### odt.dat in labelled format `--odtlabelled`

```
> supnet -g "(a,(b,(c,d))); ((a,b),(c,d))" -r1 -R2 -z100 --HC --odtlabelled -v0
Cost:0 Steps:62 Climbs:3 TopNetworks:5 Class:TreeChild TimeConsistency:0 HCruns:1 HCTime:0.000192642 Naive:63 Naivetime:5.65052e-05 DTcnt:244
```

```
> cat odt.dat
optcost=0
time=0.000192881
hctime=0.000192642
mergetime=2.38419e-07
topnets=5
class=0
timeconsistency=0
improvements=3
steps=62
bbruns=0
startnets=1
memoryMB=6
dtcnt=244
naivetime=5.65052e-05
naivecnt=63
bbnaivecnt=0
bbnaivetime=0
bbdpcnt=0
bbdptime=0
```

### More examples

Larger instance; tree-child search:
```  
> supnet -g "((i,c),((a,d),(f,(b,((g,(e,j)),h)))));((i,h),((c,f),((a,(d,e)),(j,(g,b)))));(((f,(b,d)),((j,g),(e,i))),(c,(h,a)));(((f,(i,(j,d))),(c,b)),(((h,a),g),e));((((h,e),((c,f),a)),(d,b)),((g,i),j))" -R8 -q1 --HC --hcrunstats
   i: ((((((((a)#1,((((#1,((((((g)#5,d))#8,((f)#6,j)))#3,h)),#5))#2,((((e,#2))#7,b),#8))))#4,(#6,c)),i),#4),#3),#7) cost=36
   >: ((((((a)#1,((((#1,((((((g)#5,d))#8,((f)#6,j)))#3,h)),#5))#2,((((e,#2))#7,b),#8))))#4,(#6,c)),((i,#4),#3)),#7) cost=34
   >: (((i,#4),#3),(((((a)#1,((((#1,((((((g)#5,d))#8,((f)#6,j)))#3,h)),#5))#2,((((e,#2))#7,b),#8))))#4,(#6,c)),#7)) cost=33
   >: ((((i,#7),#4),#3),((((a)#1,((((#1,((((((g)#5,d))#8,((f)#6,j)))#3,h)),#5))#2,((((e,#2))#7,b),#8))))#4,(#6,c))) cost=32
   >: ((((i,#4),#3),#7),((((a)#1,((((#1,((((((g)#5,d))#8,((f)#6,j)))#3,h)),#5))#2,((((e,#2))#7,b),#8))))#4,(#6,c))) cost=31
   >: (((i,#3),#7),((((a)#1,((((#1,((((((g)#5,d))#8,((f)#6,j)))#3,h)),#5))#2,((((e,#2))#7,b),#8))))#4,(#6,(c,#4)))) cost=30
   >: (((i,#3),#7),((((a)#1,(((#1,((((((g)#5,d))#8,((f)#6,j)))#3,h)))#2,(((((e,#5),#2))#7,b),#8))))#4,(#6,(c,#4)))) cost=28
   >: (((i,#3),#7),((((a)#1,(((#1,((((((g)#5,d))#8,((f)#6,j)))#3,h)))#2,((((#5,(e,#2)))#7,b),#8))))#4,(#6,(c,#4)))) cost=27
   >: (((i,#3),#7),(((((#1,((((((g)#5,d))#8,((f)#6,j)))#3,h)))#2,((a)#1,((((#5,(e,#2)))#7,b),#8))))#4,(#6,(c,#4)))) cost=26
   >: (((i,#3),#7),(((((#1,h))#2,((a)#1,((((#5,(((((((g)#5,d))#8,((f)#6,j)))#3,e),#2)))#7,b),#8))))#4,(#6,(c,#4)))) cost=25
   >: (((((#1,h))#2,(i,#3)),#7),((((a)#1,((((#5,(((((((g)#5,d))#8,((f)#6,j)))#3,e),#2)))#7,b),#8)))#4,(#6,(c,#4)))) cost=24
   >: ((((#1,h))#2,((i,#3),#7)),((((a)#1,((((#5,(((((((g)#5,d))#8,((f)#6,j)))#3,e),#2)))#7,b),#8)))#4,(#6,(c,#4)))) cost=23
   >: ((((#1,h))#2,(((g)#5,(i,#3)),#7)),((((a)#1,((((#5,(((((d)#8,((f)#6,j)))#3,e),#2)))#7,b),#8)))#4,(#6,(c,#4)))) cost=22
   >: ((((#1,h))#2,(((i,(g)#5),#3),#7)),((((a)#1,((((#5,(((((d)#8,((f)#6,j)))#3,e),#2)))#7,b),#8)))#4,(#6,(c,#4)))) cost=21
   >: ((((#1,h))#2,(((i,(g)#5),#3),#7)),((((a)#1,(((((((((d)#8,((f)#6,j)))#3,e),#5),#2))#7,b),#8)))#4,(#6,(c,#4)))) cost=20
   >: ((((#1,h))#2,(((f)#6,((i,(g)#5),#3)),#7)),((((a)#1,(((((((((d)#8,j))#3,e),#5),#2))#7,b),#8)))#4,(#6,(c,#4)))) cost=19
   >: ((((#1,h))#2,((((f,#4))#6,((i,(g)#5),#3)),#7)),((((a)#1,(((((((((d)#8,j))#3,e),#5),#2))#7,b),#8)))#4,(#6,c))) cost=17
1. Cost:17 Steps:2281 Climbs:17 TopNetworks:8 Class:TreeChild TimeConsistency:0
 New optimal cost: 17
Cost:17 Steps:2281 Climbs:17 TopNetworks:8 Class:TreeChild TimeConsistency:0 HCruns:1 HCTime:0.744191 Naive:2282 Naivetime:0.738144 DTcnt:584192
```

Recommended with large HC-runs using quasi-consensus rand networks
```  
> supnet -g "(a,(b,(c,d))); ((a,b),(c,d))" -q1000 -R3 --HC -v1 | tail -20
   i: (((d)#2,((((((b)#1,a),#1))#3,c),#2)),#3) cost=1
   >: ((((d)#2,(((((b)#1,a))#3,c),#2)),#1),#3) cost=0
   i: (((((a,#2))#3,d))#1,((b)#2,((c,#1),#3))) cost=1
   >: (((((a,#2))#3,d))#1,(((c,#1),(b)#2),#3)) cost=0
   i: (((d)#2,(#3,((b)#3,(((c,#2))#1,a)))),#1) cost=1
   >: (((d)#2,((b)#3,(((#3,(c,#2)))#1,a))),#1) cost=0
   i: (((((b)#3,a),#3))#1,((((c)#2,d),#2),#1)) cost=1
   >: ((((b)#3,a))#1,(((((c)#2,d),#3),#2),#1)) cost=0
   i: ((((((c)#1,a))#2,(d,#1)))#3,((b,#2),#3)) cost=1
   >: ((((((c)#1,a))#2,(d,#1)))#3,(#2,(b,#3))) cost=0
   i: ((((#3,((c)#1,((b)#3,a))))#2,(#1,d)),#2) cost=1
   >: ((#3,((((c)#1,((b)#3,a)))#2,(#1,d))),#2) cost=0
   i: (#2,(#1,(((#3,((c)#1,((b)#3,a))))#2,d))) cost=1
   >: (#2,(#3,(#1,((((c)#1,((b)#3,a)))#2,d)))) cost=0
   i: (((b)#3,((((#1,c))#2,((d)#1,a)),#3)),#2) cost=1
   >: (((b)#3,((((#1,c),#3))#2,((d)#1,a))),#2) cost=0
   i: (((c)#1,(((a)#2,(((#2,b))#3,d)),#1)),#3) cost=1
   >: (((c)#1,((a)#2,(((#2,b))#3,(d,#1)))),#3) cost=0
   i: ((((c)#2,a))#1,(#3,(#2,(((b,#1))#3,d)))) cost=0
Cost:0 Steps:36650 Climbs:1626 TopNetworks:254 Class:TreeChild TimeConsistency:0 HCruns:1000 HCTime:0.170419 Naive:37650 Naivetime:0.0531425 DTcnt:242997
```

Using stopping criterion (--hcstopinit=1000) with quasi-consensus rand networks.

```  
> supnet -g "(a,(b,(c,d))); ((a,b),(c,d))" -q-1 -R3 --HC --hcrunstats --hcstopinit=1000 | tail -10
   >: (((#1,c))#3,(#2,(((((a)#2,d))#1,b),#3))) cost=0
   i: (((((#3,((d)#3,c)))#1,a))#2,((b,#2),#1)) cost=0
   i: ((((a)#2,(#2,b)))#1,(#1,(#3,((c)#3,d)))) cost=1
   >: ((((a)#2,b))#1,(#2,(#1,(#3,((c)#3,d))))) cost=0
   i: ((((d)#1,(((((b)#2,a),#2))#3,c)),#1),#3) cost=1
   >: (((((d)#1,((((b)#2,a))#3,c)),#2),#1),#3) cost=0
   i: ((((((((a)#1,d))#2,b),#1))#3,(c,#2)),#3) cost=1
   >: (((a)#1,(((((d)#2,b),#1))#3,(c,#2))),#3) cost=0
   i: (((((((a)#1,(#1,b)))#2,d))#3,(c,#3)),#2) cost=1
Cost:0 Steps:48504 Climbs:2131 TopNetworks:254 Class:TreeChild TimeConsistency:0 HCruns:1297 HCTime:0.222054 Naive:49801 Naivetime:0.0698125 DTcnt:318725
```

## More examples

Print min total cost 10 random gene trees vs random tree-child network with 5 reticulations over 8 species; print the initial network.

```
> supnet -r10 -A8 --pnetworks  | supnet -G- -r1 -A8 -R5 --pnetworks --odtnaivecost
(((((h)#4,(#2,f)))#1,((((e)#3,((c)#5,g)))#2,((a,((#5,d),#3)),#4))),(#1,b))
67
```


