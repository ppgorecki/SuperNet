# SuperNet

Inferring phylogenetic networks from from collections of gene trees using Deep Coalescence.

## Usage

```
supnet ...
```

## Reading and printing input files

### Gene tree(s)

Gene trees from a string `-g` with a semicolon as separator. Print gene trees via `-eg`.

```
> supnet -g '(a,(b,c));(a,(a,d))' -eg
(a,(b,c))
(a,(a,d))
```

Gene trees from a file `-G`. 

```
> supnet -G examples/gtrees.txt -eg
(a,(b,c))
(a,(a,d))
```

Use `-` for stdin.

```
> cat examples/gtrees.txt | supnet -G- -eg
(a,(b,c))
(a,(a,d))
```

### Species tree(s) 

Similarly, use `-s/-S` and `-es` to read and print species tree sets.

### Networks

Use `-n/-N` and `-en` to read and print networks.

### (Not implemented yet) Matching gene tree names with species names

- `-l [aDELIM|bDELIM|pPOS]`: rules for matching gene tree names with species tree names;
    - `a`: after delimiter
    - `b`: before delimiter
    - `p`: from position

## Random trees and networks and quasi consensus generators

- `-r NUM`: generate NUM random species trees/networks
- `-q NUM`: generate NUM quasi-consensus species trees/networks using a given set of gene trees
- `-R NUM`: inserts NUM reticulation nodes in all networks/trees from -q, -r, -n, -N (default 0)
- `-A SPECIESNUM`: define SPECIESNUM species a,b,...,z

### Random trees
    
Print 3 random species trees over a..e (5)
```
>  supnet -A5 -r3 -en
(a,(((e,c),d),b))
(e,((d,a),(c,b)))
((((a,d),e),b),c)
```

Print two quasi consensus trees.
```
> supnet -g '(a,((b,c),d));(a,(b,d))' -q2 -en
(d,(a,(c,b)))
(a,((b,c),d))
```

### Random networks 

By default, the generator produces tree-child networks. For other clasess use:
- `-e1`: generate class 1 networks, i.e., each internal node has at most one reticulation child
- `-e2`: general phylogenetic networks (no restrictions)

See also `-e` and `-o` option.

Print one random tree-child network trees over a..e
```
>  supnet -A5 -r1 -en
((a,b),((d,e),c))
```

Print one random class 1 network trees over abc with 10 reticulations.
```
>  supnet -A3 -R5 -r1 -en1
((#5)#3,(#4,((((b,#1))#2)#5,(#3,((a)#1,((#2)#4,c))))))
```

Print one quasi consensus network with two reticulations.
```
> supnet -G examples/gtrees.txt -q1 -R2 -en
((((((d)#1,c),b))#2,(a,#1)),#2)
```

### Fixed seed for random generator `-z`

Use `-z SEED` to set seed for random generator (srand).


### Inserting reticulations

Insert 2 reticulations into a network; tree-child network in output (default)
```  
> supnet -R2 -n '(a,((d)#1,(b,(c,#1))))' -en
(#3,((a)#3,((((b)#2,d))#1,(#2,(c,#1)))))
```

Insert 10 reticulations into a network; general network `-e2`
```  
> supnet -R10 -n '(a,((d)#1,(b,(c,#1))))' -en2
(#11,(#7,(((((#8)#4,a))#7,(d)#3),((((c)#10,((((#10)#6)#11,#3))#9),(#2)#1),(#4,((((((b)#8,(#9)#2),#1))#5,#6),#5))))))
```

## Display trees of networks

Print display trees based on reticulation switching `-et`. Trees maybe non-unique.
```
>  supnet -n '((((c)#B,b))#A,(#A,(#B,a)))' -et
((c,a),b)
((c,b),a)
((c,a),b)
((c,b),a)
```

Print all display trees with their ids '-eT'.
```
>  supnet -n '((((c)#B,b))#A,(#A,(#B,a)))' -eT
0 ((c,a),b)
1 ((c,b),a)
2 ((c,a),b)
3 ((c,b),a)
```

## Cost functions `-C COST`, `-ec`, `-eC`

Use `-C COST` to set cost function from {DL,D,L,DC,RF,DCE}; the default is DC. 
Print the cost with `-ec` (only the cost) or `-eC` (trees + their cost).

Print DC cost
```
>  supnet -g '(b,(a,c))' -s '(a,(b,c))' -CDC -ec 
1
```

## Verbose option `-v [123][34]`:

   - `0` - quiet mode in HC, BB and ODT (see also `-o[sS]`)
   - `1` - in HC print visited network after each improvement (strictly)
   - `2` - in HC print visited networks if the cost is equal to the current or improved
   - `3` - in HC print all visited networks
   - `4` - in HC,BB,ODT print basic info
   - `5` - in HC,BB,ODT print more detailed info

## Optimal Display Tree (ODT) heuristics `-o`

ODT is solved by a hill climbing heuristic algorithm that optimizes a cost function for networks. The algorithm will start from a set of predefined networks or generate random networks partially based on input gene trees if specified. The algorithm will make use of two moves: "tail" and "nni" to improve the networks.

The algorithm will keep track of various statistics. It will print the optimal cost found and write all optimal networks to an "odt.log" file. Summary statistics, including the optimal cost, total time taken, time taken for hill climbing, time taken for the merge step, number of networks, number of improvements made, number of steps taken, and the number of starting networks will be saved to an "odt.dat" file.

Sketch of the algorithm:

```
Input: a set of gene trees GSet and a cost function (now only DC).

For each predefined network N:
  Run hill climbing HC(GSet, N) and store the result if the resuling cost is not worst than the current one.
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
    if `-oe` is set or N has less than '-t THRESHOLD` reticulations:
       compute the cost by naive algorithm (enumeration of all possible display trees)
    otherwise:        
       compute the cost by BB algorithm (using multiple calls of DP)
```

By default each HC step starts from the set of predefined networks (see `-n`, `-N` or random networks); otherwise by using quasi-random networks.

Main option `-o` with suboptions `[TNt123sq]+`:
   - `T` - use TailMoves (default)
   - `N` - use NNI instead of TailMoves
   - `s` - print summary when a single HC run is completed
   - `S` - as above but only when better networks is found
   - `q` - do not save odt.log/odt.dat 
   - `e` - use a naive approach to compute ODT cost 

By default searches and initial networks are limited to tree-child. 
Use `-e1` or `-e2` to extend the search to relaxed or general networks, resp.

Additional options:
  - `-t THRESHOLD` - run naive odt computation when reticulation number is `< THRESHOLD`; default is 13
  - `-K` NUM - stopping criterion: stop when there is no new network after NUM HC runs
  - `-E float`: set naive sampling level (exponential distribution)    
  - `-O` ODTFILE - the name of odt.log 
  - `-D` DATFILE - the name of odt.dat; see also -ea 

#### Minimalistic run with a single network

Minimalistic hill climbing (HC) run using Tail Moves (default) with a single random network (and one HC run) `-r1` with two reticulations `-R2`. Print cost and save result to odt.log. 

```
> supnet -g "(a,(b,(c,d))); ((a,b),(c,d))" -r1 -R2 -z100
```

#### Stats in odt with `-oS` and `-os` used with multiple starting networks.

Print stats when HC is completed with improvement `-oS`. Do not print detailed stats for individual HC runs `-v0`.
```
> supnet -g "(a,(b,(c,d))); ((a,b),(c,d)); ((a,c),(b,d)); ((a,d),e)" -r10 -R1 -oS -v0 -z 1003
1. Cost:7 Steps:128 Climbs:6 TopNetworks:6 Class:TreeChild
NewOptCost: 7
2. Cost:6 Steps:108 Climbs:4 TopNetworks:2 Class:TreeChild
NewOptCost: 6
3. Cost:5 Steps:146 Climbs:4 TopNetworks:11 Class:TreeChild
NewOptCost: 5
4. Cost:4 Steps:56 Climbs:3 TopNetworks:4 Class:TreeChild
NewOptCost: 4
5. Cost:3 Steps:198 Climbs:7 TopNetworks:1 Class:TreeChild
NewOptCost: 3
Cost:3 Steps:1080 Climbs:44 TopNetworks:1 Class:TreeChild HCruns:10 HCTime:0.00315738 Naive:1090 Naivetime:0.000870228 DTcnt:2180
```

Print stats after completed HC run for each initial network `-os`. 
```
> supnet -g "(a,(b,(c,d))); ((a,b),(c,d)); ((a,c),(b,d)); ((a,d),e)" -r10 -R1 -os -v0 -z 1003
1. Cost:7 Steps:128 Climbs:6 TopNetworks:6 Class:TreeChild
NewOptCost: 7
2. Cost:6 Steps:108 Climbs:4 TopNetworks:2 Class:TreeChild
NewOptCost: 6
3. Cost:5 Steps:146 Climbs:4 TopNetworks:11 Class:TreeChild
NewOptCost: 5
4. Cost:4 Steps:56 Climbs:3 TopNetworks:4 Class:TreeChild
NewOptCost: 4
5. Cost:3 Steps:198 Climbs:7 TopNetworks:1 Class:TreeChild
NewOptCost: 3
6. Cost:5 Steps:80 Climbs:4 TopNetworks:11 Class:TreeChild
7. Cost:4 Steps:70 Climbs:3 TopNetworks:6 Class:TreeChild
8. Cost:9 Steps:52 Climbs:2 TopNetworks:12 Class:TreeChild
10. Cost:4 Steps:105 Climbs:6 TopNetworks:6 Class:TreeChild
Cost:3 Steps:1080 Climbs:44 TopNetworks:1 Class:TreeChild HCruns:10 HCTime:0.00337648 Naive:1090 Naivetime:0.000950575 DTcnt:2180
```

#### HC, BB and DP algorithms verbose `-v`

Print only improvements in HC `-v1`:
```  
> supnet -g "(a,(b,(c,d))); ((a,b),(c,d)); ((a,c),(b,d))" -r1 -R1 -o_ -v1 -z 1001
   i: ((b)#1,(c,((d,#1),a))) cost=5
   >: ((b)#1,((d,c),(#1,a))) cost=3
   >: (((b)#1,a),((d,c),#1)) cost=2
Cost:2 Steps:58 Climbs:3 TopNetworks:3 Class:TreeChild HCruns:1 HCTime:0.000189781 Naive:59 Naivetime:4.00543e-05 DTcnt:118
```

Print improvements and equal cost networks `-v2`:
```  
> supnet -g "(a,(b,(c,d))); ((a,b),(c,d)); ((a,c),(b,d))" -r1 -R1 -o_ -v2 -z 1001
   i: ((b)#1,(c,((d,#1),a))) cost=5
   >: ((b)#1,((d,c),(#1,a))) cost=3
   =: ((d,(b)#1),(c,(#1,a))) cost=3
   =: ((b)#1,(((#1,d),c),a)) cost=3
   =: ((b)#1,((#1,(d,c)),a)) cost=3
   >: (((b)#1,a),((d,c),#1)) cost=2
   =: (((b)#1,a),(c,(d,#1))) cost=2
   =: (((b)#1,a),((d,#1),c)) cost=2
   =: (a,((((b)#1,d),c),#1)) cost=2
   =: (((b)#1,a),((d,#1),c)) cost=2
Cost:2 Steps:58 Climbs:3 TopNetworks:3 Class:TreeChild HCruns:1 HCTime:0.000197411 Naive:59 Naivetime:4.1008e-05 DTcnt:118
```

Print every visited network `-v3`:
```  
> supnet -g "(a,(b,(c,d))); ((a,b),(c,d)); ((a,c),(b,d))" -r1 -R1 -o_ -v3 -z 1001
   i: ((b)#1,(c,((d,#1),a))) cost=5
   <: ((b)#1,(c,(#1,(d,a)))) cost=6
   >: ((b)#1,((d,c),(#1,a))) cost=3
   <: ((b)#1,(c,((d,#1),a))) cost=5
   <: ((b)#1,(c,(#1,(d,a)))) cost=6
   <: ((b)#1,(c,(d,(#1,a)))) cost=5
   =: ((d,(b)#1),(c,(#1,a))) cost=3
   <: ((b)#1,(d,(c,(#1,a)))) cost=4
   <: (((d,b))#1,(c,(#1,a))) cost=4
   <: ((b)#1,(d,((#1,c),a))) cost=5
   <: ((b)#1,(d,(#1,(a,c)))) cost=5
   <: ((b)#1,(d,((#1,a),c))) cost=4
   <: (((b)#1,c),(d,(#1,a))) cost=5
   <: ((b)#1,((d,(#1,a)),c)) cost=5
   <: (((b,c))#1,(d,(#1,a))) cost=6
   =: ((b)#1,(((#1,d),c),a)) cost=3
   <: ((b)#1,((d,(#1,c)),a)) cost=4
   =: ((b)#1,((#1,(d,c)),a)) cost=3
   <: ((b)#1,(#1,((d,c),a))) cost=4
   <: ((b)#1,(((d,a),c),#1)) cost=6
   <: ((b)#1,((d,(c,a)),#1)) cost=5
   <: ((b)#1,(((d,c),a),#1)) cost=4
   >: (((b)#1,a),((d,c),#1)) cost=2
   <: (((d,(b)#1),a),(c,#1)) cost=5
   <: (((b)#1,(d,a)),(c,#1)) cost=6
   =: (((b)#1,a),(c,(d,#1))) cost=2
   <: ((d,((b)#1,a)),(c,#1)) cost=5
   <: (((b)#1,a),(d,(c,#1))) cost=3
   <: ((((d,b))#1,a),(c,#1)) cost=4
   <: ((((b)#1,c),a),(d,#1)) cost=4
   <: (((b)#1,(a,c)),(d,#1)) cost=4
   <: (((b)#1,a),(d,(#1,c))) cost=3
   <: ((((b)#1,a),c),(d,#1)) cost=3
   =: (((b)#1,a),((d,#1),c)) cost=2
   <: ((((b,c))#1,a),(d,#1)) cost=5
   =: (a,((((b)#1,d),c),#1)) cost=2
   <: (a,((d,((b)#1,c)),#1)) cost=3
   <: (a,(((b)#1,(d,c)),#1)) cost=3
   <: (a,((b)#1,((d,c),#1))) cost=3
   <: ((b)#1,(((d,a),c),#1)) cost=6
   <: ((b)#1,((d,(c,a)),#1)) cost=5
   <: ((b)#1,(((d,c),a),#1)) cost=4
   <: ((b)#1,((d,c),(#1,a))) cost=3
   <: ((b)#1,(((d,c),#1),a)) cost=3
   <: (((b,a))#1,((d,c),#1)) cost=3
   <: ((((d,c),(b)#1),a),#1) cost=3
   <: (((b)#1,((d,c),a)),#1) cost=4
   <: (((d,c),((b)#1,a)),#1) cost=3
   <: (((((d,c),b))#1,a),#1) cost=3
   =: (((b)#1,a),((d,#1),c)) cost=2
   <: (((b)#1,a),(d,(c,#1))) cost=3
   <: (((b)#1,(a,#1)),(d,c)) cost=3
   <: ((((b)#1,a),#1),(d,c)) cost=3
   <: (((((b)#1,a),d),c),#1) cost=5
   <: ((d,(((b)#1,a),c)),#1) cost=4
   <: ((((b)#1,a),(d,c)),#1) cost=3
   <: ((d,c),(((b)#1,a),#1)) cost=3
   <: (((b)#1,((d,c),#1)),a) cost=3
   <: ((b)#1,(a,((d,c),#1))) cost=3
Cost:2 Steps:58 Climbs:3 TopNetworks:3 Class:TreeChild HCruns:1 HCTime:0.000253439 Naive:59 Naivetime:4.07696e-05 DTcnt:118
```

To show details on DP & BB algorithms use `-v` with 4, 5 and 6. See examples below.


#### Run naive odt computation threshold `-t THRESHOLD`

Run naive odt computation, when the number of reticulations in a network is `< THRESHOLD`; otherwise run BB & DP algorithms.

Higher `-t NUM`, more naive algorithm computations.

```
> supnet -g "(a,(b,(c,d)))" -n "(#1,((b,((d,c))#1),a))" -t2 -o_ -v35
HC start: usenaive_oe=0 runnaiveleqrt_t=2
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
Cost:0 Steps:30 Climbs:1 TopNetworks:6 Class:TreeChild HCruns:1 HCTime:0.000154734 Naive:31 Naivetime:1.45435e-05 DTcnt:57
Optimal networks saved: odt.log
Stats data save to: odt.dat
```

Lower `-t`, more BB & DP computations
```
> supnet -g "(a,(b,(c,d)))" -n "(#1,((b,((d,c))#1),a))" -t1 -o_ -v35
HC start: usenaive_oe=0 runnaiveleqrt_t=1
ODT-cost-BB 1 tree(s); rt=1 cost=0
   i: (#1,((b,((d,c))#1),a)) cost=0
ODT-cost-BB 1 tree(s); rt=1 cost=1
   <: (#1,(((d,b),(c)#1),a)) cost=1
ODT-cost-BB 1 tree(s); rt=1 cost=0
   =: (#1,((b,(d,(c)#1)),a)) cost=0
ODT-cost-BB 1 tree(s); rt=1 cost=1
   <: (#1,((d,(b,(c)#1)),a)) cost=1
ODT-cost-BB 1 tree(s); rt=1 cost=2
   <: (#1,((b,(c)#1),(d,a))) cost=2
ODT-cost-BB 1 tree(s); rt=1 cost=1
   <: ((d,#1),((b,(c)#1),a)) cost=1
ODT-cost-BB 1 tree(s); rt=1 cost=3
   <: (#1,(d,((b,(c)#1),a))) cost=3
ODT-cost-BB 1 tree(s); rt=1 cost=1
   <: (#1,(((b,c),(d)#1),a)) cost=1
ODT-cost-BB 1 tree(s); rt=1 cost=0
   =: (#1,((b,((d)#1,c)),a)) cost=0
ODT-cost-BB 1 tree(s); rt=1 cost=1
   <: (#1,(((b,(d)#1),c),a)) cost=1
ODT-cost-BB 1 tree(s); rt=1 cost=2
   <: (#1,((b,(d)#1),(a,c))) cost=2
ODT-cost-BB 1 tree(s); rt=1 cost=1
   <: ((#1,c),((b,(d)#1),a)) cost=1
ODT-cost-BB 1 tree(s); rt=1 cost=3
   <: (#1,(((b,(d)#1),a),c)) cost=3
ODT-cost-BB 1 tree(s); rt=1 cost=1
   <: (#1,((((b,d),c))#1,a)) cost=1
ODT-cost-BB 1 tree(s); rt=1 cost=1
   <: (#1,(((d,(b,c)))#1,a)) cost=1
ODT-cost-BB 1 tree(s); rt=1 cost=1
   <: (#1,(((d,c))#1,(b,a))) cost=1
ODT-cost-BB 1 tree(s); rt=1 cost=0
   =: ((b,#1),(((d,c))#1,a)) cost=0
ODT-cost-BB 1 tree(s); rt=1 cost=1
   <: (#1,(b,(((d,c))#1,a))) cost=1
ODT-cost-BB 1 tree(s); rt=1 cost=0
   =: (#1,(((b,(d,c)))#1,a)) cost=0
ODT-cost-BB 1 tree(s); rt=1 cost=1
   <: (#1,(b,(a,((d,c))#1))) cost=1
ODT-cost-BB 1 tree(s); rt=1 cost=1
   <: (#1,((b,a),((d,c))#1)) cost=1
ODT-cost-BB 1 tree(s); rt=1 cost=0
   =: (((b,((d,c))#1),#1),a) cost=0
ODT-cost-BB 1 tree(s); rt=1 cost=2
   <: (#1,(b,(((d,a),c))#1)) cost=2
ODT-cost-BB 1 tree(s); rt=1 cost=2
   <: (#1,(b,((d,(c,a)))#1)) cost=2
ODT-cost-BB 1 tree(s); rt=1 cost=1
   <: (#1,((b,a),((d,c))#1)) cost=1
ODT-cost-BB 1 tree(s); rt=1 cost=1
   <: (#1,(b,(((d,c))#1,a))) cost=1
ODT-cost-BB 1 tree(s); rt=1 cost=0
   =: ((#1,a),(b,((d,c))#1)) cost=0
ODT-cost-BB 1 tree(s); rt=1 cost=1
   <: (#1,(b,(((d,c),a))#1)) cost=1
ODT-cost-BB 1 tree(s); rt=1 cost=0
   =: (((#1,b),((d,c))#1),a) cost=0
ODT-cost-BB 1 tree(s); rt=1 cost=0
   =: ((#1,(b,((d,c))#1)),a) cost=0
ODT-cost-BB 1 tree(s); rt=1 cost=0
   =: ((b,((d,c))#1),(#1,a)) cost=0
HC run completed: (#1,((b,((d,c))#1),a)) cost=0
Cost:0 Steps:30 Climbs:1 TopNetworks:6 Class:TreeChild HCruns:1 HCTime:0.000185013 BB:31 Naive:18 Naivetime:6.67572e-06 DTcnt:18 BBnaivetime:4.26769e-05 BBnaivecnt:18 BBdptime:5.19753e-05 BBdpcnt:31
Optimal networks saved: odt.log
Stats data save to: odt.dat
```

### NNI vs Tail Moves

ODT via Tail Moves (default, `-oT`)
```
> supnet -g "(a,(b,(c,d))); ((a,b),(c,d)); ((a,c),(b,d))" -r1 -R1 -oT -v3 -z 1003
   i: (#1,((((c,d))#1,b),a)) cost=2
   <: (#1,((c,((d)#1,b)),a)) cost=4
   <: (#1,(((d)#1,b),(c,a))) cost=4
   =: ((c,#1),(((d)#1,b),a)) cost=2
   <: (#1,(c,(((d)#1,b),a))) cost=5
   <: (#1,(((c,(d)#1),b),a)) cost=3
   <: (#1,(((d)#1,(c,b)),a)) cost=5
   <: (#1,((((c)#1,b),d),a)) cost=4
   <: (#1,(((c)#1,b),(a,d))) cost=6
   <: ((#1,d),(((c)#1,b),a)) cost=3
   <: (#1,((((c)#1,b),a),d)) cost=6
   =: (#1,((((c)#1,d),b),a)) cost=2
   <: (#1,(((c)#1,(b,d)),a)) cost=4
   <: (((((c,d))#1,b),#1),a) cost=3
   <: (#1,((((c,a),d))#1,b)) cost=5
   <: (#1,(((c,(d,a)))#1,b)) cost=6
   <: ((#1,a),(((c,d))#1,b)) cost=3
   <: (#1,((((c,d))#1,a),b)) cost=3
   <: (#1,(((c,d))#1,(b,a))) cost=3
   <: (#1,((((c,d),a))#1,b)) cost=4
   <: ((#1,(((c,d))#1,b)),a) cost=3
   <: ((((c,d))#1,b),(#1,a)) cost=3
   <: ((((c,d))#1,(#1,b)),a) cost=3
   <: (#1,(b,(((c,d))#1,a))) cost=3
   <: (#1,(((c,d))#1,(b,a))) cost=3
   <: (#1,((((c,b),d))#1,a)) cost=5
   <: (#1,(((c,(d,b)))#1,a)) cost=4
   <: (#1,(((c,d))#1,(a,b))) cost=3
   <: ((#1,b),(((c,d))#1,a)) cost=3
   <: (#1,((((c,d))#1,a),b)) cost=3
   <: (#1,((((c,d),b))#1,a)) cost=3
Cost:2 Steps:30 Climbs:1 TopNetworks:3 Class:TreeChild HCruns:1 HCTime:0.000142336 Naive:31 Naivetime:2.26498e-05 DTcnt:62
```

ODT via NNI moves `-oN`
```
> supnet -g "(a,(b,(c,d))); ((a,b),(c,d)); ((a,c),(b,d))" -r1 -R1 -oN -v3 -z 1003
   i: (#1,((((c,d))#1,b),a)) cost=2
   <: ((((c,d))#1,b),(a,#1)) cost=3
   <: (a,(#1,(((c,d))#1,b))) cost=3
   <: (#1,((b,a),((c,d))#1)) cost=3
   <: (#1,((a,((c,d))#1),b)) cost=3
Cost:2 Steps:4 Climbs:1 TopNetworks:1 Class:TreeChild HCruns:1 HCTime:4.31538e-05 Naive:5 Naivetime:5.72205e-06 DTcnt:10
```

### Network searches: `-e0` (default), `-e1`, `-e2`

Supnet has tree types of network types:
- `-e0` default, tree-child networks, where a non-leaf node has at least one non-reticulation child
- `-e1` relaxed networks, where a non-leaf node has at most one non-reticulation child
- `-e2` general networks

These values are reported in `[Cc]lass` field in output and odt.dat file.

The last one cannot be executed with DP algorithm.

Default `-e0`: tree-child networks and Tail Moves:
```
> supnet -g "(a,(b,(c,d))); ((a,b),(c,d)); ((a,c),(b,d))" -r1 -R3 -o_ -v2 -z 12
   i: (#3,((a)#2,(((c,#2))#3,((b)#1,(#1,d))))) cost=3
   =: (#3,((a)#2,((c)#3,(((b,#2))#1,(#1,d))))) cost=3
   =: (#3,((a)#2,((((#1,c),#2))#3,((b)#1,d)))) cost=3
   >: (#3,(#1,((a)#2,(((c,#2))#3,((b)#1,d))))) cost=2
   =: (#3,((#1,((a)#2,((c)#3,((b)#1,d)))),#2)) cost=2
   =: (#3,(#1,((a)#2,((c)#3,(((b,#2))#1,d))))) cost=2
   >: (#3,((a)#2,(#1,(((c,#2))#3,((b)#1,d))))) cost=1
   =: (#3,((a)#2,(#1,((c)#3,(((b,#2))#1,d))))) cost=1
   =: (#3,(((#1,a))#2,(((c,#2))#3,((b)#1,d)))) cost=1
   =: ((a)#2,(#1,(#3,(((c,#2))#3,((b)#1,d))))) cost=1
   =: ((a)#2,(#1,(((c,#2))#3,((b)#1,(#3,d))))) cost=1
   =: ((a)#2,(#3,(#1,(((c,#2))#3,((b)#1,d))))) cost=1
   =: ((a)#2,(#1,(((c,#2))#3,(#3,((b)#1,d))))) cost=1
   =: ((a)#2,(#1,(((c,#2))#3,(((#3,b))#1,d)))) cost=1
   =: (#3,((a)#2,(#1,(((((b)#1,c),#2))#3,d)))) cost=1
   =: (#3,((a)#2,(#1,((b)#1,(((c,#2))#3,d))))) cost=1
   =: (#3,((b)#1,((a)#2,(#1,(((c,#2))#3,d))))) cost=1
   =: (#3,((a)#2,((b)#1,(#1,(((c,#2))#3,d))))) cost=1
   >: (#3,((((b)#1,a))#2,(#1,(((c,#2))#3,d)))) cost=0
   =: ((((b)#1,a))#2,(#3,(#1,(((c,#2))#3,d)))) cost=0
   =: ((#3,(#1,(((c,#2))#3,d))),(((b)#1,a))#2) cost=0
Cost:0 Steps:68 Climbs:4 TopNetworks:2 Class:TreeChild HCruns:1 HCTime:0.000372887 Naive:69 Naivetime:0.00014472 DTcnt:550
```

```
> cat odt.log
(#3,((((b)#1,a))#2,(#1,(((c,#2))#3,d))))
((((b)#1,a))#2,(#3,(#1,(((c,#2))#3,d))))
```

Class 1 (relaxed) networks. 
```
> supnet -g "(a,(b,(c,d))); ((a,b),(c,d)); ((a,c),(b,d))" -r1 -R3 -o_ -v2 -z 12 -e1
   i: ((c,a),((#3)#2,(((b)#1)#3,(#2,(#1,d))))) cost=4
   >: (a,((#3)#2,(((b)#1)#3,(#2,((c,#1),d))))) cost=3
   =: (a,((#3)#2,(((b)#1)#3,(#2,(#1,(c,d)))))) cost=3
   >: (a,((#3)#2,(((b)#1)#3,((c,#2),(#1,d))))) cost=2
   =: (a,((#3)#2,(((b)#1)#3,(#2,(c,(#1,d)))))) cost=2
   =: (a,((#3)#2,(c,(((b)#1)#3,(#2,(#1,d)))))) cost=2
   =: (a,((#3)#2,((c,((b)#1)#3),(#2,(#1,d))))) cost=2
   =: (a,((#3)#2,(((b)#1)#3,(c,(#2,(#1,d)))))) cost=2
   =: (a,(((c,#3))#2,(((b)#1)#3,(#2,(#1,d))))) cost=2
   =: (a,((#3)#2,(((b)#1)#3,(c,(#1,(d,#2)))))) cost=2
   >: ((a,#2),((#3)#2,(((b)#1)#3,(c,(#1,d))))) cost=1
   =: (#2,((a,(#3)#2),(((b)#1)#3,(c,(#1,d))))) cost=1
   =: (#2,(((a,#3))#2,(((b)#1)#3,(c,(#1,d))))) cost=1
   =: ((a,#2),((#3)#2,((((b)#1)#3,d),(c,#1)))) cost=1
   =: (((a,#2),(#3)#2),(((b)#1)#3,(c,(#1,d)))) cost=1
   =: ((#3)#2,((a,#2),(((b)#1)#3,(c,(#1,d))))) cost=1
   =: ((a,#2),((#3)#2,(c,(((b)#1)#3,(#1,d))))) cost=1
   =: ((a,#2),((#3)#2,((c,((b)#1)#3),(#1,d)))) cost=1
   =: ((a,#2),((#3)#2,((((b)#1)#3,c),(#1,d)))) cost=1
   =: ((a,#2),((#3)#2,((((b)#1)#3,(#1,d)),c))) cost=1
   =: ((((#3)#2,a),#2),(((b)#1)#3,(c,(#1,d)))) cost=1
   =: ((a,#2),(((b)#1)#3,(c,(#1,((#3)#2,d))))) cost=1
   =: (((#3)#2,(a,#2)),(((b)#1)#3,(c,(#1,d)))) cost=1
   =: ((a,#2),(((b)#1)#3,(((#3)#2,c),(#1,d)))) cost=1
   =: ((a,#2),(((b)#1)#3,(c,((#3)#2,(#1,d))))) cost=1
   =: ((a,#2),(((b)#1)#3,((#3)#2,(c,(#1,d))))) cost=1
   =: (((a,#2),(((b)#1)#3,(c,(#1,d)))),(#3)#2) cost=1
   =: (((((b)#1)#3,a),#2),((#3)#2,(c,(#1,d)))) cost=1
   =: ((a,#2),((#3)#2,(c,(#1,(((b)#1)#3,d))))) cost=1
   =: ((((b)#1)#3,(a,#2)),((#3)#2,(c,(#1,d)))) cost=1
   =: ((a,#2),(((b)#1)#3,((#3)#2,(c,(#1,d))))) cost=1
   =: ((a,#2),((#3)#2,((((b)#1)#3,c),(#1,d)))) cost=1
   =: ((a,#2),((#3)#2,(c,(((b)#1)#3,(#1,d))))) cost=1
   =: ((a,#2),(((#3)#2,(c,(#1,d))),((b)#1)#3)) cost=1
Cost:1 Steps:103 Climbs:4 TopNetworks:14 Class:Relaxed HCruns:1 HCTime:0.000551701 Naive:104 Naivetime:0.00020957 DTcnt:832
```

```
> cat odt.log
((a,#2),((#3)#2,(((b)#1)#3,(c,(#1,d)))))
(#2,((a,(#3)#2),(((b)#1)#3,(c,(#1,d)))))
(#2,(((a,#3))#2,(((b)#1)#3,(c,(#1,d)))))
((a,#2),((#3)#2,((((b)#1)#3,d),(c,#1))))
(((a,#2),(#3)#2),(((b)#1)#3,(c,(#1,d))))
((a,#2),((#3)#2,(c,(((b)#1)#3,(#1,d)))))
((a,#2),((#3)#2,((c,((b)#1)#3),(#1,d))))
((a,#2),(((b)#1)#3,(c,(#1,((#3)#2,d)))))
((a,#2),(((b)#1)#3,(((#3)#2,c),(#1,d))))
((a,#2),(((b)#1)#3,(c,((#3)#2,(#1,d)))))
((a,#2),(((b)#1)#3,((#3)#2,(c,(#1,d)))))
(((((b)#1)#3,a),#2),((#3)#2,(c,(#1,d))))
((a,#2),((#3)#2,(c,(#1,(((b)#1)#3,d)))))
((((b)#1)#3,(a,#2)),((#3)#2,(c,(#1,d))))
```

General networks run (`-e2`)
```
> supnet -g "(a,(b,(c,d))); ((a,b),(c,d)); ((a,c),(b,d))" -r1 -R3 -o_ -v2 -z 12 -e2
   i: ((((b)#3,#2),a),((#3)#1,(((c)#2,#1),d))) cost=0
   =: ((((b)#3,a),#2),((#3)#1,(((c)#2,#1),d))) cost=0
   =: (((b)#3,(#2,a)),((#3)#1,(((c)#2,#1),d))) cost=0
   =: ((((b)#3,#2),a),((#3)#1,(((c)#2,d),#1))) cost=0
   =: ((((b)#3,#2),a),((#3)#1,((c)#2,(#1,d)))) cost=0
   =: ((((b)#3,#2),a),((#3)#1,(#1,((c)#2,d)))) cost=0
   =: ((((b)#3,#2),a),((#3)#1,((c)#2,(#1,d)))) cost=0
   =: (((((b)#3,#2),#1),a),((#3)#1,((c)#2,d))) cost=0
   =: ((((b)#3,#2),(a,#1)),((#3)#1,((c)#2,d))) cost=0
   =: ((((b)#3,#2),a),((#3)#1,((c)#2,(d,#1)))) cost=0
   =: (((((b)#3,#2),a),#1),((#3)#1,((c)#2,d))) cost=0
   =: ((((b)#3,#2),a),(((#3)#1,((c)#2,d)),#1)) cost=0
   =: ((((b)#3,#2),a),((#3)#1,(((c)#2,d),#1))) cost=0
   =: ((((b)#3,#2),a),((#3)#1,(((c,#1))#2,d))) cost=0
   =: ((#2,((b)#3,a)),((#3)#1,(((c)#2,#1),d))) cost=0
   =: (((b)#3,(#2,a)),((#3)#1,(((c)#2,#1),d))) cost=0
   =: (((b)#3,(a,#2)),((#3)#1,(((c)#2,#1),d))) cost=0
   =: ((((b)#3,a),#2),((#3)#1,(((c)#2,#1),d))) cost=0
Cost:0 Steps:89 Climbs:1 TopNetworks:9 Class:General HCruns:1 HCTime:0.000433207 Naive:90 Naivetime:0.000182152 DTcnt:693
```

```
> cat odt.log
((((b)#3,#2),a),((#3)#1,(((c)#2,#1),d)))
((((b)#3,a),#2),((#3)#1,(((c)#2,#1),d)))
(((b)#3,(#2,a)),((#3)#1,(((c)#2,#1),d)))
((((b)#3,#2),a),((#3)#1,(((c)#2,d),#1)))
((((b)#3,#2),a),((#3)#1,((c)#2,(#1,d))))
(((((b)#3,#2),#1),a),((#3)#1,((c)#2,d)))
((((b)#3,#2),(a,#1)),((#3)#1,((c)#2,d)))
(((((b)#3,#2),a),#1),((#3)#1,((c)#2,d)))
((((b)#3,#2),a),((#3)#1,(((c,#1))#2,d)))
```

Note that using a broader network class does not guarantee a better cost.

### odt.dat in labelled format `-ea`

```
> supnet -g "(a,(b,(c,d))); ((a,b),(c,d))" -r1 -R2 -z100 -o_ -ea -v0
Cost:0 Steps:68 Climbs:3 TopNetworks:6 Class:TreeChild HCruns:1 HCTime:0.000213146 Naive:69 Naivetime:6.17504e-05 DTcnt:254
```

```
> cat odt.dat
optcost=0
time=0.000213146
hctime=0.000213146
mergetime=0
topnets=6
class=0
improvements=3
steps=68
bbruns=0
startnets=1
memoryMB=6
dtcnt=254
naivetime=6.17504e-05
naivecnt=69
bbnaivecnt=0
bbnaivetime=0
bbdpcnt=0
bbdptime=0
```

### More ODT examples

Print cost, improvements and stats `-ost`; tree-child search; quasi consensus initial network `-q` with 3 random reticulations `-R3`.

```  
> supnet -g "(a,(b,(c,(d,e)))); ((a,b),(c,(e,a))); ((b,c),(d,a))" -q1 -R3 -os -v1
   i: ((#2,(#3,(((a,b))#1,((((c)#2,e))#3,d)))),#1) cost=6
   >: ((#2,((((#3,a),b))#1,((((c)#2,e))#3,d))),#1) cost=5
   >: ((#2,((((#3,b))#1,((((c)#2,e))#3,d)),a)),#1) cost=4
   >: (((((#3,b))#1,(#2,((((c)#2,e))#3,d))),a),#1) cost=3
1. Cost:3 Steps:90 Climbs:4 TopNetworks:3 Class:TreeChild
NewOptCost: 3
Cost:3 Steps:90 Climbs:4 TopNetworks:3 Class:TreeChild HCruns:1 HCTime:0.000480413 Naive:91 Naivetime:0.000240803 DTcnt:728
```

Larger instance; tree-child search:
```  
> supnet -g "((i,c),((a,d),(f,(b,((g,(e,j)),h)))));((i,h),((c,f),((a,(d,e)),(j,(g,b)))));(((f,(b,d)),((j,g),(e,i))),(c,(h,a)));(((f,(i,(j,d))),(c,b)),(((h,a),g),e));((((h,e),((c,f),a)),(d,b)),((g,i),j))" -R8 -q1 -os
   i: ((((#2,(((((e)#6,a),#3))#1,((c)#5,b))))#7,((((j)#4,h))#3,((#7,(#5,(((i,#4))#8,((d)#2,f)))),#8))),(#6,(#1,g))) cost=43
   >: (((((((#2,((e)#6,a)),#3))#1,((c)#5,b)))#7,((((j)#4,h))#3,((#7,(#5,(((i,#4))#8,((d)#2,f)))),#8))),(#6,(#1,g))) cost=39
   >: ((((((((e)#6,(#2,a)),#3))#1,((c)#5,b)))#7,((((j)#4,h))#3,((#7,(#5,(((i,#4))#8,((d)#2,f)))),#8))),(#6,(#1,g))) cost=38
   >: ((((((((e)#6,(#2,a)),#3))#1,((c)#5,b)))#7,((((j)#4,h))#3,((#7,(((i,#4))#8,((d)#2,(#5,f)))),#8))),(#6,(#1,g))) cost=37
   >: ((((((((e)#6,(#2,a)),#3))#1,((c)#5,b)))#7,((((j)#4,h))#3,(#7,(((i,#4))#8,((d)#2,(#5,f)))))),(#6,(#1,(g,#8)))) cost=32
   >: ((((((((e)#6,(#2,a)),#3))#1,((c)#5,b)))#7,((((j)#4,h))#3,((((#7,i),#4))#8,((d)#2,(#5,f))))),(#6,(#1,(g,#8)))) cost=30
   >: ((((((((e)#6,(#2,a)),#3))#1,((c)#5,b)))#7,((((j)#4,h))#3,(((i,#4))#8,((d)#2,(#5,f))))),(#6,(#1,((#7,g),#8)))) cost=27
   >: ((((((e)#6,(#2,a)),#3))#1,((((c)#5,b))#7,((((j)#4,h))#3,(((i,#4))#8,((d)#2,(#5,f)))))),(#6,(#1,((#7,g),#8)))) cost=26
   >: ((((((e)#6,(#2,a)),#3))#1,((((c)#5,b))#7,(((((((j)#4,h))#3,i),#4))#8,((d)#2,(#5,f))))),(#6,(#1,((#7,g),#8)))) cost=25
   >: ((((((((((j)#4,h))#3,e))#6,(#2,a)),#3))#1,((((c)#5,b))#7,(((i,#4))#8,((d)#2,(#5,f))))),(#6,(#1,((#7,g),#8)))) cost=24
   >: ((((((((((j)#4,h))#3,e))#6,(#2,a)),#3))#1,((((c)#5,b))#7,(((i,#4))#8,((d)#2,(#5,f))))),(#1,((#7,(#6,g)),#8))) cost=23
   >: (((((c)#5,b))#7,(((i,#4))#8,((d)#2,(((((((((j)#4,h))#3,e))#6,(#2,a)),#3))#1,(#5,f))))),(#1,((#7,(#6,g)),#8))) cost=22
   >: (((((c)#5,b))#7,(((i,#4))#8,((d)#2,(((((((h)#3,e))#6,(#2,a)),#3))#1,(#5,f))))),(#1,((#7,(#6,((j)#4,g))),#8))) cost=21
   >: (((((c)#5,b))#7,((i)#8,(((d,#4))#2,(((((((h)#3,e))#6,(#2,a)),#3))#1,(#5,f))))),(#1,((#7,(#6,((j)#4,g))),#8))) cost=20
   >: ((((((d,#4))#2,((c)#5,b)))#7,((i)#8,(((((((h)#3,e))#6,(#2,a)),#3))#1,(#5,f)))),(#1,((#7,(#6,((j)#4,g))),#8))) cost=19
1. Cost:19 Steps:2360 Climbs:15 TopNetworks:16 Class:TreeChild
NewOptCost: 19
Cost:19 Steps:2360 Climbs:15 TopNetworks:16 Class:TreeChild HCruns:1 HCTime:0.776569 Naive:2361 Naivetime:0.770168 DTcnt:604416
```

Print improvements; skip odt.log:
```  
> supnet -g "(a,(b,(c,d))); ((a,b),(c,d))" -r1 -R3 -oNq
   i: (((b)#1,(((#1,c))#2,((d)#3,(a,#3)))),#2) cost=4
   >: (((b)#1,((a,#3),(((#1,c))#2,(d)#3))),#2) cost=2
   >: (((a,#3),((((#1,c))#2,(d)#3),(b)#1)),#2) cost=1
Cost:1 Steps:19 Climbs:3 TopNetworks:6 Class:TreeChild HCruns:1 HCTime:0.000133753 Naive:20 Naivetime:4.1008e-05 DTcnt:160
```

Recommended with large HC-runs using quasi-consensus rand networks
```  
> supnet -g "(a,(b,(c,d))); ((a,b),(c,d))" -q1000 -R3 -oS | tail -20
   i: (((((#2,((((a)#2,d))#3,b)))#1,c),#3),#1) cost=1
   i: ((((#1,b))#3,(#2,((a)#1,((d)#2,c)))),#3) cost=1
   >: (((#1,b))#3,(#2,((a)#1,(((d)#2,c),#3)))) cost=0
   i: ((((#3,((c)#3,d)))#1,((b)#2,(a,#2))),#1) cost=1
   >: ((((#3,(((c)#3,d),#2)))#1,((b)#2,a)),#1) cost=0
   i: (#1,((((#3,((((c)#2,b))#3,a)))#1,d),#2)) cost=1
   >: (#1,(#3,(((((((c)#2,b))#3,a))#1,d),#2))) cost=0
   i: ((((#2,((((b)#2,c))#1,a)))#3,(#1,d)),#3) cost=1
   >: (((b)#2,(((#2,((c)#1,a)))#3,(#1,d))),#3) cost=0
   i: (((((d)#2,(#2,c)))#3,(#1,((b)#1,a))),#3) cost=1
   >: (((((d)#2,(#1,(#2,c))))#3,((b)#1,a)),#3) cost=0
   i: (#3,((((((d)#1,b))#2,(#2,a)))#3,(#1,c))) cost=1
   >: (#3,(#2,((((((d)#1,b))#2,a))#3,(#1,c)))) cost=0
   i: ((#1,(((#3,((((d)#3,b))#1,c)))#2,a)),#2) cost=1
   >: ((#1,((((((d)#3,b))#1,(#3,c)))#2,a)),#2) cost=0
   i: (((b,#2))#1,(((a)#2,(#3,((c)#3,d))),#1)) cost=1
   >: (((b,#2))#1,((a)#2,(#3,(((c)#3,d),#1)))) cost=0
   i: ((((((c)#3,(#2,a)))#1,((b)#2,d)),#3),#1) cost=1
   >: (((((c)#3,(#2,a)))#1,((b)#2,(d,#3))),#1) cost=0
Cost:0 Steps:37567 Climbs:1649 TopNetworks:254 Class:TreeChild HCruns:1000 HCTime:0.179766 Naive:38567 Naivetime:0.056159 DTcnt:248146
```

Using stopping criterion (-K1000) with quasi-consensus rand networks
```  
> supnet -g "(a,(b,(c,d))); ((a,b),(c,d))" -q-1 -R3 -oS -K1000 | tail -20
   i: (#3,(((((a)#2,((c)#1,(d,#1))))#3,b),#2)) cost=0
   i: (#3,((((((c,#2))#3,b))#1,((d)#2,a)),#1)) cost=0
   i: ((((c)#1,(((d,#1))#2,((a)#3,b))),#3),#2) cost=0
   i: (((#1,d))#3,(#3,(#2,((a)#2,((c)#1,b))))) cost=1
   >: (((#1,d))#3,(#2,(#3,((a)#2,((c)#1,b))))) cost=0
   i: ((((((b)#2,d))#1,(((c,#1))#3,a)),#2),#3) cost=1
   >: ((((d)#1,((((b)#2,(c,#1)))#3,a)),#2),#3) cost=0
   i: (#2,((c)#1,(((((d,#1))#2,a))#3,(#3,b)))) cost=1
   >: ((c)#1,(((((d,#1))#2,a))#3,(#2,(#3,b)))) cost=0
   i: ((#1,(((#3,d))#2,((c)#3,((b)#1,a)))),#2) cost=1
   >: ((((#1,(#3,d)))#2,((c)#3,((b)#1,a))),#2) cost=0
   i: ((a)#2,(#3,(((((d)#1,(#1,c)))#3,b),#2))) cost=0
   i: (((#1,a))#2,(#2,((((b)#1,d))#3,(c,#3)))) cost=1
   >: (((#1,a))#2,(#2,((b)#1,((d)#3,(c,#3))))) cost=0
   i: ((#3,((((((a)#2,c))#3,(#2,b)))#1,d)),#1) cost=1
   i: (#2,(((a)#1,((((c)#3,(#1,b)))#2,d)),#3)) cost=1
   >: (#2,((a)#1,((((c)#3,(#1,b)))#2,(d,#3)))) cost=0
   i: (((c,#2))#1,((((((a)#3,d))#2,b),#3),#1)) cost=1
   >: (((c,#2))#1,(#3,(((((a)#3,d))#2,b),#1))) cost=0
Cost:0 Steps:47412 Climbs:2075 TopNetworks:254 Class:TreeChild HCruns:1287 HCTime:0.232986 Naive:48699 Naivetime:0.0721416 DTcnt:312161
```

## Evaluation parameter `-e`

Usage: `-e [gsrD12...]+`

    - `g`: print a gene tree
    - `s`: print a species tree
    - `n`: print a network
    - `t`: print all display trees
    - `T`: print all display trees with their ids       
    - `c`: print cost between two trees (G,S)\n"       
    - `C`: print cost between two trees with trees (G,S)\n"      
    - `i`: print species dictionary
    - `a`: odtdatfile in labeled format
    - `_`: print subtrees of species trees (-s/-S)
    - `L`: for each v in V(N), print the number of nodes reachable from v (only from networks)
    - `l`: for each v in V(N), print the number of leaves reachable from v (only from networks)

    - `p`: pairwise comparison of all networks
    - `u`: print unique networks from the input; summary stats printed on stderr
    - `U`: similar to u plus counts of networks
    - `x`: two networks are equal if their shapes are isomorphic in u,U,p (i.e., ignore leaf labels)

    - `d`: run DP algorithm only (approx DC) per each (G,N)
    - `b`: run BB algorithm only (exact DC); see also jJ subopt.      
    - `o`: min total cost of all gene trees vs network using a naive approach


### Networks classes set by parameter `-e`
    - `1`: class 1 (relaxed) networks in network generators and hill climbing moves
    - `2`: general phylogenetic networks in network generators and hill climbing moves

    Default class is 0 = tree-child.


### Branch and Bound algorithm for DC (with -eb)

 - `-ebJ` - gen bb.tsv with stats
 - `-ebj` - gen bb.dot with bb tree search
 - `-ebk` - print time per each pair

### Preserve root when searching the species tree/network space and in quasi-consensus `-er`: 

Print 5 quasi consensus tree with preserved split of the root.
```
> supnet -g '(a,((b,c),d));(a,(b,e))' -s'((a,b),(c,(d,e)))' -q5 -enr
((b,a),(e,(c,d)))
(((e,d),c),(b,a))
((a,b),((e,c),d))
((e,(d,c)),(b,a))
((e,(d,c)),(a,b))
```

### Print unique networks from odt.log `-eu`

```  
> supnet -g '(a,(b,(c,d))); ((a,b),(c,d))' -r1 -R3 -o2 && supnet -eu -N odt.log
   i: ((c)#3,((d)#2,((((b,#2))#1,(#3,a)),#1))) cost=3
   >: ((c)#3,((d)#2,(((((#3,b),#2))#1,a),#1))) cost=2
   >: ((c)#3,(((#3,d))#2,((((b,#2))#1,a),#1))) cost=0
Cost:0 Steps:32 Climbs:3 TopNetworks:8 Class:TreeChild HCruns:1 HCTime:0.00017333 Naive:33 Naivetime:5.10216e-05 DTcnt:202
((c)#3,(((#3,d))#2,((((b,#2))#1,a),#1)))
((c)#3,(((((b,#2))#1,a),((#3,d))#2),#1))
(((c)#3,((((b,#2))#1,a),#1)),((#3,d))#2)
(((#3,((c)#3,d)))#2,((((b,#2))#1,a),#1))
(((#3,d))#2,((((((c)#3,b),#2))#1,a),#1))
(((#3,d))#2,(((c)#3,(((b,#2))#1,a)),#1))
(((#3,d))#2,((((b,#2))#1,((c)#3,a)),#1))
(((#3,d))#2,(((((c)#3,(b,#2)))#1,a),#1))
unique=8 all=8
```

## Dot representation

```
supnet -n '(((b)#A,a),(#A,c))' -d
```

```
>  supnet -n '((a)#A,(#A,(c,b)))' -d | dot -Tpdf > n.pdf
```

```
>  supnet -n '((((c)#B,b))#A,(#A,(#B,a)))' -d | dot -Tpdf > n.pdf
```

```
>  supnet  -r1 -R10  -A20 -d | dot -Tpdf > n.pdf
```

## Debug tree/network info 

Detailed tree/network info (debug level).
```
>  supnet -g '(a,((b,a),c))' -eD
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

Print debug info on network structure
```
  supnet -n '((a)#A,(#A,(c,b)))' -enD
```

```
>  supnet -n '((((c)#B,b))#A,(#A,(#B,a)))' -enD
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

## More examples

Print min total cost 10 random gene trees vs random tree-child network with 5 reticulations over 8 species; print the initial network

```
> supnet -r10 -A8 -en  | supnet -G- -r1 -A8 -R5 -eon
((((h)#1,((a)#4,(c,(((g,#4))#5,e)))),(#1,((b)#2,(((#5,d))#3,(#2,f))))),#3)
62
```

Display tree usage stats:
```  
> supnet -N odt.log -et | sort | uniq -c | sort -k1 -n
      1 (((c,d),a),b)
      1 (c,((d,a),b))
      2 (((c,b),d),a)
      2 ((c,(d,b)),a)
      3 (((c,a),b),d)
      3 ((c,a),(d,b))
      3 ((c,(b,a)),d)
      3 (c,(d,(b,a)))
      4 (((c,b),a),d)
      7 (c,((d,b),a))
     17 ((c,d),(b,a))
     18 (((c,d),b),a)
```

Pairwise comparison of random dags (general dags)
```  
> supnet -r10 -R1 -A2 -en2p
(#1,((a)#1,b))
(((a)#1,b),#1)
(((a)#1,b),#1)
(((b)#1,a),#1)
(#1,((b)#1,a))
((a)#1,(b,#1))
(#1,((a)#1,b))
((a)#1,(b,#1))
(((b)#1,a),#1)
(#1,((b)#1,a))

(#1,((a)#1,b))	(((a)#1,b),#1)	1
(#1,((a)#1,b))	(((a)#1,b),#1)	1
(#1,((a)#1,b))	(((b)#1,a),#1)	0
(#1,((a)#1,b))	(#1,((b)#1,a))	0
(#1,((a)#1,b))	((a)#1,(b,#1))	1
(#1,((a)#1,b))	(#1,((a)#1,b))	1
(#1,((a)#1,b))	((a)#1,(b,#1))	1
(#1,((a)#1,b))	(((b)#1,a),#1)	0
(#1,((a)#1,b))	(#1,((b)#1,a))	0

(((a)#1,b),#1)	(((a)#1,b),#1)	1
(((a)#1,b),#1)	(((b)#1,a),#1)	0
(((a)#1,b),#1)	(#1,((b)#1,a))	0
(((a)#1,b),#1)	((a)#1,(b,#1))	1
(((a)#1,b),#1)	(#1,((a)#1,b))	1
(((a)#1,b),#1)	((a)#1,(b,#1))	1
(((a)#1,b),#1)	(((b)#1,a),#1)	0
(((a)#1,b),#1)	(#1,((b)#1,a))	0

(((a)#1,b),#1)	(((b)#1,a),#1)	0
(((a)#1,b),#1)	(#1,((b)#1,a))	0
(((a)#1,b),#1)	((a)#1,(b,#1))	1
(((a)#1,b),#1)	(#1,((a)#1,b))	1
(((a)#1,b),#1)	((a)#1,(b,#1))	1
(((a)#1,b),#1)	(((b)#1,a),#1)	0
(((a)#1,b),#1)	(#1,((b)#1,a))	0

(((b)#1,a),#1)	(#1,((b)#1,a))	1
(((b)#1,a),#1)	((a)#1,(b,#1))	0
(((b)#1,a),#1)	(#1,((a)#1,b))	0
(((b)#1,a),#1)	((a)#1,(b,#1))	0
(((b)#1,a),#1)	(((b)#1,a),#1)	1
(((b)#1,a),#1)	(#1,((b)#1,a))	1

(#1,((b)#1,a))	((a)#1,(b,#1))	0
(#1,((b)#1,a))	(#1,((a)#1,b))	0
(#1,((b)#1,a))	((a)#1,(b,#1))	0
(#1,((b)#1,a))	(((b)#1,a),#1)	1
(#1,((b)#1,a))	(#1,((b)#1,a))	1

((a)#1,(b,#1))	(#1,((a)#1,b))	1
((a)#1,(b,#1))	((a)#1,(b,#1))	1
((a)#1,(b,#1))	(((b)#1,a),#1)	0
((a)#1,(b,#1))	(#1,((b)#1,a))	0

(#1,((a)#1,b))	((a)#1,(b,#1))	1
(#1,((a)#1,b))	(((b)#1,a),#1)	0
(#1,((a)#1,b))	(#1,((b)#1,a))	0

((a)#1,(b,#1))	(((b)#1,a),#1)	0
((a)#1,(b,#1))	(#1,((b)#1,a))	0

(((b)#1,a),#1)	(#1,((b)#1,a))	1

```

Print unique random networks with counts
```  
> supnet -r100000 -R1 -A3 -eU
4512	(((b)#1,a),(#1,c))
3905	((a)#1,(#1,(b,c)))
5746	((#1,((c)#1,b)),a)
5045	((((c)#1,a),b),#1)
4441	((#1,a),((c)#1,b))
4935	((((b)#1,a),c),#1)
4902	((((b)#1,c),a),#1)
4478	((b,#1),((a)#1,c))
5009	(#1,(((a)#1,c),b))
3863	((((c,b))#1,a),#1)
3858	((((a,c))#1,b),#1)
5561	(c,((b)#1,(a,#1)))
5553	((#1,((a)#1,b)),c)
3897	(((c)#1,(a,b)),#1)
5669	(((b)#1,(c,#1)),a)
5587	(b,((c)#1,(#1,a)))
4856	(#1,(((a)#1,b),c))
4907	((c)#1,((#1,b),a))
3842	((((a,b))#1,c),#1)
3840	((b)#1,((c,a),#1))
5594	(b,(#1,((a)#1,c)))
unique=21 all=0
```

Print unique random shapes of networks with counts
```  
> supnet -r100000 -R1 -A3 -eUx
11732	(#1,(((c,b))#1,a))
33471	(c,((b)#1,(#1,a)))
11591	(((b)#1,(c,a)),#1)
30021	(#1,(((c)#1,a),b))
13185	((c,#1),((a)#1,b))
unique=5 all=0
```

Generate the picture of BB-tree search:

```  
> supnet -r1 -A15 -en  | supnet -G- -r1 -A15 -R12 -ebj; dot bb.dot -Tpdf > bb.pdf
-4
```

## Development variants

### Display tree cache 

When compiling with macro DTCACHE `make supnet_dtcache` option `-c MAXTREECACHESIZE` will limit the number of stored nodes in the cache. The default is `1000000` nodes.

```
make supnet_dtcache
supnet_dtcache -c 5000000 ...
```

### Limit the number of improvements

The number of HC steps in improvements 
```
make supnet_dtcache
supnet_dtcache -c 5000000 ...
```

