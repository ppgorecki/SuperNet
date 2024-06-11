#VALGRIND=valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes --tool=memcheck
VALGRIND=valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes --tool=memcheck 
CPPFLAGS = -std=c++17 -Ofast -fomit-frame-pointer
#CPPFLAGS = -std=c++17 -g 
CC = g++ -Ofast -fomit-frame-pointer
MAKEFLAGS += -j10 # parallel

.ONESHELL:

.PHONY: all supnet clean_progx progx

all: supnet

OBJS=tools.o clusters.o dag.o rtree.o bb.o bbstats.o network.o hillclimb.o supnet.o iso.o dp.o contrnet.o bitcluster.o treespace.o hcstats.o odtstats.o topsort.o randnets.o testers.o neteditop.o

SRC=tools.cpp clusters.cpp dag.cpp rtree.cpp dp.cpp bb.cpp bbstats.cpp network.cpp hillclimb.cpp supnet.cpp iso.cpp contrnet.cpp bitcluster.cpp treespace.cpp odtstats.cpp hcstats.cpp topsort.cpp randnets.cpp testers.cpp neteditop.o



supnet: $(OBJS)
	$(CC) $(LFLAGS) -o $@ $^

clean_nodtcache:
	rm -f network.o

supnet_nodtcache: clean_nodtcache $(OBJS)
	$(MAKE) _supnet_nodtcache

_supnet_nodtcache: network.o
_supnet_nodtcache: CPPFLAGS += -D NODTCACHE
_supnet_nodtcache: $(OBJS) 
	$(CC) $(LFLAGS) -o supnet_nodtcache ${OBJS}

tools.o: tools.cpp tools.h bb.h clusters.h network.h dag.h rtree.h \
 bitcluster.h treespace.h stats.h dagset.h
clusters.o: clusters.cpp clusters.h tools.h rtree.h dag.h network.h bb.h \
 treespace.h bitcluster.h stats.h dagset.h
dag.o: dag.cpp tools.h dag.h clusters.h topsort.h
rtree.o: rtree.cpp tools.h clusters.h rtree.h dag.h network.h bb.h \
 treespace.h bitcluster.h stats.h dagset.h costs.h
dp.o: dp.cpp rtree.h tools.h clusters.h dag.h network.h bb.h treespace.h \
 bitcluster.h stats.h dagset.h contrnet.h dp.h
bb.o: bb.cpp bb.h tools.h rtree.h clusters.h dag.h network.h treespace.h \
 bitcluster.h stats.h dagset.h costs.h contrnet.h dp.h
bbstats.o: bbstats.cpp bb.h tools.h rtree.h clusters.h dag.h network.h \
 treespace.h bitcluster.h stats.h dagset.h contrnet.h dp.h
network.o: network.cpp rtree.h tools.h clusters.h dag.h network.h bb.h \
 treespace.h bitcluster.h stats.h dagset.h dp.h costs.h
hillclimb.o: hillclimb.cpp hillclimb.h rtree.h tools.h clusters.h dag.h \
 network.h bb.h treespace.h bitcluster.h stats.h dagset.h costs.h \
 randnets.h
supnet.o: supnet.cpp bb.h tools.h bitcluster.h clusters.h contrnet.h \
 network.h dag.h rtree.h treespace.h stats.h dagset.h costs.h hillclimb.h \
 randnets.h topsort.h testers.h
iso.o: iso.cpp tools.h dag.h clusters.h
contrnet.o: contrnet.cpp rtree.h tools.h clusters.h dag.h network.h bb.h \
 treespace.h bitcluster.h stats.h dagset.h contrnet.h
bitcluster.o: bitcluster.cpp bitcluster.h tools.h clusters.h
treespace.o: treespace.cpp treespace.h bitcluster.h tools.h clusters.h \
 rtree.h dag.h network.h bb.h stats.h dagset.h
odtstats.o: odtstats.cpp stats.h bb.h tools.h dagset.h rtree.h clusters.h \
 dag.h network.h treespace.h bitcluster.h
hcstats.o: hcstats.cpp stats.h bb.h tools.h dagset.h rtree.h clusters.h \
 dag.h network.h treespace.h bitcluster.h hillclimb.h costs.h randnets.h
topsort.o: topsort.cpp dag.h tools.h clusters.h topsort.h
randnets.o: randnets.cpp network.h dag.h tools.h clusters.h bb.h rtree.h \
 bitcluster.h treespace.h stats.h dagset.h randnets.h
testers.o: testers.cpp testers.h network.h dag.h tools.h clusters.h bb.h \
 rtree.h bitcluster.h treespace.h stats.h dagset.h contrnet.h randnets.h \
 neteditop.h


gdb: clean  
	g++ -g -o supnet ${SRC}

gdb_nodtcache: clean  
	g++ -D NODTCACHE -g -o supnet_nodtcache ${SRC}

supnet_bfs:   
	g++ -D USE_QUEUE_BFS -O3 -o supnet_bfs ${SRC}

supnet_minrt:  
	g++ -D USE_PRIORITY_QUEUE_MINRT -O3 -o supnet_minrt ${SRC}

supnet_minlb:  
	g++ -D USE_PRIORITY_QUEUE_MINLB -O3 -o supnet_minlb ${SRC}

clean:
	rm -f supnet supnet_bfs supnet_minlb supnet_minrt *.o  *.old *~ x *.log *.gcov *.gcda *.gcno

o1g: clean 
	g++ -g -O1 -o supnet ${SRC}

o3: clean 
	g++ -g -O3 -o supnet ${SRC}
	
tgz: 
	tar czvf supnet.tgz ${SRC} *.h Makefile README

check-syntax:
	gcc -o nul -S ${CHK_SOURCES}

valgrinddp: gdb 
	${VALGRIND} supnet -r100 -R2 -A3 -g "(a,(b,c))" -ed

valgrinddpbb: gdb
	${VALGRIND} supnet -g "(a,(b,(c,d))); ((a,b),(c,d))" -r1 -R2 -oTe3

valgrinddc: gdb
	${VALGRIND} supnet -n "(#2,((#1)#2,((d,((a)#1,b)),c)))" -g "(a,(b,(c,d)))" -eY

valgrinddc2: gdb
	${VALGRIND} supnet -r1 -R3 -A4 -g "(a,(b,(c,d)))" -eY

valgrindhc: gdb
	${VALGRIND} supnet -N /tmp/hcnet.txt -G /tmp/hcgtr.txt -oT1t -CDC

valgrind80: gdb
	${VALGRIND} supnet -g '(((((t5,(((((t31,t64),t19),(t55,t51)),((((t62,(t13,t43)),(t59,t76)),t44),(t53,t10))),(((t4,((t20,t40),t71)),t39),t37))),(t11,t54)),(t77,(t69,t28))),(((((t45,(t6,t79)),t15),(t65,(t32,(t27,t73)))),((((((t16,t35),t17),((((t56,t23),(t61,t47)),(t60,t3)),(t48,t52))),t8),(t67,(t9,t70))),(((t80,t46),t29),(t72,t22)))),t58)),((((((t7,t7),((t1,(t41,t49)),t2)),(t36,t21)),(((t34,t34),(((t30,(t18,(t66,t33))),(t26,t57)),(t38,t12))),((t68,t50),t24))),(((((t14,t74),t42),(t42,(t74,t14))),(((((t18,(t57,t26)),t30),(t12,t38)),t34),((t68,t50),t24))),((((((t26,t57),(t12,t38)),((t66,t33),t30)),t18),((t50,t68),t24)),((t14,t74),t42)))),((t7,(t2,(t1,(t49,t41)))),((((t41,t49),t1),t2),t7))))' -n '(((t5,(((t11,t54),(#M,((((t59,t76),((#C,(#F,(t37,((t4,((((#S,t40),(t46)#B),(t20)#L))#A),t39)))))#M),(((t62,(t53,(t10)#Q)),(t44,(#P,t43))))#R),((t64,(#T,t19)),((t55,(t51)#H))#I)))),((t69,t28),t77))),(#Q,((((#G,t29),(#A,((#K,(t58,(t71)#S)),(#L,((t72,(t30)#E),t22))))),(((((t27,(t73)#F),t32),t65),((t15,((t6,t79),t45)),#O)),(((t8,((t35,(t13)#P),t16)),(t17,(#D,(((((t3,#H),t60),(((t61,(t52)#D),t47))#O),((t23,t56))#J),(t48,(t25)#C))))),(t67,(t70,(t9)#G))))),#N))),((t36,t21),(#R,((t7,(t2,(#J,(t1,(t41,t49))))),(((t34,((t80,(#B,(t75,(t78,t63)))))#K),(((t38,t12),((t57,t26),(#E,(#I,((t33,t66),t18))))),((t24,(t50,t68)))#N)),(t42,(((t31)#T,t74),t14)))))))' -t5 -ebk -CDC

valgrind85: gdb	
	${VALGRIND} --log-file=valgrind85.v2.q1.log supnet -G "gene_trees_45x.txt" -z 14 -O gene_trees_45_valgrind.odt -q1 -R5 -ot1 -CDC -m1

valgrind23: gdb_dtcache	
	${VALGRIND} --log-file=valgrind85.v2.q1.log supnet_dtcache -G "gene_trees_45.txt" -z 14 -O gene_trees_45_valgrind.odt -q1 -R10 -ot1 -CDC  -z10 -m2

valgrind45si: gdb	
	${VALGRIND} --log-file=valgrind45si.log supnet -G "gene_trees_45si.txt" -z 14 -O gene_trees_45_valgrind.odt -q1 -R10 -ot1 -CDC  -z10 -m1

valgrind45n: gdb	
	${VALGRIND} --log-file=valgrind45.log supnet -G "gene_trees_45.txt" -z 14 -O gene_trees_45_valgrind.odt -q1 -R10 -ot1 -CDC  -z10 -m1

valgrindgu: 	
	${VALGRIND} --log-file=valgrindgu.log supnet --guidetree "((a,b,c),(d,e,f))" -q1000 -R2

valgrinddg: 	
	${VALGRIND} --log-file=valgrinddg.log supnet -g "((i,c),((a,d),(f,(b,((g,(e,j)),h)))));((i,h),((c,f),((a,(d,e)),(j,(g,b)))));(((f,(b,d)),((j,g),(e,i))),(c,(h,a)));(((f,(i,(j,d))),(c,b)),(((h,a),g),e));((((h,e),((c,f),a)),(d,b)),((g,i),j))" -R8 -q1 --HC --hcrunstats

# profiler: 
# 	g++ -Wall -pg -O3 -o supnet ${SRC}	
# 	supnet -r100 -R2 -A10 -g "((a,b),(c,d))" -ed
# 	gprof supnet

# profiler2:
# 	g++ -D DTCACHE -Wall -pg -O3 -o supnet_dtcache ${SRC}	
# 	supnet_dtcache -G wheat_trees_clean -R6 -q10 --hcstopinit=200 -v1 --HC --hcrunstatsext  --relaxed --hcsavewhenimproved --odtlabelled --guidetree "(1,2,3,4);(5,6,7,8);(9,10,11,12);(13,14,15);(16,17,18);(19,20,21);(22,23,24);(25,26);(27,28,29);(30,31,32,33);(34,35);(36,37,38,39);(40,41,42,43)" --outfiles R6_wheat_relaxed_guidetree  --displaytreesampling="0.03 0.0625 0.125 0.25 0.5 1" --hcsamplerstats
# 	gprof supnet_dtcache

# profilerlines: 
# 	rm -f *gcov *gcda *gcno
# 	g++ -Wall -pg -fprofile-arcs -ftest-coverage -o supnet ${SRC}
# 	#supnet -N /tmp/hcnet.txt -G /tmp/hcgtr.txt -oT1t -CDC
# 	supnet -r1000 -R2 -A10 -g "((a,b),(c,d))" -ed
# 	gcov supnet.cpp


profilerlines: 
	mkdir -p $@
	rm -f *gcov *gcda *gcno 
	cp wheat_trees_clean ${SRC} *.h 	$@	
	cd $@
	g++ -Wall -pg -fprofile-arcs -ftest-coverage -o supnet ${SRC}
	supnet -G wheat_trees_clean -R6 -q1 --hcstopinit=200 -v1 --HC --hcrunstatsext  --relaxed --hcsavewhenimproved --odtlabelled --guidetree "(1,2,3,4);(5,6,7,8);(9,10,11,12);(13,14,15);(16,17,18);(19,20,21);(22,23,24);(25,26);(27,28,29);(30,31,32,33);(34,35);(36,37,38,39);(40,41,42,43)" --displaytreesampling="0.03 0.0625 0.125 0.25 0.5 1" --hcsamplerstats	
	gcov supnet-bb.cpp -l
	echo $@ completed


profiler:
	mkdir -p $@
	rm -f *gcov *gcda *gcno 
	cp wheat_trees_clean ${SRC} *.h $@	
	cd $@
	g++ -Wall -pg -O3 -o supnet ${SRC}	
	supnet -G wheat_trees_clean -R6 -q1 --hcstopinit=200 -v1 --HC --hcrunstatsext  --relaxed --hcsavewhenimproved --odtlabelled --guidetree "(1,2,3,4);(5,6,7,8);(9,10,11,12);(13,14,15);(16,17,18);(19,20,21);(22,23,24);(25,26);(27,28,29);(30,31,32,33);(34,35);(36,37,38,39);(40,41,42,43)" --displaytreesampling="0.03 0.0625 0.125 0.25 0.5" --hcsamplerstats
	gprof supnet
	echo $@ completed

profiler10:
	mkdir -p $@
	rm -f *gcov *gcda *gcno 
	cp wheat_trees_clean ${SRC} *.h $@	
	cd $@
	g++ -Wall -pg -O3 -o supnet ${SRC}	
	supnet -G wheat_trees_clean -R6 -q10 --hcstopinit=200 -v1 --HC --hcrunstatsext  --relaxed --hcsavewhenimproved --odtlabelled --guidetree "(1,2,3,4);(5,6,7,8);(9,10,11,12);(13,14,15);(16,17,18);(19,20,21);(22,23,24);(25,26);(27,28,29);(30,31,32,33);(34,35);(36,37,38,39);(40,41,42,43)" --displaytreesampling="0.03 0.0625 0.125 0.25 0.5" --hcsamplerstats
	gprof supnet
	echo $@ completed	
	

dep: 
	gcc -MM ${SRC}