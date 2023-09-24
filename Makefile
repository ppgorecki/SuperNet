#VALGRIND=valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes --tool=memcheck
VALGRIND=valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes --tool=memcheck 
CPPFLAGS = -Ofast -fomit-frame-pointer -g 
CC = g++ -Ofast -fomit-frame-pointer -g 
MAKEFLAGS += -j10 # parallel

.PHONY: all supnet clean_progx progx

all: supnet

OBJS=tools.o clusters.o dag.o rtree.o bb.o bbstats.o network.o hillclimb.o supnet.o iso.o dp.o contrnet.o bitcluster.o treespace.o hcstats.o odtstats.o topsort.o randnets.o

SRC=tools.cpp clusters.cpp dag.cpp rtree.cpp dp.cpp bb.cpp bbstats.cpp network.cpp hillclimb.cpp supnet.cpp iso.cpp contrnet.cpp bitcluster.cpp treespace.cpp odtstats.cpp hcstats.cpp topsort.cpp randnets.cpp

clusters.o: clusters.cpp clusters.h tools.h rtree.h dag.h network.h bb.h \
 treespace.h bitcluster.h

tools.o: tools.cpp tools.h bb.h

rtree.o: rtree.cpp tools.h clusters.h rtree.h dag.h network.h bb.h \
 treespace.h bitcluster.h costs.h

network.o: network.cpp rtree.h tools.h clusters.h dag.h network.h bb.h \
 treespace.h bitcluster.h dp.h costs.h

dag.o: dag.cpp tools.h dag.h

iso.o: iso.cpp tools.h dag.h

hillclimb.o: hillclimb.cpp hillclimb.h rtree.h tools.h clusters.h dag.h \
 network.h bb.h treespace.h bitcluster.h dagset.h costs.h

dp.o: dp.cpp rtree.h tools.h clusters.h dag.h network.h bb.h treespace.h \
 bitcluster.h contrnet.h dp.h


bb.o: bb.cpp bb.h tools.h rtree.h clusters.h dag.h network.h treespace.h \
 bitcluster.h costs.h contrnet.h dp.h

contrnet.o: contrnet.h network.h dag.h tools.h bb.h rtree.h clusters.h \
 bitcluster.h treespace.h

treespace.o: treespace.cpp treespace.h bitcluster.h tools.h clusters.h \
 rtree.h dag.h network.h bb.h

supnet.o: supnet.cpp tools.h bb.h clusters.h rtree.h dag.h network.h \
 treespace.h bitcluster.h contrnet.h hillclimb.h dagset.h costs.h

supnet: $(OBJS)
	$(CC) $(LFLAGS) -o $@ $^

#supnet_dtcache: 
#g++ -D DTCACHE -O3 -o -$@ ${SRC}

clean_dtcache:
	rm -f network.o
supnet_dtcache: clean_dtcache $(OBJS)
	$(MAKE) _supnet_dtcache
_supnet_dtcache: network.o
_supnet_dtcache: CPPFLAGS += -D DTCACHE
_supnet_dtcache: $(OBJS) 
	$(CC) $(LFLAGS) -o supnet_dtcache ${OBJS}
	
gdb: clean  
	g++ -g -o supnet ${SRC}

gdb_dtcache: clean  
	g++ -D DTCACHE -g -o supnet_dtcache ${SRC}

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

profiler: 
	g++ -Wall -pg -O3 -o supnet ${SRC}	
	supnet -r100 -R2 -A10 -g "((a,b),(c,d))" -ed
	gprof supnet

profilerlines: 
	rm -f *gcov *gcda *gcno
	g++ -Wall -pg -fprofile-arcs -ftest-coverage -o supnet ${SRC}
	supnet -N /tmp/hcnet.txt -G /tmp/hcgtr.txt -oT1t -CDC
	#supnet -r1000 -R2 -A10 -g "((a,b),(c,d))" -ed
	gcov supnet


