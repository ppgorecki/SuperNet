#VALGRIND=valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes --tool=memcheck
VALGRIND=valgrind --show-leak-kinds=all --track-origins=yes --tool=memcheck 
CPPFLAGS = -O3 
CC = g++ 
LFLAGS =  
MAKEFLAGS += -j10 # parallel

all: supnet

SRC=tools.cpp clusters.cpp dag.cpp rtree.cpp dp.cpp bb.cpp bbstats.cpp network.cpp hillclimb.cpp supnet.cpp iso.cpp contrnet.cpp

clusters.o: clusters.h clusters.cpp tools.h
tools.o : tools.cpp tools.h clusters.h 
rtree.o : rtree.cpp rtree.h tools.h clusters.h dag.h
network.o : network.cpp network.h dag.h rtree.h tools.h
dag.o: dag.cpp dag.h tools.cpp tools.h dagset.h
iso.o: iso.cpp dag.cpp dag.h tools.cpp tools.h
hillclimb.o: hillclimb.h hillclimb.cpp tools.h dagset.h
dp.o: dp.cpp network.cpp network.h rtree.h rtree.cpp
bb.o: bb.cpp dp.cpp network.cpp network.h rtree.h rtree.cpp bbstats.cpp
contrnet.o: contrnet.cpp contrnet.h network.h
supnet: tools.o clusters.o dag.o rtree.o bb.o bbstats.o network.o hillclimb.o supnet.o iso.o dp.o contrnet.o bbstats.o
	$(CC) $(LFLAGS) -o $@ $^

gdb: clean  
	g++ -g -o supnet ${SRC}

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
	tar czvf supnet.tgz ${SRC} Makefile README

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


