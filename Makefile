VALGRIND=valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes --tool=memcheck
CPPFLAGS = -O3 
CC = g++ 
LFLAGS =  

all: supnet

SRC=tools.cpp clusters.cpp dag.cpp rtree.cpp network.cpp dp.cpp hillclimb.cpp supnet.cpp iso.cpp 

clusters.o: clusters.h clusters.cpp tools.h
tools.o : tools.cpp tools.h clusters.h 
rtree.o : rtree.cpp rtree.h tools.h clusters.h dag.h
network.o : network.cpp network.h dag.h rtree.h tools.h
dag.o: dag.cpp dag.h tools.cpp tools.h dagset.h
iso.o: iso.cpp dag.cpp dag.h tools.cpp tools.h
hillclimb.o: hillclimb.h hillclimb.cpp tools.h dagset.h
dp.o: dp.cpp network.cpp network.h rtree.h rtree.cpp

supnet: tools.o clusters.o dag.o rtree.o network.o hillclimb.o supnet.o iso.o dp.o
	$(CC) $(LFLAGS) -o $@ $^

gdb: clean  
	g++ -g -o supnet ${SRC}

clean:
	rm -f *.o $(TARGET) *.old *~ x *.log *.gcov *.gcda *.gcno

o1g: clean 
	g++ -g -O1 -o supnet ${SRC}

o3: clean 
	g++ -g -O3 -o supnet ${SRC}
	
tgz: 
	tar czvf supnet.tgz ${SRC} Makefile README

check-syntax:
	gcc -o nul -S ${CHK_SOURCES}

valgrinddp: gdb 
	${VALGRIND} supnet -r100 -R2 -A3 -g "(a,(b,c))"-ed

profiler: 
	g++ -Wall -pg -O3 -o supnet ${SRC}	
	supnet -r100 -R2 -A10 -g "((a,b),(c,d))" -ed
	gprof supnet

profilerlines: 
	rm *gcov *gcda *gcno
	g++ -Wall -pg -fprofile-arcs -ftest-coverage -o supnet ${SRC}
	supnet -r1000 -R2 -A10 -g "((a,b),(c,d))" -ed
	gcov supnet


