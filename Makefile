CPPFLAGS = -O3
CC = g++ 
LFLAGS =  

all: supnet

clusters.o: clusters.h clusters.cpp tools.h
tools.o : tools.cpp tools.h clusters.h 
rtree.o : rtree.cpp rtree.h tools.h clusters.h dag.h
network.o : network.cpp network.h dag.h rtree.h tools.h
dag.o: dag.cpp dag.h tools.cpp tools.h
iso.o: iso.cpp dag.cpp dag.h tools.cpp tools.h
hillclimb.o: hillclimb.h hillclimb.cpp tools.h

supnet: tools.o clusters.o dag.o rtree.o network.o hillclimb.o supnet.o iso.o
	$(CC) $(LFLAGS) -o $@ $^

gdb: clean  
	g++ -g -o supnet tools.cpp rtree.cpp dag.cpp network.cpp clusters.cpp hillclimb.cpp supnet.cpp iso.cpp

clean:
	rm -f *.o $(TARGET) *.old *~ x *.log

o1g: clean 
	g++ -g -O1 -o supnet tools.cpp supnet.cpp 

o3: clean 
	g++ -g -O3 -o supnet tools.cpp supnet.cpp 
	
tgz : 
	tar czvf supnet.tgz supnet.cpp tools.h tools.cpp Makefile 

tgzlin : 
	tar czvf supnet2.tgz supnet2 README 

check-syntax:
	gcc -o nul -S ${CHK_SOURCES}


tools.o : tools.cpp tools.h