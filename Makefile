CPPFLAGS = -O3 -g 
CC = g++ 
LFLAGS =  

all: supnet

tools.o : tools.cpp tools.h

supnet: supnet.o tools.o
	$(CC) $(LFLAGS) -o $@ $^

gdb: clean  
	g++ -g -o tools.cpp supnet.cpp 

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
