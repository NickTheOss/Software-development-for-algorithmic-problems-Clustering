CPP=g++

objects=cluster.o hashtable.o

SOURCE = cluster.cpp hashtable.cpp

HEADER = cluster.h  hashtable.h

OUT = cluster

FLAGS = -g -c 

all : $(objects)
	$(CPP) -g $(objects) -o $(OUT) 

cluster.o : cluster.cpp
	$(CPP) $(FLAGS) cluster.cpp

	
clean:
	rm -f cluster $(objects)

