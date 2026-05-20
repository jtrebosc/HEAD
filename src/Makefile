C++ = g++
CFLAGS = -O3 -march=native -fpermissive
CLIBS = -lgsl -lgslcblas -lm -fopenmp

all: HEAD_processing HEAD_2D_processing

HEAD_processing: HEAD_processing.cpp simpson_rw.o
	$(C++) HEAD_processing.cpp -o ../bin/HEAD_processing $(CLIBS) $(CFLAGS)

HEAD_2D_processing: HEAD_2D_processing.cpp simpson_rw.o simpson_rw.hpp
	$(C++) HEAD_2D_processing.cpp simpson_rw.o -o ../bin/HEAD_2D_processing $(CLIBS) $(CFLAGS)

simpson_rw.o:simpson_rw.cpp  simpson_rw.hpp
	g++ -std=c++11 -c simpson_rw.cpp simpson_rw.hpp

clean:
	rm *.o *.gch
