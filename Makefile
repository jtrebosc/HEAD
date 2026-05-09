HEAD_processing: HEAD_processing.cpp simpson_rw.o
	g++ HEAD_processing.cpp -o HEAD_processing -lgsl -lgslcblas -lm -O3 -march=native -fpermissive -fopenmp
simpson_rw.o:simpson_rw.cpp
	g++ -std=c++11 -c simpson_rw.cpp

