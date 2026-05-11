HEAD_2D_processing: HEAD_2D_processing.cpp simpson_rw.o simpson_rw.hpp
	g++ HEAD_2D_processing.cpp simpson_rw.o -o HEAD_2D_processing -lgsl -lgslcblas -lm -O3 -march=native -fpermissive -fopenmp
simpson_rw.o:simpson_rw.cpp  simpson_rw.hpp
	g++ -std=c++11 -c simpson_rw.cpp simpson_rw.hpp

