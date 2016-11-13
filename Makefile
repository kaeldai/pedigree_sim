CXX = g++

CFLAGS = -std=c++11 -Wno-deprecated-declarations 
INC_DIRS = -I include -I /usr/local/include/eigen -I /usr/local/include/eigen/unsupported -I ~/.local/gsl_22/include -I ~/.local/boost_1_60_0_gcc/include
LIB_DIRS = -L ~/.local/gsl_22/lib -L ~/.local/boost_1_60_0_gcc/lib
LIBS = -lhts -lgsl -lgslcblas -lm -lboost_program_options

all: src/simulator.cpp
	$(CXX) $(CFLAGS) $(INC_DIRS) $(LIB_DIRS) -o build/simulator src/simulator.cpp $(LIBS)

run_trio:
	./build/simulator --mu 1e-5 --mu-somatic 1e-5 --region "chr21" --prefix "ceph_trio" --ped tests/trio.ped tests/grch38_chr21.fa

run_septet:
	./build/simulator --mu 1e-5 --mu-somatic 1e-5 --region "chr21" --prefix "ceph_septet" --ped tests/septet.ped tests/grch38_chr21.fa

run_ceph_full:
	./build/simulator --mu 1e-5 --mu-somatic 1e-5 --region "chr21" --prefix "ceph_full" --ped tests/ceph1463.ped tests/grch38_chr21.fa
