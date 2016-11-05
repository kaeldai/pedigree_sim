CXX = g++
CFLAGS = -std=c++11 -Wno-deprecated-declarations 
INC_DIRS = -I /usr/local/eigen3 -I /usr/local/eigen3/unsupported -I ~/.local/gsl_22/include -I ~/.local/boost_1_60_0_gcc/include
LIB_DIRS = -L ~/.local/gsl_22/lib -L ~/.local/boost_1_60_0_gcc/lib
LIBS = -lhts -lgsl -lgslcblas -lm -lboost_program_options

all: simulator.cpp
	$(CXX) $(CFLAGS) $(INC_DIRS) $(LIB_DIRS) -o simulator simulator.cpp $(LIBS)

run:
	./simulator --mu 1e-5 --mu-somatic 1e-5 --region "chr22:20,000,000-21,000,000" --ped ceu.ped ~/Data/references/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa
