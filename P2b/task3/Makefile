CXX=mpicxx
INCDIR=../include
CXXFLAGS= -I$(INCDIR) --std=c++11 -g -Wall 
LFLAGS = -O2 -g -fopenmp

SRC =$(INCDIR)/tma.h
SRC+=$(wildcard $(INCDIR)/tma/*.h)
SRC+=$(wildcard $(INCDIR)/tma/*/*.h)

main: main.cpp $(SRC)
	@$(CXX) $(CXXFLAGS) $(LFLAGS) -o $@ $<

main.cpp: $(SRC)
	
run: main
	OMP_NUM_THREADS=4 mpirun -np 4 ./main 20

