CXX=mpicxx
INCDIR=../include
CXXFLAGS=-O3 -I$(INCDIR) -std=c++11 -g -Wall

SRC =$(INCDIR)/tma.h
SRC+=$(wildcard $(INCDIR)/tma/*.h)
SRC+=$(wildcard $(INCDIR)/tma/*/*.h)

main: main.cpp $(SRC)
	@$(CXX) $(CXXFLAGS) -o $@ $<

main.cpp: $(SRC)
	
run: main
	mpirun -np 4 ./main 1000

