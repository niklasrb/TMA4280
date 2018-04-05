CXX=mpicxx
INCDIR=../include
CXXFLAGS=-O3 -I$(INCDIR) -std=c++11

SRC =$(INCDIR)/tma.h
SRC+=$(wildcard $(INCDIR)/tma/*.h)
SRC+=$(wildcard $(INCDIR)/tma/*/*.h)

main: main.cpp $(SRC)
	@$(CXX) $(CXXFLAGS) -o $@ $<

main.cpp: $(SRC)
	
run: test
	@./test

clean:
	@rm -f test
