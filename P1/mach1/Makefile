CC = mpiCC
MPI = mpirun
CFLAGS = -Wall -g
LDFLAGS = -O3

OBJ = mach1.o
BIN = mach1

main: $(OBJ)
	$(CC) $(CFLAGS) -o $(BIN) $(OBJ) $(LDFLAGS)

%.o: %.cpp
	$(CC) $(CFLAGS) -c $<

test: main
	$(MPI) -np 4 $(BIN) 100

vtest: main
	$(MPI) -np 4 $(BIN) -1

utest: main 
	$(MPI) -np 4 $(BIN) 0

.PHONY: clean
clean:
	rm -rf $(BIN) $(OBJ)
