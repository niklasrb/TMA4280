CC = mpiCC
MPI = mpirun
CFLAGS = -Wall -g -fopenmp --std=c++11
LDFLAGS = -O3

OBJ1 = zeta.o
BIN1 = zeta
OBJ2 = mach.o
BIN2 = mach

main: $(OBJ1) $(OBJ2)
	$(CC) $(CFLAGS) -o $(BIN1) $(OBJ1) $(LDFLAGS)
	$(CC) $(CFLAGS) -o $(BIN2) $(OBJ2) $(LDFLAGS)

%.o: %.cpp
	$(CC) $(CFLAGS) -c $<


.PHONY: clean
clean:
	rm -rf $(BIN1) $(OBJ1) $(BIN2) $(OBJ2)
