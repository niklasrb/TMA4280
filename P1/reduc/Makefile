CC = mpiCC
MPI = mpirun
CFLAGS = -Wall -g
LDFLAGS = -O3

OBJ1 = zeta1.o
BIN1 = zeta1
OBJ2 = mach1.o
BIN2 = mach1

main: $(OBJ1) $(OBJ2)
	$(CC) $(CFLAGS) -o $(BIN1) $(OBJ1) $(LDFLAGS)
	$(CC) $(CFLAGS) -o $(BIN2) $(OBJ2) $(LDFLAGS)

%.o: %.cpp
	$(CC) $(CFLAGS) -c $<


.PHONY: clean
clean:
	rm -rf $(BIN1) $(OBJ1) $(BIN2) $(OBJ2)
