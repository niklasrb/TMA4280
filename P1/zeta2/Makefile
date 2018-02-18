CC = g++
CFLAGS = -Wall -g -fopenmp
LDFLAGS = -O3

OBJ = zeta2.o
BIN = zeta2

main: $(OBJ)
	$(CC) $(CFLAGS) -o $(BIN) $(OBJ) $(LDFLAGS)

%.o: %.cpp
	$(CC) $(CFLAGS) -c $<

.PHONY: clean
clean:
	rm -rf $(BIN) $(OBJ)
