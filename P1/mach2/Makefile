CC = g++
CFLAGS = -Wall -g -fopenmp --std=c++11
LDFLAGS = -O3

OBJ = mach2.o
BIN = mach2

main: $(OBJ)
	$(CC) $(CFLAGS) -o $(BIN) $(OBJ) $(LDFLAGS)

%.o: %.cpp
	$(CC) $(CFLAGS) -c $<

.PHONY: clean
clean:
	rm -rf $(BIN) $(OBJ)
