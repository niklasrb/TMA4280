CC = g++
CFLAGS = -Wall -g
LDFLAGS = -O3

OBJ = zeta0.o
BIN = zeta0

main: $(OBJ)
	$(CC) $(CFLAGS) -o $(BIN) $(OBJ) $(LDFLAGS)

%.o: %.cpp
	$(CC) $(CFLAGS) -c $<

vtest: main
	./$(BIN) -1

utest: main 
	./$(BIN) 0

.PHONY: clean
clean:
	rm -rf $(BIN) $(OBJ)
