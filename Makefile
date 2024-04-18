.PHONY = all clean

DIR = lib/
BIN = bin/
SRC = src/
CC = g++
FFTW = -L./lib/fftw-3.3.10/.libs -l:libfftw3.a -lm 
FLINT = -lflint -I/usr/local/lib -L/usr/local/lib
PY = -I/usr/include/python3.10 -lpython3.10
all: $(BIN)randgen.o $(BIN)rtest

$(BIN)randgen.o: $(SRC)randgen.cpp $(SRC)randgen.h 
	$(CC) -Wall -g -c -O3 $(SRC)randgen.cpp -o $(BIN)randgen.o $(FLINT) $(PY)

$(BIN)rtest: $(SRC)rtest.cpp $(BIN)randgen.o
	$(CC) -Wall -g -O3 $(SRC)rtest.cpp $(BIN)randgen.o -o $(BIN)rtest $(FLINT) $(PY) $(FFTW)


clean:
	@echo "Cleaning up..."
	rm -f bin/*