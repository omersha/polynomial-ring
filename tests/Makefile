PROJ=tests
CC=g++

CFLAGS=--std=c++17 -Wall -O3 -m64 
INC=-I ../

$(PROJ): tests.cpp tests.h
	$(CC) $(CFLAGS) tests.cpp $(INC) -o tests

.PHONY: clean

clean:
	rm $(PROJ)


