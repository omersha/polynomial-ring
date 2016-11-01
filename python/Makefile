PROJ=polynomials
CC=g++

COMPILE_FLAGS=--std=c++17 -Wall -O3 -c -m64 -fPIC
LINK_FLAGS=-shared -Wl,-soname,$(PROJ).so

INC=-I ../

$(PROJ): python.cpp python.h
	$(CC) $(COMPILE_FLAGS) python.cpp $(INC)
	$(CC) $(LINK_FLAGS) python.o -o $(PROJ).so

.PHONY: clean

clean:
	rm $(PROJ).so python.o


