# Makefile for simulation

LDFLAGS=-lm -lomp
LIBS=-L/opt/homebrew/opt/libomp/lib
CFLAGS=-Wall -Wno-deprecated -Ofast -march=native -std=c11 -Xclang -fopenmp
INCLUDE=-I/opt/homebrew/opt/libomp/include

all: galaxy.exe clean

galaxy.exe: barnes_hut.o
	gcc -o galaxy.exe barnes_hut.o $(LDFLAGS) $(LIBS)

barnes_hut.o: barnes_hut.c
	gcc -c barnes_hut.c $(CFLAGS) $(INCLUDE)

clean:
	rm -f *.o *.gal
