# Makefile for simulation

LDFLAGS=-lm
CFLAGS = -g -Wall -Wno-deprecated -O3 -march=native -std=c11

all: galaxy.exe clean

galaxy.exe: barnes_hut.o
	gcc -o galaxy.exe barnes_hut.o $(LDFLAGS)

barnes_hut.o: barnes_hut.c
	gcc -c barnes_hut.c $(CFLAGS)

clean:
	rm -f *.o *.gal
