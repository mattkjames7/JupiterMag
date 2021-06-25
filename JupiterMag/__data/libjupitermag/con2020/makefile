# CFLAGS for CC
CFLAGS=-lm -fopenmp -fPIC

# Compiler
CCo=g++ -c $(CFLAGS)
CC=g++ $(CFLAGS)

all: obj lib

obj:
	$(CCo) con2020.cc -o con2020.o 
	$(CCo) bessel.cc -o bessel.o 
	$(CCo) trap.cc -o trap.o 

lib:
	$(CC) libcon2020.cc *o -shared -o libcon2020.so

clean:
	rm -v *.o *.so *.bin
	
