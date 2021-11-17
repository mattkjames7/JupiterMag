# CFLAGS for CC
CFLAGS=-g -lm -fopenmp -fPIC -std=c++11 -W

# Compiler
CCo=g++ -c $(CFLAGS)
CC=g++ $(CFLAGS)

all: obj lib test

obj:
	$(CCo) con2020.cc -o con2020.o 
	$(CCo) bessel.cc -o bessel.o 
	$(CCo) trap.cc -o trap.o 

lib:
	$(CC) libcon2020.cc *.o -shared -o libcon2020.so
	
test:
	$(CC) test.cc *.o -o test

clean:
	rm -v *.o 
	rm -v *.so
	rm test
