# CFLAGS for CC
CFLAGS=-lm -fopenmp -g -std=c++11 

# Compiler
CCo=g++ -fPIC -c $(CFLAGS)
CC=g++ -fPIC $(CFLAGS)


all: obj lib 

obj:
	cd con2020;	make all
	cd internal; make all

	
	$(CCo) model.cc -o model.o

	
lib:
	$(CC) -shared libjupitermag.cc *.o \
			internal/*.o \
			con2020/*.o \
			-o libjupitermag.so


clean:
	cd internal; make clean
	cd con2020; make clean
	rm *.o
	rm -v  *.so
