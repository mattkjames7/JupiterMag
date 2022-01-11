# CFLAGS for CC
CFLAGS=-lm -fopenmp -g -std=c++17

# Compiler
CCo=g++ -fPIC -c $(CFLAGS)
CC=g++ -fPIC $(CFLAGS)


all: obj lib test

obj:
	cd con2020;	make all
	cd spline; make all

	
	$(CCo) model.cc -o model.o
	$(CCo) trace.cc -o trace.o
	$(CCo) interptraceclosestpos.cc -o interptraceclosestpos.o

	
lib:
	$(CC) -shared libjupitermag.cc *.o \
			con2020/*.o \
			spline/*.o \
			-o libjupitermag.so -Wl,-rpath='$$ORIGIN/libinternalfield' -L ./libinternalfield -linternalfield

test:
	$(CC) test.cc *.o \
			con2020/*.o \
			spline/*.o \
			-o test	-Wl,-rpath='$$ORIGIN/libinternalfield' -L ./libinternalfield -linternalfield


clean:
	cd con2020; make clean
	cd spline; make clean
	rm *.o
	rm -v  *.so
	rm -v test
