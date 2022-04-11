# CFLAGS for CC
CFLAGS=-lm -fPIC -std=c++17 -Wextra

# Compiler
CCo=g++ -c $(CFLAGS)
CC=g++ $(CFLAGS)

all: obj lib 

obj:
	$(CCo) con2020.cc -o con2020.o 
	$(CCo) bessel.cc -o bessel.o 
	$(CCo) trap.cc -o trap.o 
	$(CCo) libcon2020.cc -o libcon2020.o

lib:
	$(CC) *.o -shared -o libcon2020.so
	
test:
	$(CC) test.cc *.o -o test

clean:
	rm -v *.o 
	rm -v *.so
	#rm test
