# CFLAGS for CC
CFLAGS=-lm -fPIC -std=c++17 -Wextra

# Compiler
CCo=g++ -c $(CFLAGS)
CC=g++ $(CFLAGS)
CCWo=x86_64-w64-mingw32-g++-win32 -c $(CFLAGS)
CCW=x86_64-w64-mingw32-g++-win32 $(CFLAGS)

all: obj lib 

obj:
	$(CCo) con2020.cc -o con2020.o 
	$(CCo) bessel.cc -o bessel.o 
	$(CCo) polyeval.cc -o polyeval.o 
	$(CCo) trap.cc -o trap.o 
	$(CCo) smoothd.cc -o smoothd.o 
	$(CCo) libcon2020.cc -o libcon2020.o

lib:
	$(CC) *.o -shared -o libcon2020.so
	
windows:
	$(CCWo) con2020.cc -o ../build/windows/con2020.o 
	$(CCWo) bessel.cc -o ../build/windows/bessel.o 
	$(CCWo) polyeval.cc -o ../build/windows/polyeval.o 
	$(CCWo) trap.cc -o ../build/windows/trap.o 
	$(CCWo) smoothd.cc -o ../build/windows/smoothd.o 
	$(CCWo) libcon2020.cc -o ../build/windows/libcon2020.o

test:
	$(CC) test.cc *.o -o test

clean:
	rm -v *.o 
	rm -v *.so
	#rm test
