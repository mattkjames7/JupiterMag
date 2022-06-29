# CFLAGS for CC
CFLAGS=-lm -std=c++17 -Wextra -g

# Compiler
CCo=g++ -fPIC -c $(CFLAGS)
CC=g++ -fPIC $(CFLAGS)
CCWo=x86_64-w64-mingw32-g++-win32 -c $(CFLAGS)
CCW=x86_64-w64-mingw32-g++-win32 $(CFLAGS)

all: obj lib

obj:
	$(CCo) spline.cc -o spline.o
	$(CCo) libspline.cc -o libspline.o

windows:
	$(CCWo) spline.cc -o ../build/windows/spline.o
	$(CCWo) libspline.cc -o ../build/windows/libspline.o
		
lib: 
	$(CC) -shared -o libspline.so libspline.o spline.o
	
clean:
	rm -v *.o *.so
