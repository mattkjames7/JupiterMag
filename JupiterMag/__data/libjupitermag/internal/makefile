# CFLAGS for CC
CFLAGS=-lm -fopenmp -fPIC -std=c++11 

# Compiler
CCo=g++ -c $(CFLAGS)
CC=g++ $(CFLAGS)

all: obj lib

obj:
	#this bit converts the internale field model coefficients into 
	#something addressable within the library
	python3 EncodeCoefficients.py
	ld -r -b binary -o jrm09coeffs.o jrm09coeffs.bin
	ld -r -b binary -o vip4coeffs.o vip4coeffs.bin
	$(CCo) internal.cc -o internal.o 
	$(CCo) libinternal.cc -o libinternal.o

lib:
	$(CC) *.o -shared -o libinternal.so

clean:
	rm -v *.o *.so *.bin
	
