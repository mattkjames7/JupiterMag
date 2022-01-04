# CFLAGS for CC
CFLAGS=-lm -fopenmp -fPIC -std=c++17 -g

# Compiler
CCo=g++ -c $(CFLAGS)
CC=g++ $(CFLAGS)

all: obj lib

obj:
	#this bit converts the internale field model coefficients into 
	#something addressable within the library
	python3 EncodeCoefficients.py
	ld -r -b binary -o jrm09coeffs.o jrm09coeffs.bin
	ld -r -b binary -o jrm33coeffs.o jrm33coeffs.bin
	ld -r -b binary -o vip4coeffs.o vip4coeffs.bin
	ld -r -b binary -o gsfc13evcoeffs.o gsfc13evcoeffs.bin
	ld -r -b binary -o gsfc15evcoeffs.o gsfc15evcoeffs.bin
	ld -r -b binary -o gsfc15evscoeffs.o gsfc15evscoeffs.bin
	ld -r -b binary -o isaaccoeffs.o isaaccoeffs.bin
	ld -r -b binary -o jpl15evcoeffs.o jpl15evcoeffs.bin
	ld -r -b binary -o jpl15evscoeffs.o jpl15evscoeffs.bin
	ld -r -b binary -o o4coeffs.o o4coeffs.bin
	ld -r -b binary -o o6coeffs.o o6coeffs.bin
	ld -r -b binary -o p11acoeffs.o p11acoeffs.bin
	ld -r -b binary -o shacoeffs.o shacoeffs.bin
	ld -r -b binary -o u17evcoeffs.o u17evcoeffs.bin
	ld -r -b binary -o v117evcoeffs.o v117evcoeffs.bin
	ld -r -b binary -o vipalcoeffs.o vipalcoeffs.bin
	ld -r -b binary -o vit4coeffs.o vit4coeffs.bin
	$(CCo) internal.cc -o internal.o 
	$(CCo) internalmodel.cc -o internalmodel.o 
	$(CCo) libinternal.cc -o libinternal.o

lib:
	$(CC) *.o -shared -o libinternal.so

clean:
	rm -v *.o *.so *.bin
	
