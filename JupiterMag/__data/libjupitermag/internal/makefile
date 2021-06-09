# CFLAGS for CC
CFLAGS=-lm -fopenmp

# Compiler
CCo=g++ -fPIC -c $(CFLAGS)
CC=g++ $(CFLAGS)

all: obj lib

obj:
	#this bit converts the internale field model coefficients into 
	#something addressable within the library
	ld -r -b binary -o jrm09coeffs.o jrm09_coeffs.bin
	ld -r -b binary -o vip4coeffs.o vip4_coeffs.bin

lib:
