# CFLAGS for CC
CFLAGS=-lm -std=c++17 -Wextra -g

# Compiler
CCo=g++ -fPIC -c $(CFLAGS)
CC=g++ -fPIC $(CFLAGS)
CCWo=x86_64-w64-mingw32-g++-win32 -c $(CFLAGS)
CCW=x86_64-w64-mingw32-g++-win32 $(CFLAGS)

all: obj lib

obj:
	cd con2020;	make all
	cd spline; make all
	cd libinternalfield; make all
	
	$(CCo) model.cc -o model.o
	$(CCo) trace.cc -o trace.o
	$(CCo) interptraceclosestpos.cc -o interptraceclosestpos.o

	
lib:
	$(CC) -shared libjupitermag.cc *.o \
			con2020/*.o \
			spline/*.o \
			libinternalfield/libinternalfield.a \
			-o libjupitermag.so 
			#-o libjupitermag.so -Wl,-rpath='$$ORIGIN/libinternalfield' -L ./libinternalfield -linternalfield

test:
	$(CC) test.cc *.o \
			con2020/*.o \
			spline/*.o \
			-o test	-Wl,-rpath='$$ORIGIN/libinternalfield' -L ./libinternalfield -linternalfield


windows:
	mkdir -pv build/windows
	-rm -v build/windows/*.o
	cd con2020;	make windows
	cd spline; make windows
	cd libinternalfield; make winobj
	cp -v libinternalfield/libinternalfield/build/windows/*.o build/windows/
	$(CCWo) model.cc -o build/windows/model.o
	$(CCWo) trace.cc -o build/windows/trace.o
	$(CCWo) interptraceclosestpos.cc -o build/windows/interptraceclosestpos.o	
	$(CCW) -shared libjupitermag.cc build/windows/*.o -o libjupitermag.dll


clean:
	cd con2020; make clean
	cd spline; make clean
	cd libinternalfield; make clean
	rm -f *.o
	rm -vf *.so
	rm -vf test
