cd libinternalfield
call compile.bat
cd ..

cd con2020
call compile.bat
cd ..

cd spline 
call compile.bat
cd ..



g++ -fPIC -c -lm -std=c++17 -Wextra model.cc -o model.o
g++ -fPIC -c -lm -std=c++17 -Wextra trace.cc -o trace.o
g++ -fPIC -c -lm -std=c++17 -Wextra interptraceclosestpos.cc -o interptraceclosestpos.o

	

g++ -fPIC -lm -std=c++17 -Wextra -shared libjupitermag.cc *.o con2020/*.o spline/*.o libinternalfield/libinternalfield/*.o -o libjupitermag.dll
