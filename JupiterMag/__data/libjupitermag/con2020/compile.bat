g++ -c -lm -fPIC -std=c++17 -Wextra con2020.cc -o con2020.o 
g++ -c -lm -fPIC -std=c++17 -Wextra bessel.cc -o bessel.o 
g++ -c -lm -fPIC -std=c++17 -Wextra trap.cc -o trap.o 
g++ -c -lm -fPIC -std=c++17 -Wextra polyeval.cc -o polyeval.o 
g++ -c -lm -fPIC -std=c++17 -Wextra libcon2020.cc -o libcon2020.o


g++ -lm -fPIC -std=c++17 -Wextra *.o -shared -o libcon2020.dll
