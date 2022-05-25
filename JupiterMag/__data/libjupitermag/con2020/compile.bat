g++ -c -lm -fPIC -std=c++17 -Wextra -g con2020.cc -o con2020.o 
g++ -c -lm -fPIC -std=c++17 -Wextra -g bessel.cc -o bessel.o 
g++ -c -lm -fPIC -std=c++17 -Wextra -g trap.cc -o trap.o 
g++ -c -lm -fPIC -std=c++17 -Wextra -g libcon2020.cc -o libcon2020.o


g++ -lm -fPIC -std=c++17 -Wextra -g *.o -shared -o libcon2020.dll
