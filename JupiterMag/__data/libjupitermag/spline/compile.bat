
g++ -c -lm -fPIC -std=c++17 -Wextra -g spline.cc -o spline.o
g++ -c -lm -fPIC -std=c++17 -Wextra -g libspline.cc -o libspline.o
	

g++ -lm -fPIC -std=c++17 -Wextra -g -shared -o libspline.so libspline.o spline.o

