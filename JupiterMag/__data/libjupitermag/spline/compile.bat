
g++ -c -lm -fPIC -std=c++17 -Wextra spline.cc -o spline.o
g++ -c -lm -fPIC -std=c++17 -Wextra libspline.cc -o libspline.o
	

g++ -lm -fPIC -std=c++17 -Wextra -shared -o libspline.so libspline.o spline.o

