all : 3DKBA

auxiliary_function.o : auxiliary_function.hh auxiliary_function.cc
	g++ -c -O3 auxiliary_function.cc
	
miniapp.o : miniapp.hh miniapp.cc
	g++ -c -O3 -fopenmp miniapp.cc
	
main.o : main.cc
	g++ -c -O3 main.cc
	
3DKBA : main.o miniapp.o auxiliary_function.o
	g++ -O3 -fopenmp -o 3DKBA main.o miniapp.o auxiliary_function.o

clean :
	rm *.o 3DKBA
