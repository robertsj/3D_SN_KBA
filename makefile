CC = icpc

ifeq ($(CC), icpc)
	OPENMP = -openmp
else
	OPENMP = -fopenmp
endif

all : 3DKBA

auxiliary_function.o : auxiliary_function.hh auxiliary_function.cc
	$(CC) -c -O3 auxiliary_function.cc
	
miniapp.o : miniapp.hh miniapp.cc
	$(CC) -c -O3 -g $(OPENMP) miniapp.cc
	
main.o : main.cc
	$(CC) -c -O3 main.cc
	
3DKBA : main.o miniapp.o auxiliary_function.o
	$(CC) -O3 $(OPENMP) main.o miniapp.o auxiliary_function.o -o 3DKBA

clean :
	rm *.o 3DKBA