all: cholesky

cholesky: main.cpp
	g++ -std=c++1y -fopenmp -fpermissive -O3 -Wall main.cpp -lgomp -o cholesky

cholesky_papi: main.cpp
	g++ -std=c++1y -fopenmp -fpermissive -O0 -g3 -Wall main.cpp -lgomp -lpapi -o cholesky

clean:
	rm -v cholesky
