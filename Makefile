
all: cholesky

cholesky: main.cu
	nvcc  -O3 -Xcompiler -std=c++14 -Xcompiler -fopenmp -Xcompiler -fpermissive -Xcompiler -Wall -Xptxas -v -gencode arch=compute_50,code=sm_50 main.cu -lgomp -lrbio -lcholmod -lmetis -o cholesky

gnu:
	cp -v main.cu main.cpp
	g++-9 -O3 -std=c++14 -fopenmp -w -D_OPENMP=201811 main.cpp -lgomp -lmetis -lrbio -lcholmod -o main
	rm -v main.cpp
	
clean:
	rm -v cholesky main
