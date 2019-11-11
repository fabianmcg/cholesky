all: cholesky

cholesky: main.cu
	nvcc  -O3 -Xcompiler -std=c++14 -Xcompiler -fopenmp -Xcompiler -fpermissive -Xcompiler -Wall -Xptxas -v -gencode arch=compute_50,code=sm_50 main.cu -lgomp -lrbio -o cholesky

clean:
	rm -v cholesky
