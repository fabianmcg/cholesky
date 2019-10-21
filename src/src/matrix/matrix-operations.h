#ifndef __MOP_H__
#define __MOP_H__

#include <chrono>
#include <random>

#include "matrix.h"

namespace __core__ {
template <typename T> Matrix<T> identity(size_t n) {
	Matrix<T> I(n);
    for(size_t i=0;i<n;++i)
        I(i,i)=1;
    return I;
}
template <typename T> Matrix<T> randomTriangularMatrix(size_t n,size_t m,T a=0.,T b=1.) {
	size_t seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator(seed);
	std::uniform_real_distribution<double>distribution(-1.0,0.0);
	Matrix<T> A(n,m);
	for(size_t i=0;i<n;++i)
		for(size_t j=0;j<m;++j)
			if(i<=j)
				A(i,j)=(-distribution(generator))*(b-a)+a;
	return A;
}
template <typename T> Matrix<T> randomTriangularMatrix(size_t n,T a=0.,T b=1.) {
    return randomTriangularMatrix<T>(n,n,a,b);
}
template <typename T> Matrix<T> randomPDMatrix(size_t n,T a=0.,T b=1.) {
	Matrix<T> B=randomTriangularMatrix<T>(n,a,b);
    return B*B.transpose();
}
}
#endif
 
