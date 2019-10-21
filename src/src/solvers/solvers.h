#ifndef __SOLVERS_H__
#define __SOLVERS_H__

#include "../matrix/matrix.h"
#include "cholesky.h"
#include "backward.h"
#include "forward.h"

namespace __core__ {
template <typename T> Matrix<T> cholesky(const Matrix<T>& A) {
	Matrix<T> R=A;
	cholesky(*R,R.n());
	return R;
}
template <typename T> Matrix<T> backward(const Matrix<T>& L,const Matrix<T>& b) {
	Matrix<T> x=b;
    backward(*x,*L,*b,b.n(),b.m());
    return x;
}
template <typename T> Matrix<T> forward(Matrix<T>& U,Matrix<T>& b) {
	Matrix<T> x=b;
    forward(*x,*U,*b,b.n(),b.m());
    return x;
}
}
#endif
