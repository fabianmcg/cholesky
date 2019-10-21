#ifndef __SOLVERS_BACKWARD_H__
#define __SOLVERS_BACKWARD_H__

#include <cmath>

namespace __core__ {
template <typename T> void backward(T* x,const T* L,const T* b,size_t n,size_t m) {
#pragma omp parallel for num_threads(4)
    for(size_t k=0;k<m;++k) {
        for(size_t i=0;i<n;++i) {
            x[i*m+k]=b[i*m+k];
            for(size_t j=0;j<i;++j)
                x[i*m+k]-=L[i*n+j]*x[j*m+k];
            x[i*m+k]=x[i*m+k]/L[i*n+i];
        }
    }
}
}
#endif
