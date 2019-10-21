#ifndef __SOLVERS_FORWARD_H__
#define __SOLVERS_FORWARD_H__

#include <cmath>

namespace __core__ {
template <typename T> void forward(T* x,const T* U,const T* b,size_t n,size_t m) {
#pragma omp parallel for num_threads(4)
    for(long k=0;k<m;++k) {
        for(long i=n-1;i>=0;i=i-1) {
            x[i*m+k]=b[i*m+k];
            for(long j=n-1;j>i;--j)
                x[i*m+k]-=U[i*n+j]*x[j*m+k];
            x[i*m+k]=x[i*m+k]/U[i*n+i];
        }
    }
}
}
#endif
