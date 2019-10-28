#ifndef __SOLVERS_CHOLESKY_H__
#define __SOLVERS_CHOLESKY_H__

#include <cerrno>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <iostream>
#include <type_traits>

#include <papi.h>

namespace __core__ {

template <int ver,typename T=void> typename std::enable_if<ver==0,void>::type cholesky(T* R,size_t n) {
    for(size_t k=0;k<n;++k) {
        for(size_t j=k+1;j<n;++j)
            for(size_t i=j;i<n;++i)
                R[j*n+i]-=R[k*n+i]*R[k*n+j]/R[k*n+k];
        T tmp=sqrt(std::abs(R[k*n+k]));
        for(size_t j=k;j<n;++j)
            R[k*n+j]=R[k*n+j]/tmp;
    }
    for(size_t i=0;i<n;++i)
        for(size_t j=0;j<n;++j)
            if(i>j)
                R[i*n+j]=0;
}
template <int ver,typename T=void> typename std::enable_if<ver==1,void>::type cholesky(T* R,size_t n) {
    for(size_t k=0;k<n;++k) {
#pragma omp parallel for num_threads(4)
        for(size_t j=k+1;j<n;++j)
            for(size_t i=j;i<n;++i)
                R[j*n+i]-=R[k*n+i]*R[k*n+j]/R[k*n+k];
        T tmp=sqrt(std::abs(R[k*n+k]));
#pragma omp parallel for num_threads(4)
        for(size_t j=k;j<n;++j)
            R[k*n+j]=R[k*n+j]/tmp;
    }
    for(size_t i=0;i<n;++i)
        for(size_t j=0;j<n;++j)
            if(i>j)
                R[i*n+j]=0;
}
template <int ver,typename T=void> typename std::enable_if<ver==2,void>::type  cholesky(T* R,size_t n) {
#pragma omp parallel num_threads(4)
	{
#pragma omp single
		{
			for(size_t k=0;k<n;++k) {
				for(size_t j=k+1;j<n;++j)
					for(size_t i=j;i<n;++i) {
#pragma omp task depend(in:R[k*n+i]) depend(in:R[k*n+j]) depend(in:R[k*n+k]) depend(out: R[j*n+i])
						R[j*n+i]-=R[k*n+i]*R[k*n+j]/R[k*n+k];
					}
				for(size_t j=k+1;j<n;++j)
#pragma omp task depend(out:R[k*n+j]) depend(in:R[k*n+k])
					R[k*n+j]=R[k*n+j]/sqrt(std::abs(R[k*n+k]));
#pragma omp task depend(out:R[k*n+k])
				R[k*n+k]=R[k*n+k]/sqrt(std::abs(R[k*n+k]));
			}
		}
	}
    for(size_t i=0;i<n;++i)
        for(size_t j=0;j<n;++j)
            if(i>j)
                R[i*n+j]=0;
}
}
#endif
