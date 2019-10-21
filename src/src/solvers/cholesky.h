#ifndef __SOLVERS_CHOLESKY_H__
#define __SOLVERS_CHOLESKY_H__

#include <cerrno>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <iostream>

#include <papi.h>

namespace __core__ {

//template <typename T> void cholesky(T* R,size_t n) {
//    for(size_t k=0;k<n;++k) {
//        for(size_t j=k+1;j<n;++j)
//            for(size_t i=j;i<n;++i)
//                R[j*n+i]-=R[k*n+i]*R[k*n+j]/R[k*n+k];
//        T tmp=sqrt(std::abs(R[k*n+k]));
//        for(size_t j=k;j<n;++j)
//            R[k*n+j]=R[k*n+j]/tmp;
//    }
//    for(size_t i=0;i<n;++i)
//        for(size_t j=0;j<n;++j)
//            if(i>j)
//                R[i*n+j]=0;
//}

static void __check_papi__(std::string file,int line,std::string call,int retval) {
	std::cerr<<file<<"\tFAILED\nLine #"<<line<<std::endl;
    if ( retval == PAPI_ESYS )
        std::cerr<<"System error in: "<<call<<std::endl;
    else if ( retval > 0 )
    	std::cerr<<"Error calculating: "<<call<<std::endl;
    else
    	std::cerr<<"Error in "<< call<<": "<<PAPI_strerror(retval)<<std::endl;
    std::cerr<<std::endl;
    exit(1);
}

template <typename T> void cholesky(T* R,size_t n) {
//	constexpr int Events[]={PAPI_DP_OPS};
//	constexpr int Events[]={PAPI_L1_LDM,PAPI_L2_LDM,PAPI_L3_LDM,PAPI_LD_INS};
//	constexpr int Events[]={PAPI_L1_TCM,PAPI_L2_TCM,PAPI_L3_TCM,PAPI_TOT_INS};
	constexpr int Events[]={PAPI_LD_INS,PAPI_SR_INS,PAPI_TOT_INS};
	constexpr int NUM_EVENTS=sizeof(Events)/sizeof(int);
	long_long values[NUM_EVENTS];
	int retval;
	if((retval=PAPI_start_counters(const_cast<int*>(Events),NUM_EVENTS))!=PAPI_OK)
		__check_papi__(__FILE__, __LINE__, "PAPI_start_counters", retval);
	for(size_t k=0;k<n;++k) {
		for(size_t j=k+1;j<n;++j)
			for(size_t i=j;i<n;++i)
				R[j*n+i]-=R[k*n+i]*R[k*n+j]/R[k*n+k];
		T tmp=sqrt(std::abs(R[k*n+k]));
		for(size_t j=k;j<n;++j)
			R[k*n+j]=R[k*n+j]/tmp;
	}

	if((retval=PAPI_read_counters(values,NUM_EVENTS))!=PAPI_OK)
		__check_papi__(__FILE__, __LINE__, "PAPI_read_counters", retval);

	for(size_t i=0;i<n;++i)
		for(size_t j=0;j<n;++j)
			if(i>j)
				R[i*n+j]=0;
	for(size_t i=0;i<NUM_EVENTS;++i) {
		PAPI_event_info_t info;
		PAPI_get_event_info(Events[i],&info);
		std::cerr<<"Event[ "<<info.short_descr<<" ]:\t"<<values[i]<<std::endl;
	}
	std::cerr<<std::endl;
    PAPI_shutdown();
}

//template <typename T> void cholesky(T* R,size_t n) {
//	int Events[]={PAPI_DP_OPS,PAPI_L1_LDM,PAPI_L2_LDM,PAPI_L3_LDM,PAPI_L1_TCM,PAPI_L2_TCM,PAPI_L3_TCM,PAPI_LD_INS};
//	float real_time, proc_time, mflops;
//	long long flpins;
//	int retval;
//	if((retval=PAPI_ipc( &real_time, &proc_time, &flpins, &mflops))<PAPI_OK)
//		__check_papi__(__FILE__, __LINE__, "PAPI_flops", retval);
//
//	for(size_t k=0;k<n;++k) {
//		for(size_t j=k+1;j<n;++j)
//			for(size_t i=j;i<n;++i)
//				R[j*n+i]-=R[k*n+i]*R[k*n+j]/R[k*n+k];
//		T tmp=sqrt(std::abs(R[k*n+k]));
//		for(size_t j=k;j<n;++j)
//			R[k*n+j]=R[k*n+j]/tmp;
//	}
//
//	if((retval=PAPI_ipc( &real_time, &proc_time, &flpins, &mflops))<PAPI_OK)
//		__check_papi__(__FILE__, __LINE__, "PAPI_flops", retval);
//
//	for(size_t i=0;i<n;++i)
//		for(size_t j=0;j<n;++j)
//			if(i>j)
//				R[i*n+j]=0;
//	printf("Real_time:\t%f\nProc_time:\t%f\nTotal flpins:\t%lld\nMFLOPS:\t\t%f\n",
//    real_time, proc_time, flpins, mflops);
//    printf("%s\tPASSED\n", __FILE__);
//    PAPI_shutdown();
//}
}
#endif
