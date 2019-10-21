#ifndef __TIMING_H__
#define __TIMING_H__

#include <chrono>
#include <cmath>
#include <iostream>

namespace __core__ {
struct cpu_timer {
	std::chrono::high_resolution_clock::time_point t1;
	std::chrono::high_resolution_clock::time_point t2;
	__inline__ __attribute__((always_inline)) void start() {
		t1=std::chrono::high_resolution_clock::now();
	}
	__inline__ __attribute__((always_inline)) void stop() {
		t2=std::chrono::high_resolution_clock::now();
	}
	__inline__ __attribute__((always_inline)) double elapsed_time()const {
		return static_cast<double>(std::chrono::duration_cast<std::chrono::nanoseconds>( t2 - t1 ).count())*pow(10.,-9.);
	}

};
std::ostream &operator<<(std::ostream &oss,const cpu_timer &timer) {
	oss<<timer.elapsed_time();
	return oss;
}
template <typename T> __inline__ __attribute__((always_inline)) T& operator<<(T &t,const cpu_timer &timer) {
	t=timer.elapsed_time();
	return t;
}
}
#endif
 
