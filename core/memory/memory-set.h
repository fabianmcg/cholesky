#ifndef __MEMORY_SET_H__
#define __MEMORY_SET_H__

#include <type_traits>
#include <future>

#include <cuda.h>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#include <chrono>

#include "../macro-definitions.h"
#include "../enum-definitions.h"
#include "../debug/debug.h"
#include "../meta/meta.h"
#include "util.h"

namespace __core__ {
namespace __memory__ {
//linear memory
template <typename memory_T,SyncBehaviorType sync_behavior=SYNC,typename T=void,enable_IT<(memory_T::location==HOST)&&(sync_behavior==SYNC)> = 0>
void __memset__(T *ptr,int value,std::size_t size,int device=0,cudaStream_t stream=DEFAULT_STREAM) {
	constexpr std::size_t sod=__memory_private__::__sizeof__<T>();
	memset(reinterpret_cast<void*>(ptr),value,size*sod);
}
template <typename memory_T,SyncBehaviorType sync_behavior=SYNC,typename T=void,enable_IT<(memory_T::location==HOST)&&(sync_behavior==ASYNC)> = 0>
std::future<void*> __memset__(T *ptr,int value,std::size_t size,int device=0,cudaStream_t stream=DEFAULT_STREAM) {
	constexpr std::size_t sod=__memory_private__::__sizeof__<T>();
	return std::async(std::launch::async,memset,reinterpret_cast<void*>(ptr),value,size*sod);
}
template <typename memory_T,SyncBehaviorType sync_behavior=ASYNC,typename T=void,enable_IT<memory_T::location!=HOST> = 0>
void __memset__(T *ptr,int value,std::size_t size,int device=-1,cudaStream_t stream=DEFAULT_STREAM) {
	constexpr std::size_t sod=__memory_private__::__sizeof__<T>();
	if(device>=0) {
		cuda_error(cudaSetDevice(device),API_ERROR);
	}
	if(sync_behavior==ASYNC) {
		cuda_error(cudaMemsetAsync(reinterpret_cast<void*>(ptr),value,size*sod,stream),MEMORY_ERROR);
	}
	else {
		cuda_error(cudaMemset(reinterpret_cast<void*>(ptr),value,size*sod),MEMORY_ERROR);
	}
}

//2D memory
template <typename memory_T,SyncBehaviorType sync_behavior=ASYNC,enable_IT<memory_T::location==HOST> = 0>
auto __memset2D__(void *ptr,std::size_t pitch,std::size_t xsize,std::size_t ysize,int value,int device=-1,cudaStream_t stream=DEFAULT_STREAM) {
	return __memset__<memory_T>(ptr,value,pitch*ysize);
}
template <typename memory_T,SyncBehaviorType sync_behavior=ASYNC,enable_IT<memory_T::location!=HOST> = 0>
void __memset2D__(void *ptr,std::size_t pitch,std::size_t xsize,std::size_t ysize,int value,int device=-1,cudaStream_t stream=DEFAULT_STREAM) {
	if(device>=0) {
		cuda_error(cudaSetDevice(device),API_ERROR);
	}
	if(sync_behavior==ASYNC) {
		cuda_error(cudaMemset2DAsync(ptr,pitch,value,xsize,ysize,stream),MEMORY_ERROR);
	}
	else {
		cuda_error(cudaMemset2D(ptr,pitch,value,xsize,ysize),MEMORY_ERROR);
	}
}

//3D memory
template <typename memory_T,SyncBehaviorType sync_behavior=ASYNC,enable_IT<memory_T::location==HOST> = 0>
auto __memset3D__(cudaPitchedPtr &ptr,int value,const cudaExtent &extent,int device=-1,cudaStream_t stream=DEFAULT_STREAM) {
	return __memset__<memory_T>(ptr.ptr,value,ptr.pitch*extent.depth*ptr.ysize);
}
template <typename memory_T,SyncBehaviorType sync_behavior=ASYNC,enable_IT<memory_T::location!=HOST> = 0>
void __memset3D__(cudaPitchedPtr &ptr,int value,const cudaExtent &extent,int device=-1,cudaStream_t stream=DEFAULT_STREAM) {
	if(device>=0) {
		cuda_error(cudaSetDevice(device),API_ERROR);
	}
	if(sync_behavior==ASYNC) {
		cuda_error(cudaMemset3DAsync(ptr,value,extent,stream),MEMORY_ERROR);
	}
	else {
		cuda_error(cudaMemset3D(ptr,value,extent),MEMORY_ERROR);
	}
}
}
}
#endif
