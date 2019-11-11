#ifndef __MEMORY_FREE_H__
#define __MEMORY_FREE_H__

#include <cstdlib>
#include <type_traits>

#include <cuda.h>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>

#include "../macro-definitions.h"
#include "../enum-definitions.h"
#include "../debug/debug.h"
#include "../meta/meta.h"

namespace __core__ {
namespace __memory__ {
//++++++++++++++++++++++++++++++++++++++++++++++++++ default versions ++++++++++++++++++++++++++++++++++++++++++++++++++
//*************************************************** unified memory ***************************************************
template <typename memory_T,typename T,enable_IT<eq_CE(memory_T::location,HOST)> = 0> __forceinline__ void __free__(T* ptr) {
	if(ptr==nullptr)
		return;
	if(!memory_T::is_pinned)
		std::free((void*)ptr);
	else {
		cuda_error(cudaFreeHost(ptr),MEMORY_ERROR);
	}
}
template <typename memory_T,typename T,enable_IT<!eq_CE(memory_T::location,HOST)> = 0> __forceinline__ void __free__(T* ptr) {
	if(ptr==nullptr)
		return;
	cuda_error(cudaFree(ptr),MEMORY_ERROR);
}
//************************************************* not unified memory *************************************************
template <typename memory_T,typename T,enable_IT<eq_CE(memory_T::location,HOST)> = 0> __forceinline__ void __free__(T* ptr,int device) {
	if(ptr==nullptr)
		return;
	if(!memory_T::is_pinned)
		std::free((void*)ptr);
	else {
		if(device>=0) {
			cuda_error(cudaSetDevice(device),API_ERROR);
		}
		cuda_error(cudaFreeHost(ptr),MEMORY_ERROR);
	}
}
template <typename memory_T,typename T,enable_IT<!eq_CE(memory_T::location,HOST)> = 0> __forceinline__ void __free__(T* ptr,int device) {
	if(ptr==nullptr)
		return;
	if(device>=0) {
		cuda_error(cudaSetDevice(device),API_ERROR);
	}
	cuda_error(cudaFree(ptr),MEMORY_ERROR);
}
}
}
#endif
