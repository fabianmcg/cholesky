#ifndef __CUDA_UTIL_H__
#define __CUDA_UTIL_H__

#include <cuda.h>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#include <driver_types.h>

#include "../enum-definitions.h"
#include "../debug/debug.h"

namespace __core__ {
namespace __util__ {
namespace __cuda__ {
__forceinline__  int device_count() {
	int i=0;
	cuda_error(cudaGetDeviceCount(&i),API_ERROR);
	return i;
}
inline int valid_device(const int i) {
	if(i<0)
		return 0;
	else if(i<device_count())
		return i;
	return 0;
}
inline int get_device() {
	int i=0;
	cuda_error(cudaGetDevice(&i),API_ERROR);
	return i;
}
__forceinline__ bool visible_devices(int first_device,int second_device) {
	if(first_device==second_device)
		return true;
	int is_visible=0;
	cuda_error(cudaDeviceCanAccessPeer(&is_visible,first_device,second_device),API_ERROR);
	return is_visible==1;
}
__forceinline__ cudaPointerAttributes get_ptr_attributes(void *ptr) {
	cudaPointerAttributes R;
	cuda_error(cudaPointerGetAttributes(&R,ptr),API_ERROR);
	return R;
}
int get_ptr_dev(void *ptr) {
	cudaPointerAttributes R;
	cuda_error(cudaPointerGetAttributes(&R,ptr),API_ERROR);
	return R.device;
}
bool is_ptr_at_dev(void *ptr,int dev) {
	cudaPointerAttributes R;
	cuda_error(cudaPointerGetAttributes(&R,ptr),API_ERROR);
	return dev==R.device;
}
}
}
}
#endif
