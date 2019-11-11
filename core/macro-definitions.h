#ifndef __MACRO_DEFINITIONS_H__
#define __MACRO_DEFINITIONS_H__

//functional macros
#define STRIP_PARENTHESIS_H(...) __VA_ARGS__
#define INVOKE(expr) expr
#define STRIP_PARENTHESIS(args) INVOKE(STRIP_PARENTHESIS_H args)

#define CHOOSE_MACRO(_0,_1,_2,_3,_4,_5,_6,_7, name, ...) name
#define BAD_ARG_COUNT(...) static_assert(false,"Bad quantity of arguments in macro, line: "+to_string(__LINE__)+", in file: "+__FILE__)
#define JOIN(_0,_1) _0##_1

#define UFOR_WA(counter_type,counter_name,begin,end,arguments_types,arguments_types_and_names,arguments,function_body) { typedef struct { \
	static inline __attribute__((always_inline)) __attribute__((optimize(3))) __host__ __device__ void fn(counter_type counter_name,STRIP_PARENTHESIS(arguments_types_and_names)){\
	STRIP_PARENTHESIS(function_body)} \
	} unrolled_fn;\
	const_for<unrolled_fn,counter_type,begin,end-1>::template iterator<STRIP_PARENTHESIS(arguments_types)>(STRIP_PARENTHESIS(arguments)); }

#define UFOR_NA(counter_type,counter_name,begin,end,function_body) { typedef struct { \
	static inline __attribute__((always_inline)) __attribute__((optimize(3))) __host__ __device__ void fn(counter_type counter_name){\
	STRIP_PARENTHESIS(function_body)} \
	} unrolled_fn;\
	const_for<unrolled_fn,counter_type,begin,end-1>::iterator(); }

#define UFOR(...) CHOOSE_MACRO(__VA_ARGS__,UFOR_WA,BAD_ARG_COUNT,BAD_ARG_COUNT,UFOR_NA,BAD_ARG_COUNT,BAD_ARG_COUNT,BAD_ARG_COUNT)(__VA_ARGS__)

#define __macro_for_WA__(counter_type,counter_name,begin,end,arguments_types,arguments_types_and_names,arguments,function_body) \
	for(counter_type counter_name=(counter_type)begin;counter_name<(counter_type)end;++counter_name) { \
	STRIP_PARENTHESIS(function_body) }
#define __macro_for_NA__(counter_type,counter_name,begin,end,function_body) \
	for(counter_type counter_name=(counter_type)begin;counter_name<(counter_type)end;++counter_name) { \
	STRIP_PARENTHESIS(function_body) }

#define MACRO_FOR(...) CHOOSE_MACRO(__VA_ARGS__,__macro_for_WA__,BAD_ARG_COUNT,BAD_ARG_COUNT,__macro_for_NA__,BAD_ARG_COUNT,BAD_ARG_COUNT,BAD_ARG_COUNT)(__VA_ARGS__)

//compiler attributes
#ifndef __forceinline__
#define __forceinline__ __inline__ __attribute__((always_inline))
#endif

#ifndef _inline_
#define __inline0__
#define __inline1__ inline
#define __inline2__ __forceinline__
#define _inline_helper_(level) __inline##level##__
#define _inline_(level) _inline_helper_(level)
#endif

#ifndef __optimize__
#define __optimize__  __attribute__((optimize(3)))
#endif

#ifndef ___optimize___
#define ___optimize___(x)  __attribute__((optimize(x)))
#endif

#ifndef __forceflatten__
#define __forceflatten__ __attribute__((flatten))
#endif

#ifndef __flatten__
#define __noflatten__
#define __flatten_helper__(x) __##x##flatten__
#define __flatten__(x) __flatten_helper__(x)
#endif

#define __host_device__ __host__ __device__

#define __attr_opfiffhd__ __optimize__ __forceinline__ __forceflatten__ __host_device__
#define __attr_opfiff__ __optimize__ __forceinline__ __forceflatten__
#define __attr_opfihd__ __optimize__ __forceinline__ __host_device__
#define __attr_opffhd__ __optimize__ __forceflatten__ __host_device__
#define __attr_fiffhd__ __forceinline__ __forceflatten__ __host_device__
#define __attr_opfi__ __optimize__ __forceinline__
#define __attr_opff__ __optimize__ __forceflatten__
#define __attr_ophd__ __optimize__ __host_device__
#define __attr_fiff__ __forceinline__ __forceflatten__
#define __attr_fihd__ __forceinline__ __host_device__
#define __attr_ffhd__ __forceflatten__ __host_device__

#define __unroll_gcc__ __attribute__((optimize("unroll-loops")))
#ifdef __CUDA_ARCH__
#define __unroll_cuda__ #pragma unroll
#else
#define __unroll_cuda__
#endif
#define __unroll_meta__(...) UFOR(__VA_ARGS__)

#define META_UNROLL
#ifdef META_UNROLL
#define __unroll_gpu__(...) UFOR(__VA_ARGS__)
#define __unroll_cpu__(...) UFOR(__VA_ARGS__)
#else
#ifdef PARTIAL_META_UNROLL
#define __unroll_gpu__(...) MACRO_FOR(__VA_ARGS__)
#define __unroll_cpu__(...) UFOR(__VA_ARGS__)
#else
#define __unroll_gpu__(...)	MACRO_FOR(__VA_ARGS__)
#define __unroll_cpu__(...) MACRO_FOR(__VA_ARGS__)
#endif
#endif

#ifdef __CUDACC_VER_MAJOR__
typedef cudaStream_t StreamType;
#else
typedef int StreamType;
#endif
static constexpr StreamType DEFAULT_STREAM=0;

//type definitions
#define RESTRICT_Q(x) x __restrict__
#define CRESTRICT_Q(x) const x __restrict__

#define C_A(type,x) constant_argument<type,__constant_argument__<type,x>>::value
#define CTA(type,x) constant_argument<type,__constant_argument__<type,template_argument<type>(x)>>
#endif
