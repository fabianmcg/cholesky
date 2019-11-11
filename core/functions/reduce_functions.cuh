#ifndef __REDUCE_FUNCTIONS_CUH__
#define __REDUCE_FUNCTIONS_CUH__

#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#include <cooperative_groups.h>

#include "../meta/meta.h"

namespace __core__
{
//template <typename T,uint A> Vector<T,2,A> __shfl_down_sync(unsigned mask,Vector<T,2,A>& v,int lane)
//{
//	Vector<T,2,A> result;
//	result[0]=__shfl_down_sync(mask,v[0],lane);
//	result[1]=__shfl_down_sync(mask,v[1],lane);
//	return result;
//}
//template <typename T,uint A> Vector<T,3,A> __shfl_down_sync(unsigned mask,Vector<T,3,A>& v,int lane)
//{
//	Vector<T,3,A> result;
//	result[0]=__shfl_down_sync(mask,v[0],lane);
//	result[1]=__shfl_down_sync(mask,v[1],lane);
//	result[2]=__shfl_down_sync(mask,v[2],lane);
//	return result;
//}
template <typename fn_T,bool unrolled,typename T,enable_IT<unrolled==false> = 0> __forceinline__ __forceflatten__ __optimize__ __device__ T reduce_warp(T data)
{
	for(int i=16;i>0;i=i>>1)
		data=fn_T::fn(__shfl_down_sync(0xFFFFFFFF,data,i),data);
	return data;
}
template <typename fn_T,bool unrolled,typename T,enable_IT<unrolled==true> = 0> __forceinline__ __forceflatten__ __optimize__ __device__ T reduce_warp(T data)
{
	data=fn_T::fn(__shfl_down_sync(0xFFFFFFFF,data,16),data);
	data=fn_T::fn(__shfl_down_sync(0xFFFFFFFF,data,8),data);
	data=fn_T::fn(__shfl_down_sync(0xFFFFFFFF,data,4),data);
	data=fn_T::fn(__shfl_down_sync(0xFFFFFFFF,data,2),data);
	data=fn_T::fn(__shfl_down_sync(0xFFFFFFFF,data,1),data);
	return data;
}
template <typename fn_T,typename T> __forceinline__ __forceflatten__ __optimize__ __device__ T reduce_active(cooperative_groups::coalesced_group thr_group,T data)
{
	T tmp;
	for(unsigned int i=1;i<thr_group.size();i=i<<1)
	{
		tmp=fn_T::fn(thr_group.shfl_down(data,i),data);
		if(i+thr_group.thread_rank()<thr_group.size())
			data=tmp;
	}
	return data;
}
template <bool unrolled,typename fn_T,typename T,enable_IT<unrolled==false> = 0> __forceinline__ __forceflatten__ __optimize__ __device__ T reduce_warp(T data,fn_T& fn)
{
	for(int i=16;i>0;i=i>>1)
		data=fn.reduce(__shfl_down_sync(0xFFFFFFFF,data,i),data);
	return data;
}
template <bool unrolled,typename fn_T,typename T,enable_IT<unrolled==true> = 0> __forceinline__ __forceflatten__ __optimize__ __device__ T reduce_warp(T data,fn_T& fn)
{
	data=fn.reduce(__shfl_down_sync(0xFFFFFFFF,data,16),data);
	data=fn.reduce(__shfl_down_sync(0xFFFFFFFF,data,8),data);
	data=fn.reduce(__shfl_down_sync(0xFFFFFFFF,data,4),data);
	data=fn.reduce(__shfl_down_sync(0xFFFFFFFF,data,2),data);
	data=fn.reduce(__shfl_down_sync(0xFFFFFFFF,data,1),data);
	return data;
}
template <typename fn_T,typename T> __forceinline__ __forceflatten__ __optimize__ __device__ T reduce_active(cooperative_groups::coalesced_group thr_group,T data,fn_T& fn)
{
	T tmp;
	for(unsigned int i=1;i<thr_group.size();i=i<<1)
	{
		tmp=fn.reduce(thr_group.shfl_down(data,i),data);
		if(i+thr_group.thread_rank()<thr_group.size())
			data=tmp;
	}
	return data;
}
template <typename reduce_FT,typename apply_FT,bool unrolled,typename T,typename...Args> __forceflatten__ __optimize__ __device__ inline
enable_T<(unrolled==true)&&(sizeof...(Args)>0),T> reduce_warp(T result,Args...args)
{
	result=apply_FT::fn(result,args...);
	for(int i=16;i>0;i=i>>1)
		result=reduce_FT::fn(__shfl_down_sync(0xFFFFFFFF,result,i),result);
	return result;
}
template <typename reduce_FT,typename apply_FT,bool unrolled,typename T,typename...Args> __forceflatten__ __optimize__ __device__ inline
enable_T<(unrolled==false)&&(sizeof...(Args)>0),T> reduce_warp(T result,Args...args)
{
	result=apply_FT::fn(result,args...);
	result=reduce_FT::fn(__shfl_down_sync(0xFFFFFFFF,result,16),result);
	result=reduce_FT::fn(__shfl_down_sync(0xFFFFFFFF,result,8),result);
	result=reduce_FT::fn(__shfl_down_sync(0xFFFFFFFF,result,4),result);
	result=reduce_FT::fn(__shfl_down_sync(0xFFFFFFFF,result,2),result);
	result=reduce_FT::fn(__shfl_down_sync(0xFFFFFFFF,result,1),result);
	return result;
}
template <typename reduce_FT,typename apply_FT,typename T,typename...Args> __forceflatten__ __optimize__ __device__
enable_T<(sizeof...(Args)>0),T> reduce_active(cooperative_groups::coalesced_group thr_group,T result,Args...args)
{
	T tmp;
	result=apply_FT::fn(result,args...);
	for(unsigned int i=1;i<thr_group.size();i=i<<1)
	{
		tmp=reduce_FT::fn(thr_group.shfl_down(result,i),result);
		if(i+thr_group.thread_rank()<thr_group.size())
			result=tmp;
	}
	return result;
}
template <typename fn_T,int BLOCKDIM=128,bool unrolled=true,int BDIML2=log_2_CE<BLOCKDIM>(),int WCOUNT=(BLOCKDIM>>5),typename T=void> __forceinline__ __forceflatten__ __optimize__ __device__
T reduce_block(T data)
{
	__shared__ T tmp_arr[WCOUNT];
	T result=reduce_warp<fn_T,unrolled>(data);
	if(threadIdx.x&31==0)
		tmp_arr[threadIdx.x>>5]=result;
	if(threadIdx.x==0)
		for(int i=1;i<WCOUNT;++i)
			result=fn_T::fn(tmp_arr[i],result);
	return result;
}

template <typename reduce_FT,typename T,typename IV,bool unrolled=true> struct reduce_A3
{
	typedef Vector<T,3> MDT;
	static constexpr MDT init_value={IV::value,IV::value,IV::value};
	template <typename AT,typename int_T,uchar atype> static __attr_opfiffhd__ MDT read(__Array__<AT,int_T,atype> &arr,int_T i)
	{
		MDT r;
		r[0]=arr(i,0);
		r[1]=arr(i,1);
		r[2]=arr(i,2);
		return r;
	}
	template <typename AT,typename int_T,uchar atype> static __attr_opfiffhd__ void write(__Array__<AT,int_T,atype> &arr,int_T i,const MDT &val)
	{
		arr(i,0)=val[0];
		arr(i,1)=val[1];
		arr(i,2)=val[2];
	}
	static __attr_opfiffhd__ MDT reduce(const MDT& v1,const MDT& v2)
	{
		MDT r;
		r[0]=reduce_FT::fn(v1[0],v2[0]);
		r[1]=reduce_FT::fn(v1[1],v2[1]);
		r[2]=reduce_FT::fn(v1[2],v2[2]);
		return r;
	}
	static __attr_opfiffhd__ MDT reduce_warp(const MDT& v)
	{
		MDT r;
		r[0]=__core__::reduce_warp<reduce_FT,unrolled>(v(0));
		r[1]=__core__::reduce_warp<reduce_FT,unrolled>(v(1));
		r[2]=__core__::reduce_warp<reduce_FT,unrolled>(v(2));
		return r;
	}
};

template <typename reduce_FT,typename RT,typename int_T,typename RIV,typename IIV,bool unrolled=true> struct __mXx_index__
{
	typedef cpair<RT,int_T> MDT;
	static constexpr MDT init_value={RIV::value,IIV::value};
	template <typename int_AT,uchar atype> static __attr_opfiffhd__ MDT read(__Array__<RT,int_AT,atype> &arr,int_AT i)
	{
		MDT pair;
		pair.x=arr[i];
		pair.y=i;
		return pair;
	}
	template <typename int_AT,uchar atype> static __attr_opfiffhd__ MDT read(__Array__<MDT,int_AT,atype> &arr,int_AT i)
	{
		MDT pair;
		pair=arr[i];
		return pair;
	}
	template <typename int_AT,uchar atype> static __attr_opfiffhd__ void write(__Array__<MDT,int_AT,atype> &arr,int_AT i,const MDT &val)
	{
		arr[i]=val;
	}
	static __attr_opfiffhd__ MDT reduce(const MDT& v1,const MDT& v2)
	{
		MDT r;
		r.x=reduce_FT::fn(v1.x,v2.x);
		r.y=(v1.x==r.x)?v1.y:v2.y;
		return r;
	}
	template <bool ur=unrolled,enable_IT<ur==false> = 0> static __attr_opfiff__ __device__ MDT reduce_warp(const MDT& v)
	{
		MDT r=v;
		for(int i=16;i>0;i=i>>1)
		{
			RT tmp=r.x;
			int_T sdv=__shfl_down_sync(0xFFFFFFFF,r.y,i);
			r.x=reduce_FT::fn(__shfl_down_sync(0xFFFFFFFF,r.x,i),r.x);
			r.y=(tmp==r.x)?r.y:sdv;
		}
		return r;
	}
	template <bool ur=unrolled,enable_IT<ur==true> = 0> static __attr_opfiff__ __device__ MDT reduce_warp(const MDT& v)
	{
		MDT r=v;
		RT tmp=r.x;
		int_T sdv=__shfl_down_sync(0xFFFFFFFF,r.y,16);
		r.x=reduce_FT::fn(__shfl_down_sync(0xFFFFFFFF,r.x,16),r.x);
		r.y=(tmp==r.x)?r.y:sdv;
		tmp=r.x;
		sdv=__shfl_down_sync(0xFFFFFFFF,r.y,8);
		r.x=reduce_FT::fn(__shfl_down_sync(0xFFFFFFFF,r.x,8),r.x);
		r.y=(tmp==r.x)?r.y:sdv;
		tmp=r.x;
		sdv=__shfl_down_sync(0xFFFFFFFF,r.y,4);
		r.x=reduce_FT::fn(__shfl_down_sync(0xFFFFFFFF,r.x,4),r.x);
		r.y=(tmp==r.x)?r.y:sdv;
		tmp=r.x;
		sdv=__shfl_down_sync(0xFFFFFFFF,r.y,2);
		r.x=reduce_FT::fn(__shfl_down_sync(0xFFFFFFFF,r.x,2),r.x);
		r.y=(tmp==r.x)?r.y:sdv;
		tmp=r.x;
		sdv=__shfl_down_sync(0xFFFFFFFF,r.y,1);
		r.x=reduce_FT::fn(__shfl_down_sync(0xFFFFFFFF,r.x,1),r.x);
		r.y=(tmp==r.x)?r.y:sdv;
		return r;
	}
};
}
#endif
