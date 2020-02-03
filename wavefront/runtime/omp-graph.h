#ifndef __OMP_GRAPH_RUNTIME_H__
#define __OMP_GRAPH_RUNTIME_H__

#include "../data-structures/data-structures.h"

namespace __core__ {
namespace __wavefront__ {
namespace __runtime__ {
template <typename FT=void,typename VertexType=void,typename EdgeType=void,typename Allocator=void,typename IT=int,typename...Args>
void ompGraph(GraphCXS<VertexType,EdgeType,Allocator,IT>& graph,FT& function,int threadnum,Args...args) {
	IT* ptr=graph.ptr();
#pragma omp parallel  num_threads(threadnum)
	{
	#pragma omp single
		{
			for(size_t i=0;i<graph.v();++i) {
				IT pos=graph.ptr(i),size=graph.ptr(i+1)-graph.ptr(i);
				if(size>0) {
#pragma omp task depend(iterator(it=0:size), in:ptr[graph.indxs(pos+it)]) depend(out:ptr[i])
					{
						function(i,omp_get_thread_num(),args...);
					}
				}
				else
					break;
			}
		}
	}
}
}
}
}

#endif
