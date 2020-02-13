#ifndef __BOOST_GIO_GRAPH_IO_THIRD_PARTY_H__
#define __BOOST_GIO_GRAPH_IO_THIRD_PARTY_H__

#include <utility>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphml.hpp>

#include "../../wavefront/data-structures/data-structures.h"

namespace __third_party__ {
namespace __graph_io__ {
namespace __boost__ {
template <typename VertexType=void,typename EdgeType=void,typename Allocator=void,typename IT=int>
void write(const __core__::__wavefront__::GraphCXS<VertexType,EdgeType,Allocator,IT>& graph,std::ostream& ost) {
	typedef std::pair<IT,IT> Edge;
    typedef boost::adjacency_list< boost::vecS, boost::vecS, boost::directedS> GraphB;
    std::vector<Edge> edges(graph.e());
    size_t pos=0;
    for(size_t i=0;i<graph.v();++i)
    	for(IT j=graph.ptr(i);j<graph.ptr(i+1);++j)
    		edges[pos++]=Edge(i,graph.indxs(j));
    GraphB g(edges.begin(), edges.end(), graph.v());
    boost::dynamic_properties dp;
    boost::write_graphml(ost, g, dp, true);
}
}
}
}

#endif
