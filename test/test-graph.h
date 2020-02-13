#ifndef TEST_GRAPH_H_
#define TEST_GRAPH_H_

#include <cmath>
#include <iomanip>
#include <iostream>

#include "../core/core.h"
#include "../linear-algebra/linear-algebra.h"
#include "../third-party/third-party.h"
#include "../wavefront/wavefront.h"

int test_graph(int argc, char **argv) {
	using namespace std;
	using namespace __core__;
	Graph<cpuAllocator,void> G(4,3,0);
	G.initVertices();
	G.initEdges();
	G.print();
	cerr<<"\t"<<G.edge_quantity()<<"\t"<<G.vertex_quantity()<<"\t"<<G.vertices().capacity()<<"\t"<<G.edges().capacity()<<endl;
	G.insertEdge(G[0],G[1]);
	G.insertEdge(G[1],G[2]);
	G.insertEdge(G[2],G[3]);
	G.print();
	cerr<<"\t"<<G.edge_quantity()<<"\t"<<G.vertex_quantity()<<"\t"<<G.vertices().capacity()<<"\t"<<G.edges().capacity()<<endl;
	G.insertEdge(G[1],G[0]);
	G.insertEdge(G[3],G[3]);
	G.print();
	cerr<<"\t"<<G.edge_quantity()<<"\t"<<G.vertex_quantity()<<"\t"<<G.vertices().capacity()<<"\t"<<G.edges().capacity()<<endl;
	G.eraseEdge(G[1],G[0]);
	G.eraseVertex(G[3]);
	G.print();
	cerr<<"\t"<<G.edge_quantity()<<"\t"<<G.vertex_quantity()<<"\t"<<G.vertices().capacity()<<"\t"<<G.edges().capacity()<<endl;
	return 0;
}
int test_graph2(){
	using namespace std;
	using namespace __core__;
	using namespace __third_party__;
	Graph<cpuAllocator,void> G(9,12,0),GC;
	G.initVertices();
	G.initEdges();
	G.insertEdge(G[0],G[1]);
	G.insertEdge(G[0],G[2]);
	G.insertEdge(G[0],G[3]);
	G.insertEdge(G[0],G[4]);
	G.insertEdge(G[1],G[5]);
	G.insertEdge(G[2],G[5]);
	G.insertEdge(G[3],G[6]);
	G.insertEdge(G[4],G[6]);
	G.insertEdge(G[5],G[7]);
	G.insertEdge(G[6],G[7]);
	G.insertEdge(G[7],G[8]);
	G.print();
	GraphCXS<int,double,cpuAllocator> gcxs;
	convert(gcxs,G);
	gcxs.graph().print(std::cerr)<<std::endl;
	vector<int> ind;
	auto pr=[](int v,GraphCXS<int,double,cpuAllocator>& graph) {
		cerr<<v<<endl;
	};
	topologicalSort(gcxs,pr,ind);
	std::ofstream Gout=std::move(open_file<1>("g.graphml"));
	write(gcxs,Gout);
	close_file(Gout);
	return 0;
}
#endif
