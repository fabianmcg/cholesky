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

#endif
