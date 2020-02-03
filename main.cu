#include <cmath>
#include <iomanip>
#include <iostream>

#include "core/core.h"
#include "linear-algebra/linear-algebra.h"
#include "third-party/third-party.h"
#include "wavefront/wavefront.h"

using namespace std;
using namespace __core__;
using namespace __third_party__;


int main(int argc, char **argv) {
	cerr << "************************************************************"<< endl;
	cerr << "*                  Execution began                         *"<< endl;
	cerr << "************************************************************"<< endl;
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
//	gcxs.graph().print(std::cerr)<<std::endl;
	auto fn=[](int i,int tid){
#pragma omp critical
		{
		cerr<<i<<"\t"<<tid<<endl;
		std::flush(std::cerr);
		}
	};
	ompGraph(gcxs,fn,4);
	cerr << "************************************************************"<< endl;
	cerr << "*                  Execution ended                         *"<< endl;
	cerr << "************************************************************"<< endl;
	return 0;
}
