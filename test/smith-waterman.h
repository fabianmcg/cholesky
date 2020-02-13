#ifndef SMITH_WATERMAN_H_
#define SMITH_WATERMAN_H_



//	auto G=smithWatermanDG(3,3);
//	GraphCXS<int,double,cpuAllocator> gcxs;
//	convert(gcxs,G);
//	gcxs.graph().print(std::cerr,gcxs.graph().nzv(),", ","{","}",[](std::ostream &__ost__,auto r,auto c,auto v) -> void { __ost__<<"{"<<r+1<<","<<c+1<<"}->"<<1; })<<endl;

	//	std::ofstream Gout=std::move(open_file<1>("g.graphml"));
	//	write(gcxs,Gout);
	//	close_file(Gout);
	//	auto SG=smithWatermanSDG(8,8,2,2);
	//	SG.graph().print(cerr,SG.graph().nzv(),", ","{","}",[](std::ostream &__ost__,auto r,auto c,auto v) -> void { __ost__<<"{"<<r+1<<","<<c+1<<"}->"<<1; })<<endl;
	//	std::vector<int> order,indegree;
	//	order.reserve(SG.v());
	//	auto push=[&order](int v,decltype(SG)& g){ order.push_back(v); };
	//	topologicalSort(SG,push,indegree);
	//	order.pop_back();
	//	cerr<<order<<endl;
	//	auto fn=[&SG](int v,int tid){
	//#pragma omp critical
	//		cerr<<v<<" "<<tid<<endl;
	//	};
	//	ompGraph<OMPUserOrder>(SG,fn,order,2);
#endif
