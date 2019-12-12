#ifndef __TREE_TESTS_H__
#define __TREE_TESTS_H__

#include "../core/core.h"
#include "../linear-algebra/linear-algebra.h"
#include "../third-party/third-party.h"

int buildTree() {
	using namespace std;
	using namespace __core__;
	using namespace __third_party__;
	LRTree<int,void,cpuAllocator> tree;
	auto tmp=tree.insert(1,LRTree<int,void,cpuAllocator>::nilNode);
	tree.insert(4,tree[tmp.self]);
	tree.insert(5,tree[tmp.self]);
	tree.insert(2,LRTree<int,void,cpuAllocator>::nilNode);
	tmp=tree.insert(3,LRTree<int,void,cpuAllocator>::nilNode);
	tmp=tree.insert(6,tree[tmp.self]);
	auto tmp2=tree.insert(7,tree[tmp.self]);
	tree.insert(8,tree[tmp.self]);
	tree.insert(9,tree[tmp2.self]);
	tree.print(cerr)<<endl;
	tree.moveUnder(tree.find(7).self,tree.find(1).self);
	tree.moveUnder(tmp.self,tree.find(9).self);
	tree.print(cerr)<<endl;
	cerr<<tree.find(6)<<endl;
	tree.erase(tree.find(6));
	tree.print(cerr);
	tree.erase(tree.find(3));
	tree.print(cerr);
	cerr<<tree.find(3)<<endl;
	return 0;
}
int eliminationTreeReport() {
	using namespace std;
	using namespace __core__;
	using namespace __third_party__;
	Matrix<int,cpuAllocator> A(11,11,-1,0,DEFAULT_STREAM);
	for(size_t i=0;i<11;++i)
		A(i,i)=1;
	A(2,1)=1;
	A(5,0)=1;
	A(5,3)=1;
	A(6,0)=1;
	A(7,1)=1;
	A(7,4)=1;
	A(8,5)=1;
	A(9,2)=1;
	A(9,3)=1;
	A(9,5)=1;
	A(9,7)=1;
	A(10,2)=1;
	A(10,4)=1;
	A(10,6)=1;
	A(10,7)=1;
	A(10,9)=1;
	for(size_t i=0;i<11;++i)
		for(size_t j=i+1;j<11;++j)
			A(i,j)=A(j,i);
	MatrixCXS<int,int,cpuAllocator,CRS> Acxs,ETree;
	MatrixCXS<int,int,cpuAllocator,CCS> L,AP,Acrs;
	convert(Acxs,A);
	convert(Acrs,A);
//	cerr<<Acxs<<endl;
	auto tree=elimitationTree(Acxs);
	tree.print(cerr);
	auto p=postOrderTree(tree);
	tree.print(cerr);
	convert(ETree,tree);
	cerr<<p<<endl;
	reorderMatrix(AP,Acrs,p.data(),p.data()+A.n());
	Acrs.print(cerr,Acrs.nzv(),", ","{","}",[](std::ostream &__ost__,auto r,auto c,auto v) -> void { __ost__<<"{"<<r+1<<","<<c+1<<"}->"<<v; })<<endl;
	AP.print(cerr,AP.nzv(),", ","{","}",[](std::ostream &__ost__,auto r,auto c,auto v) -> void { __ost__<<"{"<<r+1<<","<<c+1<<"}->"<<v; })<<endl;
//	cerr<<A<<endl<<endl;
//	convert(A,AP);
//	cerr<<A<<endl<<endl;
//	ETree.print(cerr,ETree.nzv(),", ","{","}",[](std::ostream &__ost__,auto r,auto c,auto v) -> void { __ost__<<"{"<<r+1<<","<<c+1<<"}->"<<v; })<<endl;
//	auto LP=colPattern(Acxs,ETree);
//	LP.print(cerr,LP.nzv(),", ","{","}",[](std::ostream &__ost__,auto r,auto c,auto v) -> void { __ost__<<"{"<<r+1<<","<<c+1<<"}->"<<v; })<<endl;
//	copyPattern(L,LP);
//	L.print(cerr,L.nzv(),", ","{","}",[](std::ostream &__ost__,auto r,auto c,auto v) -> void { __ost__<<"{"<<r+1<<","<<c+1<<"}->"<<v; })<<endl;
//	auto RP=rowPattern(LP);
//	RP.print(cerr,RP.nzv(),", ","{","}",[](std::ostream &__ost__,auto r,auto c,auto v) -> void { __ost__<<"{"<<r+1<<","<<c+1<<"}->"<<v; })<<endl;
	return 0;
}
int eliminationTree(std::string filename,std::ostream& ost=std::cerr) {
	using namespace std;
	using namespace __core__;
	using namespace __third_party__;
	MatrixCXS<double,int,cpuAllocator,CCS> A,L;
	MatrixCXS<int,int,cpuAllocator,CRS> ETree;
	MatrixCXS<int,int,cpuAllocator,CCS> LP;
	A=std::move(read<double,int>(filename));
	auto tree=elimitationTree(A);
	convert(ETree,tree);
	A.print(ost,A.nzv(),", ","{","}",[](std::ostream &__ost__,auto r,auto c,auto v) -> void { __ost__<<"{"<<r+1<<","<<c+1<<"}->"<<v; })<<endl;
	ETree.print(ost,ETree.nzv(),", ","{","}",[](std::ostream &__ost__,auto r,auto c,auto v) -> void { __ost__<<"{"<<r+1<<","<<c+1<<"}->"<<v; })<<endl;
	LP=colPattern(A,ETree);
	copyPattern(L,LP);
	L.print(ost,L.nzv(),", ","{","}",[](std::ostream &__ost__,auto r,auto c,auto v) -> void { __ost__<<"{"<<r+1<<","<<c+1<<"}->"<<v+1; })<<endl;
	auto RP=rowPattern(LP);
	RP.print(ost,RP.nzv(),", ","{","}",[](std::ostream &__ost__,auto r,auto c,auto v) -> void { __ost__<<"{"<<r+1<<","<<c+1<<"}->"<<v; })<<endl;
	return 0;
}
#endif
