#ifndef __CHOLESKY_SPARSE_PARALLEL_TEST_H__
#define __CHOLESKY_SPARSE_PARALLEL_TEST_H__

#include "../core/core.h"
#include "../linear-algebra/linear-algebra.h"
#include "../third-party/third-party.h"

double sparseCholeskyWithReorderingP(const __core__::MatrixCXS<double,int,__core__::cpuAllocator,__core__::CCS>& AF,size_t branchSize=0,double tolerance=0.2,size_t threadnum=2,std::string outName="./cholr") {
 	using namespace std;
	using namespace __core__;
	using namespace __third_party__;
	MatrixCXS<double,int,cpuAllocator,CCS> A,L,AT,AMetis;
	MatrixCXS<int,int,cpuAllocator,CRS> ETree,RP,SETree;
	MatrixCXS<int,int,cpuAllocator,CCS> LP;
	std::vector<int32_t> permutation,invPermutation;
	trimDiagonal(AT,AF);
	int metisCode=reorderPermuation(AT,permutation,invPermutation);
	error(metisCode==METIS_OK,"METIS failed, error code: "+std::to_string(metisCode));
	reorderMatrix(AMetis,AF,permutation.data(),invPermutation.data());
	auto tree=elimitationTree(AMetis);
	auto popermutation=postOrderTree(tree);
	reorderMatrix(A,AMetis,popermutation.data(),popermutation.data()+tree.size());
	tree[0].key=A.n();
	convert(ETree,tree);
	LP=colPattern(A,ETree);
	copyPattern(L,LP);
	RP=std::move(rowPattern(LP));
	MatrixCXSHandler<double,int> AH(A),LH(L);
	MatrixCXSHandler<int,int> LPH(LP),RPH(RP),EH(ETree);
	std::vector<Array<double,cpuAllocator>> c;
	std::vector<ArrayHandler<double,int>> ch;
	c.reserve(threadnum);
	for(size_t i=0;i<threadnum;++i) {
		c.push_back(std::move(Array<double,cpuAllocator>(A.n(),-1,0)));
		ch.push_back(ArrayHandler<double,int>(c[i]));
	}
//	Array<double,cpuAllocator> c(A.n(),-1,0);
//	ArrayHandler<double,int> ch(c);
	Array<int,cpuAllocator> p(A.n(),-1,0);
	ArrayHandler<int,int> ph(p);
	if(branchSize==0)
		branchSize=tree[tree[0].left_child].value/2;
	decltype(tree) STree;
	auto snodes=superNodes(tree,STree,branchSize,tolerance);
//	STree.print(std::cerr)<<std::endl;
	int tmp=STree.size()-1;
	auto potf=[&tmp](Node<int,int,int>& node,long depth){ node.key=tmp; --tmp;};
	STree.traverseLeftRight(potf);
	(*STree).key=STree.size();
	convert(SETree,STree);
//	STree.print(std::cerr)<<std::endl;
	cpu_timer timer;
	timer.start();
	choleskyLeftLookingP(LH,AH,RPH,tree,snodes,ph,ch,threadnum);
	timer.stop();
	cerr<<endl<<"Parallel sparse with reordering:"<<endl<<"\tElapsed time: "<<timer<<endl<<endl;
	c.clear();
	std::ofstream Aout=std::move(open_file<1>(outName+"A.csv"));
	std::ofstream Lout=std::move(open_file<1>(outName+"L.csv"));
	std::ofstream Tout=std::move(open_file<1>(outName+"Tree.csv"));
	std::ofstream SETout=std::move(open_file<1>(outName+"STree.csv"));
	std::ofstream STout=std::move(open_file<1>(outName+"STree.txt"));
	Aout << std::setprecision(12);
	Lout << std::setprecision(12);
	Aout<<A.n()<<","<<A.m()<<endl;
	Lout<<L.n()<<","<<L.m()<<endl;
	Tout<<ETree.n()<<","<<ETree.m()<<endl;
	SETout<<SETree.n()<<","<<SETree.m()<<endl;
	A.print(Aout,A.nzv(),"\n","","",[](std::ostream &__ost__,auto r,auto c,auto v) -> void { __ost__<<""<<r+1<<","<<c+1<<","<<v; })<<endl;
	L.print(Lout,L.nzv(),"\n","","",[](std::ostream &__ost__,auto r,auto c,auto v) -> void { __ost__<<""<<r+1<<","<<c+1<<","<<v; })<<endl;
	ETree.print(Tout,ETree.nzv(),"\n","","",[](std::ostream &__ost__,auto r,auto c,auto v) -> void { __ost__<<""<<r+1<<","<<c+1<<","<<1; })<<endl;
	SETree.print(SETout,SETree.nzv(),"\n","","",[](std::ostream &__ost__,auto r,auto c,auto v) -> void { __ost__<<""<<r+1<<","<<c+1<<","<<1; })<<endl;
	snodes.print(STout)<<endl;
	close_file(Aout);
	close_file(Lout);
	close_file(Tout);
	close_file(STout);
	close_file(SETout);
	return timer.elapsed_time();
}
#endif
