#ifndef __CHOLESKY_SPARSE_TEST_H__
#define __CHOLESKY_SPARSE_TEST_H__

#include "../core/core.h"
#include "../linear-algebra/linear-algebra.h"
#include "../third-party/third-party.h"

int sparseCholesky(std::string filename,std::string outName="./chol") {
	using namespace std;
	using namespace __core__;
	using namespace __third_party__;
	MatrixCXS<double,int,cpuAllocator,CCS> A,L;
	MatrixCXS<int,int,cpuAllocator,CRS> ETree,RP;
	MatrixCXS<int,int,cpuAllocator,CCS> LP;
	A=std::move(read<double,int>(filename));
	auto tree=elimitationTree(A);
	convert(ETree,tree);
	LP=colPattern(A,ETree);
	copyPattern(L,LP);
	RP=std::move(rowPattern(LP));
	MatrixCXSHandler<double,int> AH(A),LH(L);
	MatrixCXSHandler<int,int> LPH(LP),RPH(RP),EH(ETree);
	Array<double,cpuAllocator> c(A.n(),-1,0);
	Array<int,cpuAllocator> p(A.n(),-1,0);
	ArrayHandler<double,int> ch(c);
	ArrayHandler<int,int> ph(p);
	cpu_timer timer;
	timer.start();
	choleskyLeftLooking(LH,AH,RPH,EH,ph,ch);
	timer.stop();
	cerr<<endl<<"\tElapsed time: "<<timer<<endl<<endl;
	std::ofstream Aout=std::move(open_file<1>(outName+"-A.csv"));
	std::ofstream Lout=std::move(open_file<1>(outName+"-L.csv"));
	std::ofstream Tout=std::move(open_file<1>(outName+"-Tree.csv"));
	Aout << std::setprecision(12);
	Lout << std::setprecision(12);
	Aout<<A.n()<<","<<A.m()<<endl;
	Lout<<L.n()<<","<<L.m()<<endl;
	Tout<<ETree.n()<<","<<ETree.m()<<endl;
	A.print(Aout,A.nzv(),"\n","","",[](std::ostream &__ost__,auto r,auto c,auto v) -> void { __ost__<<""<<r+1<<","<<c+1<<","<<v; })<<endl;
	L.print(Lout,L.nzv(),"\n","","",[](std::ostream &__ost__,auto r,auto c,auto v) -> void { __ost__<<""<<r+1<<","<<c+1<<","<<v; })<<endl;
	ETree.print(Tout,ETree.nzv(),"\n","","",[](std::ostream &__ost__,auto r,auto c,auto v) -> void { __ost__<<""<<r+1<<","<<c+1<<","<<1; })<<endl;
	close_file(Aout);
	close_file(Lout);
	close_file(Tout);
	return 0;
}
int sparseCholesky(const __core__::MatrixCXS<double,int,__core__::cpuAllocator,__core__::CCS>& A,std::string outName="./chol") {
	using namespace std;
	using namespace __core__;
	using namespace __third_party__;
	MatrixCXS<double,int,cpuAllocator,CCS> L;
	MatrixCXS<int,int,cpuAllocator,CRS> ETree,RP;
	MatrixCXS<int,int,cpuAllocator,CCS> LP;
	auto tree=elimitationTree(A);
	convert(ETree,tree);
	LP=colPattern(A,ETree);
	copyPattern(L,LP);
	RP=std::move(rowPattern(LP));
	MatrixCXSHandler<double,int> AH(A),LH(L);
	MatrixCXSHandler<int,int> LPH(LP),RPH(RP),EH(ETree);
	Array<double,cpuAllocator> c(A.n(),-1,0);
	Array<int,cpuAllocator> p(A.n(),-1,0);
	ArrayHandler<double,int> ch(c);
	ArrayHandler<int,int> ph(p);
	cpu_timer timer;
	timer.start();
	choleskyLeftLooking(LH,AH,RPH,EH,ph,ch);
	timer.stop();
	cerr<<endl<<"Sparse without reordering:"<<endl<<"\tElapsed time: "<<timer<<endl<<endl;
	std::ofstream Aout=std::move(open_file<1>(outName+"-A.csv"));
	std::ofstream Lout=std::move(open_file<1>(outName+"-L.csv"));
	std::ofstream Tout=std::move(open_file<1>(outName+"-Tree.csv"));
	Aout << std::setprecision(12);
	Lout << std::setprecision(12);
	Aout<<A.n()<<","<<A.m()<<endl;
	Lout<<L.n()<<","<<L.m()<<endl;
	Tout<<ETree.n()<<","<<ETree.m()<<endl;
	A.print(Aout,A.nzv(),"\n","","",[](std::ostream &__ost__,auto r,auto c,auto v) -> void { __ost__<<""<<r+1<<","<<c+1<<","<<v; })<<endl;
	L.print(Lout,L.nzv(),"\n","","",[](std::ostream &__ost__,auto r,auto c,auto v) -> void { __ost__<<""<<r+1<<","<<c+1<<","<<v; })<<endl;
	ETree.print(Tout,ETree.nzv(),"\n","","",[](std::ostream &__ost__,auto r,auto c,auto v) -> void { __ost__<<""<<r+1<<","<<c+1<<","<<1; })<<endl;
	close_file(Aout);
	close_file(Lout);
	close_file(Tout);
	return 0;
}

int sparseCholeskyWithReordering(std::string filename,std::string outName="./cholr") {
	using namespace std;
	using namespace __core__;
	using namespace __third_party__;
	MatrixCXS<double,int,cpuAllocator,CCS> AF,A,L;
	MatrixCXS<int,int,cpuAllocator,CRS> ETree,RP;
	MatrixCXS<int,int,cpuAllocator,CCS> LP;
	AF=std::move(read<double,int>(filename));
	std::vector<int32_t> permutation,invPermutation;
	int metisCode=reorderPermuation(AF,permutation,invPermutation);
	error(metisCode==METIS_OK,"METIS failed, error code: "+std::to_string(metisCode));
	reorderMatrix(A,AF,permutation.data(),invPermutation.data());
	auto tree=elimitationTree(A);
	convert(ETree,tree);
	LP=colPattern(A,ETree);
	copyPattern(L,LP);
	RP=std::move(rowPattern(LP));
	MatrixCXSHandler<double,int> AH(A),LH(L);
	MatrixCXSHandler<int,int> LPH(LP),RPH(RP),EH(ETree);
	Array<double,cpuAllocator> c(A.n(),-1,0);
	Array<int,cpuAllocator> p(A.n(),-1,0);
	ArrayHandler<double,int> ch(c);
	ArrayHandler<int,int> ph(p);
	cpu_timer timer;
	timer.start();
	choleskyLeftLooking(LH,AH,RPH,EH,ph,ch);
	timer.stop();
	cerr<<endl<<"\tElapsed time: "<<timer<<endl<<endl;
	std::ofstream Aout=std::move(open_file<1>(outName+"-A.csv"));
	std::ofstream Lout=std::move(open_file<1>(outName+"-L.csv"));
	std::ofstream Tout=std::move(open_file<1>(outName+"-Tree.csv"));
	Aout << std::setprecision(12);
	Lout << std::setprecision(12);
	Aout<<A.n()<<","<<A.m()<<endl;
	Lout<<L.n()<<","<<L.m()<<endl;
	Tout<<ETree.n()<<","<<ETree.m()<<endl;
	A.print(Aout,A.nzv(),"\n","","",[](std::ostream &__ost__,auto r,auto c,auto v) -> void { __ost__<<""<<r+1<<","<<c+1<<","<<v; })<<endl;
	L.print(Lout,L.nzv(),"\n","","",[](std::ostream &__ost__,auto r,auto c,auto v) -> void { __ost__<<""<<r+1<<","<<c+1<<","<<v; })<<endl;
	ETree.print(Tout,ETree.nzv(),"\n","","",[](std::ostream &__ost__,auto r,auto c,auto v) -> void { __ost__<<""<<r+1<<","<<c+1<<","<<1; })<<endl;
	close_file(Aout);
	close_file(Lout);
	close_file(Tout);
	return 0;
}
int sparseCholeskyWithReordering(const __core__::MatrixCXS<double,int,__core__::cpuAllocator,__core__::CCS>& AF,std::string outName="./cholr") {
 	using namespace std;
	using namespace __core__;
	using namespace __third_party__;
	MatrixCXS<double,int,cpuAllocator,CCS> A,L,AT,AMetis;
	MatrixCXS<int,int,cpuAllocator,CRS> ETree,RP;
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
	Array<double,cpuAllocator> c(A.n(),-1,0);
	Array<int,cpuAllocator> p(A.n(),-1,0);
	ArrayHandler<double,int> ch(c);
	ArrayHandler<int,int> ph(p);
	cpu_timer timer;
	timer.start();
	choleskyLeftLooking(LH,AH,RPH,EH,ph,ch);
	timer.stop();
	cerr<<endl<<"Sparse with reordering:"<<endl<<"\tElapsed time: "<<timer<<endl<<endl;
	std::ofstream Aout=std::move(open_file<1>(outName+"-A.csv"));
	std::ofstream Lout=std::move(open_file<1>(outName+"-L.csv"));
	std::ofstream Tout=std::move(open_file<1>(outName+"-Tree.csv"));
	Aout << std::setprecision(12);
	Lout << std::setprecision(12);
	Aout<<A.n()<<","<<A.m()<<endl;
	Lout<<L.n()<<","<<L.m()<<endl;
	Tout<<ETree.n()+1<<","<<ETree.m()+1<<endl;
	A.print(Aout,A.nzv(),"\n","","",[](std::ostream &__ost__,auto r,auto c,auto v) -> void { __ost__<<""<<r+1<<","<<c+1<<","<<v; })<<endl;
	L.print(Lout,L.nzv(),"\n","","",[](std::ostream &__ost__,auto r,auto c,auto v) -> void { __ost__<<""<<r+1<<","<<c+1<<","<<v; })<<endl;
	ETree.print(Tout,ETree.nzv(),"\n","","",[](std::ostream &__ost__,auto r,auto c,auto v) -> void { __ost__<<""<<r+1<<","<<c+1<<","<<1; })<<endl;
	close_file(Aout);
	close_file(Lout);
	close_file(Tout);
	return 0;
}
#endif
