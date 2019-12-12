#ifndef __TEST_PERMUTE_H__
#define __TEST_PERMUTE_H__


#include "../core/core.h"
#include "../linear-algebra/linear-algebra.h"
#include "../third-party/third-party.h"

int permuteTest() {
	using namespace std;
	using namespace __core__;
	using namespace __third_party__;
	Matrix<int,cpuAllocator> AD(4,4,-1,0,DEFAULT_STREAM),ARD;
	AD(0,0)=1;
	AD(0,3)=2;
	AD(3,0)=3;
	AD(1,1)=4;
	AD(3,1)=7;
	AD(2,2)=5;
	AD(3,3)=6;
	cerr<<AD<<endl<<endl;
	MatrixCXS<int,int,cpuAllocator,CCS> AS,AR;
	convert(AS,AD);
	std::vector<int> permutation({3, 0, 1, 2});
	reorderMatrix(AR,AS,permutation.data());
	convert(ARD,AR);
	print(AS.ptr(),AS.n()+1)<<endl;
	print(AS.indxs(),AS.nzv())<<endl;
	print(AR.ptr(),AR.n()+1)<<endl;
	print(AR.indxs(),AR.nzv())<<endl<<endl;
	AS.print(cerr,AS.nzv(),", ","{","}",[](std::ostream &__ost__,auto r,auto c,auto v) -> void { __ost__<<"{"<<r+1<<","<<c+1<<"}->"<<v; })<<endl;
	AR.print(cerr,AR.nzv(),", ","{","}",[](std::ostream &__ost__,auto r,auto c,auto v) -> void { __ost__<<"{"<<r+1<<","<<c+1<<"}->"<<v; })<<endl<<endl;
	cerr<<ARD<<endl;
	return 0;
}

int permuteTest(std::string filename) {
	using namespace std;
	using namespace __core__;
	using namespace __third_party__;
	MatrixCXS<double,int,cpuAllocator,CCS> A;
	A=std::move(read<double,int>(filename));
	std::vector<int32_t> p,ip;
	int metisCode=reorderPermuation(A,p,ip);
	error(metisCode==METIS_OK,"METIS failed, error code: "+std::to_string(metisCode));
	cerr<<std::fixed<<std::setprecision(8);
	A.print(cerr,A.nzv(),", ","{","}",[](std::ostream &__ost__,auto r,auto c,auto v) -> void { __ost__<<"{"<<r+1<<","<<c+1<<"}->"<<v; })<<endl;
	cerr<<p<<endl;
	return 0;
}
#endif
