#include <cmath>
#include <iomanip>
#include <iostream>

#include "core/core.h"
#include "linear-algebra/linear-algebra.h"
#include "third-party/third-party.h"

using namespace std;
using namespace __core__;
using namespace __third_party__;

//#include "test/tree-tests.h"
#include "test/cholesky-sparse.h"

int main() {
	cerr<<"************************************************************"<<endl;
	cerr<<"*                  Execution began                         *"<<endl;
	cerr<<"************************************************************"<<endl;
//	sparseCholesky("./test/matrices/1138_bus.rb");
	Matrix<int,cpuAllocator> AD(4,4,-1,0,DEFAULT_STREAM),ARD;
	AD(0,0)=1;
	AD(0,3)=2;
	AD(3,0)=3;
	AD(1,1)=4;
	AD(3,1)=7;
	AD(2,2)=5;
	AD(3,3)=6;
	cerr<<AD<<endl;
	MatrixCXS<int,int,cpuAllocator,CCS> AS,AR;
	convert(AS,AD);
	AS.print(cerr,AS.nzv(),", ","{","}",[](std::ostream &__ost__,auto r,auto c,auto v) -> void { __ost__<<"{"<<r+1<<","<<c+1<<"}->"<<v; })<<endl;
	std::vector<int> permutation({3,0,1,2});
	reorderMatrix(AR,AS,permutation.data());
	convert(ARD,AR);
	cerr<<AR<<endl;
	cerr<<ARD<<endl;
	cerr<<"************************************************************"<<endl;
	cerr<<"*                  Execution ended                         *"<<endl;
	cerr<<"************************************************************"<<endl;
	return 0;
}
