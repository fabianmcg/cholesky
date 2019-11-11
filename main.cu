#include <cmath>
#include <iomanip>
#include <iostream>

#include "core/core.h"
#include "linear-algebra/linear-algebra.h"
#include "third-party/third-party.h"

using namespace std;
using namespace __core__;
using namespace __third_party__;

#define VAL M_PI

int main() {
	Matrix<double,cpuAllocator> A,R,U,AC;
	MatrixCXS<double,int,cpuAllocator,CRS> aCRS;
	string filename="./test/t2dal_e.rb";
	aCRS=std::move(read<double,int>(filename));
	convert(A,aCRS);
	R=A;
	MatrixHandler<double,size_t> RH(R);
	cpu_timer timer;
	timer.start();
	cholesky<1>(RH);
	timer.stop();
	transpose(U,R);
	multiply(AC,U,R);
	cerr<<"********************************************************************************"<<endl;
	cerr<<endl<<endl<<"\tMatrix filename: "<<filename<<"\n\t\tRows: "<<aCRS.rows()<<"\n\t\tCols: "<<aCRS.cols()<<"\n\t\tNZV: "<<aCRS.nzv()<<endl<<endl;
	cerr<<"********************************************************************************"<<endl;
	cerr<<"\tDense:"<<endl;
	cerr<<"\t\tError: "<<maxDifference(A,AC)<<endl;
	cerr<<"\t\tTiming: "<<timer<<endl;
	cerr<<"********************************************************************************"<<endl;
	MatrixMap<double> aMap,rMap,uMap,acMap;
	convert(aMap,aCRS);
	rMap=choleskyInput(aMap);
	timer.start();
	cholesky<0>(rMap);
	timer.stop();
	transpose(uMap,rMap);
	acMap=multiply(uMap,rMap);
	cerr<<"\tSparse:"<<endl;
	cerr<<"\t\tError: "<<maxDifference(aMap,acMap)<<endl;
	cerr<<"\t\tTiming: "<<timer<<endl;
	cerr<<"********************************************************************************"<<endl;
	cerr<<"\tmax(abs(Dense-Sparse)):"<<endl;
	Matrix<double,cpuAllocator> acD;
	convert(acD,acMap);
	cerr<<"\t\tError: "<<maxDifference(A,acD)<<endl;
	cerr<<"********************************************************************************"<<endl;
	return 0;
}
