#ifndef __CHOLESKY_DENSE_TIMING_H__
#define __CHOLESKY_DENSE_TIMING_H__

#include "../core/core.h"
#include "../linear-algebra/linear-algebra.h"
#include "../third-party/third-party.h"

int choleskyDenseTest(std::string filename) {
	using namespace std;
	using namespace __core__;
	using namespace __third_party__;
	Matrix<double,cpuAllocator> A,R,U,AC;
	MatrixCXS<double,int,cpuAllocator,CCS> aCCS;
	MatrixCXS<double,int,cpuAllocator,CRS> aCRS;
	aCCS=std::move(read<double,int>(filename));
	aCCS.values(1)=1;
	convert(aCRS,aCCS);
	aCRS.values(1)=2;
	convert(aCCS,aCRS);
	convert(A,aCCS);
	cerr<<std::scientific;
	cerr<<A<<endl;
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

int choleskyDense(const __core__::MatrixCXS<double,int,__core__::cpuAllocator,__core__::CCS>& A) {
	using namespace std;
	using namespace __core__;
	using namespace __third_party__;
	Matrix<double,cpuAllocator> R;
	convert(R,A);
	MatrixHandler<double,size_t> RH(R);
	cpu_timer timer;
	timer.start();
	cholesky<1>(RH);
	timer.stop();
	cerr<<endl<<"Dense:"<<endl;
	cerr<<"\tElapsed time: "<<timer<<endl<<endl;
	return 0;
}

#endif
