#include <cmath>
#include <iomanip>
#include <iostream>

#include "core/core.h"
#include "linear-algebra/linear-algebra.h"
#include "third-party/third-party.h"

using namespace std;
using namespace __core__;
using namespace __third_party__;

//#include "test/cholesky-dense-timing.h"
#include "test/cholesky-sparse.h"
#include "test/cholesky-parallel.h"
#include "test/cholmod.h"
//#include "test/tree-tests.h"

int main(int argc, char **argv) {
	cerr << "************************************************************"<< endl;
	cerr << "*                  Execution began                         *"<< endl;
	cerr << "************************************************************"<< endl;
//	std::string filename = "./test/matrices/LFAT5.rb";
//	std::string filename = "./test/matrices/LF10.rb";
//	std::string filename="./test/matrices/ex5.rb";
//	std::string filename="./test/matrices/bcsstk01.rb";
//	std::string filename="./test/matrices/bcsstk02.rb";
	std::string filename="./test/matrices/1138_bus.rb";
//	std::string filename="./test/matrices/bcsstk16.rb";
//	std::string filename="./test/matrices/mhd4800b.rb";
//	std::string filename="./test/matrices/torsion1.rb";
//	std::string filename="./test/matrices/cvxbqp1.rb";
	size_t branchSize=200,threadnum=2;
	double tolerance=0.75;
	if(argc > 1)
		filename = std::string(argv[1]);
	if(argc > 3) {
		threadnum=stol(argv[2]);
		branchSize=stol(argv[3]);
	}
	if(argc == 5)
		tolerance=stod(argv[4]);
	MatrixCXS<double, int, cpuAllocator, CCS> A, AT;
	A = std::move(read<double, int>(filename));
	cerr << "File name:\n\t"<<filename<<endl;
	cerr << "Threadnum:\n\t"<<threadnum<<endl;
	cerr << "Branch size:\n\t"<<branchSize<<endl;
	cerr << "Matrix dimensions:\n\t"<<A.n()<<"x"<<A.m()<<"\tnzv: "<<A.nzv()<<endl;
	cerr << "************************************************************"<< endl;
	auto t1=sparseCholeskyWithReordering(A,"./");
	auto t2=cholmod(A,true);
	auto t3=cholmod(A,false);
//	auto t2=sparseCholeskyWithReorderingTT(A,"./TT-");
	auto t4=sparseCholeskyWithReorderingP(A,branchSize,tolerance,threadnum,"./P-");
	cerr<<"Serial-CHOLMOD speedup:\t\t"<<t2/t1<<endl;
	cerr<<"Supernodal-simplicial speedup:\t"<<t3/t2<<endl;
	cerr<<"PTT-Serial speedup:\t\t"<<t1/t4<<endl;
	cerr<<"CHOL-PTT speedup:\t\t"<<__min__(t2,t3)/t4<<endl;
	cerr << "************************************************************"<< endl;
	cerr << "*                  Execution ended                         *"<< endl;
	cerr << "************************************************************"<< endl;
	return 0;
}
