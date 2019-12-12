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
//#include "test/tree-tests.h"

int main(int argc, char **argv) {
	cerr << "************************************************************"
			<< endl;
	cerr << "*                  Execution began                         *"
			<< endl;
	cerr << "************************************************************"
			<< endl;
//	std::string filename="./test/matrices/1138_bus.rb";
	std::string filename = "./test/matrices/LF1"
			"0.rb";
	if (argc > 1)
		filename = std::string(argv[1]);
	MatrixCXS<double, int, cpuAllocator, CCS> A, AT;
	A = std::move(read<double, int>(filename));
//	choleskyDense(A);
	sparseCholesky(A);
	sparseCholeskyWithReordering(A);
	cerr << "************************************************************"
			<< endl;
	cerr << "*                  Execution ended                         *"
			<< endl;
	cerr << "************************************************************"
			<< endl;
	return 0;
}
