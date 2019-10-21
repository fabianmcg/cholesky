#include <iomanip>
#include <iostream>

#include "src/core.h"

using namespace std;
using namespace __core__;

typedef double real_T;
int main(int argc,char *argv[]) {
	bool print=false;
	bool random=true;
    cerr<<std::fixed;
    cerr<<std::setprecision(9);
    double ct=0,bt=0,ft=0;
    size_t n=40;
    size_t nr=1;
    Matrix<real_T> A;
    if(argc>1) {
    	real_T a=0,b=1;
        n=stoi(argv[1]);
        if(argc>=3)
        	nr=stoi(argv[2]);
        if(argc>=4)
        	print=(stoi(argv[3])==1);
        if(argc==6) {
        	a=stod(argv[4]);
        	b=stod(argv[5]);
        }
        A=randomPDMatrix<real_T>(n,a,b);
    	cerr<<"************************************************************************************************"<<endl<<endl;
    	cerr<<"Test parameters:\n\tMatrix size:\t"<<n<<"x"<<n<<"\n\tNumber of repetitions:\t"<<nr<<"\n\tPrint output:\t"<<print
    			<<"\n\tInterval (a,b]:\t("<<a<<","<<b<<"]"<<endl<<endl;
    	cerr<<"************************************************************************************************"<<endl<<endl;
    }
    else if (random) {
    	real_T a=0,b=1;
        A=randomPDMatrix<real_T>(n,a,b);
    	cerr<<"************************************************************************************************"<<endl<<endl;
    	cerr<<"Test parameters:\n\tMatrix size:\t"<<n<<"x"<<n<<"\n\tNumber of repetitions:\t"<<nr<<"\n\tPrint output:\t"<<print
    			<<"\n\tInterval (a,b]:\t("<<a<<","<<b<<"]"<<endl<<endl;
    	cerr<<"************************************************************************************************"<<endl<<endl;
    }
    else {
    	cerr<<"************************************************************************************************"<<endl<<endl;
    	cerr<<"Usage: cholesky matrixSize(I) repetitions(I) printMatrix(I) a(d) b(d)"<<endl<<endl;
    	cerr<<"************************************************************************************************"<<endl<<endl;
    	A={1,2,3,0,4,5,0,0,6};
    	A.reshape(n);
        A=A*A.transpose();
    }
    A.simetrize();
    if(print) cerr<<A<<endl<<endl;
    Matrix<real_T> U(A);
    cpu_timer timer;
    for(size_t k=0;k<nr;++k) {
    	timer.start();
    	U=cholesky(A);
    	timer.stop();
    	ct+=timer.elapsed_time();
    }
    if(print) cerr<<U<<endl<<endl;
    Matrix<real_T> I=identity<real_T>(n);
    Matrix<real_T> L=U.transpose();
    Matrix<real_T> LI;
    for(size_t k=0;k<nr;++k) {
    	timer.start();
    	LI=backward(L,I);
    	timer.stop();
    	bt+=timer.elapsed_time();
    }
    Matrix<real_T> AI;
    for(size_t k=0;k<nr;++k) {
    	timer.start();
    	AI=forward(U,LI);
    	timer.stop();
    	ft+=timer.elapsed_time();
    }
    if(print) cerr<<AI<<endl<<endl;
    Matrix<real_T> IA=A*AI;
    if(print) cerr<<IA<<endl<<endl;
	cerr<<"************************************************************************************************"<<endl<<endl;
    cerr<<"Cholesky decomposition timing:\t"<<ct/nr<<endl;
    cerr<<"Backward solving timing:\t"<<bt/nr<<endl;
    cerr<<"Forward solving timing:\t"<<ft/nr<<endl;
//    cerr<<"Error:\t"<<IA.distance_max(IA,I)<<endl<<endl;
	cerr<<endl<<"************************************************************************************************"<<endl;
	return 0;
}
