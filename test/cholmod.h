#ifndef __CHOLMOD_TEST_H__
#define __CHOLMOD_TEST_H__

#include "../core/core.h"
#include "../linear-algebra/linear-algebra.h"
#include "../third-party/third-party.h"

#include <suitesparse/cholmod.h>
#include <suitesparse/SuiteSparse_config.h>

double cholmod(const __core__::MatrixCXS<double,int,__core__::cpuAllocator,__core__::CCS>& AC,bool supernodal=true) {
 	using namespace std;
	using namespace __core__;
	using namespace __third_party__;
	MatrixCXS<double,int,cpuAllocator,CCS> A,AMetis,AT,L;
	std::vector<int32_t> permutation,invPermutation;
	trimDiagonal(AT,AC);
	int metisCode=reorderPermuation(AT,permutation,invPermutation);
	error(metisCode==METIS_OK,"METIS failed, error code: "+std::to_string(metisCode));
	reorderMatrix(AMetis,AC,permutation.data(),invPermutation.data());
	auto tree=elimitationTree(AMetis);
	auto popermutation=postOrderTree(tree);
	reorderMatrix(A,AMetis,popermutation.data(),popermutation.data()+tree.size());
	cpu_timer timer;
	cholmod_common suiteSparse;
	cholmod_sparse ASP, *Lsparse;
	cholmod_factor *LSP;
	cholmod_start(&suiteSparse);

	suiteSparse.itype=CHOLMOD_INT;
    suiteSparse.final_asis = 1 ;
    if(supernodal)
    	suiteSparse.supernodal=CHOLMOD_SUPERNODAL;
    else
    	suiteSparse.supernodal=CHOLMOD_SIMPLICIAL;
//    suiteSparse.final_super = false ;
//    suiteSparse.final_ll = true ;
//    suiteSparse.final_pack = true ;
//    suiteSparse.final_monotonic=true ;
//    suiteSparse.final_resymbol=false ;
    suiteSparse.nmethods = 1 ;
    suiteSparse.method [0].ordering=CHOLMOD_NATURAL;
    suiteSparse.postorder = 0 ;

	ASP.nrow=A.n();
	ASP.ncol=A.m();
	ASP.nzmax=A.nzv();
	ASP.p=A.ptr();
	ASP.i=A.indxs();
	ASP.x=A.values();
	ASP.nz=NULL;
	ASP.z=NULL;
	ASP.stype=1;
	ASP.itype=CHOLMOD_INT;
	ASP.packed=true;
	ASP.sorted=true;
	ASP.dtype=CHOLMOD_DOUBLE;
	ASP.xtype=CHOLMOD_REAL;

    LSP = cholmod_analyze(&ASP,&suiteSparse) ;
    timer.start();
    cholmod_factorize(&ASP,LSP,&suiteSparse) ;
    timer.stop();
    if(supernodal)
    	cerr<<"CHOLMOD supernodal:"<<endl<<"\tElapsed time: "<<timer<<"\t"<<suiteSparse.cholmod_cpu_potrf_time+suiteSparse.cholmod_cpu_gemm_time+suiteSparse.cholmod_cpu_syrk_time+suiteSparse.cholmod_cpu_trsm_time+suiteSparse.cholmod_assemble_time<<endl<<endl;
    else
        cerr<<"CHOLMOD simplicial:"<<endl<<"\tElapsed time: "<<timer<<"\t"<<suiteSparse.cholmod_cpu_potrf_time<<endl<<endl;
    error(suiteSparse.status==CHOLMOD_OK,"Cholmod error, matrix is not PD");
    Lsparse=cholmod_factor_to_sparse(LSP, &suiteSparse);
    cholmod_free_factor(&LSP,&suiteSparse);
    cholmod_free_sparse(&Lsparse,&suiteSparse);
	cholmod_finish(&suiteSparse);
	return timer.elapsed_time();
}
#endif
