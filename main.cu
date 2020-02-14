#include <cmath>
#include <iomanip>
#include <iostream>

#include "core/core.h"
#include "linear-algebra/linear-algebra.h"
#include "third-party/third-party.h"
#include "wavefront/wavefront.h"

using namespace std;
using namespace __core__;
using namespace __third_party__;

int main(int argc, char **argv) {
	cerr << "************************************************************"<< endl;
	cerr << "*                  Execution began                         *"<< endl;
	cerr << "************************************************************"<< endl;
	std::string blosumfn="./test/matrices/BLOSUM62.txt";
	size_t asize=16,bsize=16,px,py,its=1;
	px=__max__(asize/4,2),py=__max__(bsize/4,2);
	bool totalTime=false,debug=true;
	int threadnum=2;
	int gap=4;
	if(argc > 2) {
		asize=stol(argv[1]);
		bsize=stol(argv[2]);
	}
	if(argc > 4) {
		px=stol(argv[3]);
		py=stol(argv[4]);
	}
	if(argc > 5)
		threadnum=std::stoi(argv[5]);
	if(argc > 6)
		totalTime=std::stoi(argv[6])>0;
	if(argc > 7)
		its=std::stoi(argv[7]);
	if( argc > 8)
		debug=std::stoi(argv[8])>0;
	if( argc > 9)
		blosumfn=std::string(argv[9]);
	if(debug)
		its=1;
	auto blosum=readBlosum(blosumfn);
	string &alphabet=get<0>(blosum);
	alphabet.erase(std::remove(alphabet.begin(),alphabet.end(),'*'),alphabet.end());
	string As=generateSequence(asize,alphabet),Bs=generateSequence(bsize,alphabet);
	Array<int,cpuAllocator> A=translateSequence(As,get<1>(blosum)),B=translateSequence(Bs,get<1>(blosum));
	cerr<<endl<<"\tA size:\t"<<asize<<"\tB size:\t"<<bsize<<endl;
	cerr<<"\tPx size:\t"<<px<<"\t\tPy size:\t"<<py<<endl;
	cerr<<"\tThreadnum:\t"<<threadnum<<endl;
	cerr<<"\tReport total time:\t"<<totalTime<<endl;
	cerr<<"\tNumber of its:\t"<<its<<endl;
	cerr<<"\tBlosum matrix:\t"<<blosumfn<<endl<<endl;
	cerr << "************************************************************"<< endl;
	Matrix<int,cpuAllocator,1> S,P,&blosumM=get<2>(blosum);
	auto score=[&blosumM](int a,int b){return blosumM(a,b);};
	double se=0;
	for(size_t i=0;i<its;++i)
		se+=smithWaterman(B,A,S,score,gap);
	se=se/its;
	cerr<<"Smith-Waterman serial duration:\n\t"<<se<<endl;
	double pe=0;
	if(debug)
		pe+=smithWatermanPD(B,A,P,score,gap,px,py,threadnum,totalTime);
	else
		for(size_t i=0;i<its;++i)
			pe+=smithWatermanP(B,A,P,score,gap,px,py,threadnum,totalTime);
	pe=pe/its;
	cerr<<"Smith-Waterman parallel duration:\n\t"<<pe<<endl;
	cerr<<"Smith-Waterman speedup:\n\t"<<se/pe<<endl;
	size_t diff=0;
	for(size_t i=0;i<S.n();++i) {
		for(size_t j=0;j<S.m();++j) {
			diff+=__abs__(S(i,j)-P(i,j));
			if(diff>0)
				break;
		}
		if(diff>0)
			break;
	}
	if((S.n()!=P.n())&&(S.m()!=P.m()))
		diff+=1;
	cerr<<"Difference: "<<diff<<endl;
	cerr << "************************************************************"<< endl;
	cerr << "*                  Execution ended                         *"<< endl;
	cerr << "************************************************************"<< endl;
	return 0;
}
