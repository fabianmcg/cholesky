#ifndef __SPARSE_MATRIX_OPERATIONS_H__
#define __SPARSE_MATRIX_OPERATIONS_H__

#include "../matrix-formats/matrix-formats.h"

namespace __core__ {
namespace __linear_algebra__ {
namespace __matrix_operations__ {
template <typename T>
void transpose(MatrixMap<T>& dmatrix,const MatrixMap<T>& smatrix) {
	dmatrix.clear();
	dmatrix.reshape(smatrix.cols(),smatrix.rows());
	for(auto it=(*smatrix).cbegin();it!=(*smatrix).cend();++it)
		dmatrix(it->first[1],it->first[0])=it->second;
}
template <typename T>
MatrixMap<T> transpose(const MatrixMap<T>& matrix) {
	MatrixMap<T> result;
	transpose(result,matrix);
	return result;
}
template <typename T>
void multiply(MatrixMap<T>& matrix,const MatrixMap<T>& A,const MatrixMap<T>& B) {
	error(A.cols()==B.rows(),"Invalid shapes for the matrices.",RUNTIME_ERROR,throw_error);
	matrix.clear();
	matrix.reshape(A.rows(),B.cols());
	for(auto it=(*A).cbegin();it!=(*A).cend();++it) {
		for(auto jt=(*B).lower_bound(size_2A({it->first[1],0}));jt!=(*B).cend();++jt) {
			if(jt->first[0]==it->first[1])
				matrix(it->first[0],jt->first[1])+=it->second*jt->second;
			else
				break;
		}
	}
}
template <typename T>
MatrixMap<T> multiply(const MatrixMap<T>& A,const MatrixMap<T>& B) {
	MatrixMap<T> result;
	multiply(result,A,B);
	return result;
}
template <typename T>
T maxDifference(const MatrixMap<T>& A,const MatrixMap<T>& B) {
	error(A.cols()==B.cols()&&A.rows()==B.rows(),"Invalid shapes for the matrices.",RUNTIME_ERROR,throw_error);
	MatrixMap<T> tmp=A;
	T result=0;
	for(auto it=(*B).cbegin();it!=(*B).cend();++it)
		tmp[it->first]-=it->second;
	for(auto it=(*tmp).cbegin();it!=(*tmp).cend();++it)
		result=__max__(result,__abs__(it->second));
	return result;
}
}
}
}
#endif
