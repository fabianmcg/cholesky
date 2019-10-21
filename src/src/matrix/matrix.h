#ifndef __MATRIX_H__
#define __MATRIX_H__

#include <cassert>
#include <initializer_list>
#include <utility>
#include <iostream>
#include <vector>

namespace __core__ {
template <typename T> class Matrix {
private:
	size_t __rows__=0;
	size_t __cols__=0;
	std::vector<T> __data__;
public:
	Matrix();
	Matrix(const Matrix& matrix);
	Matrix(Matrix&& matrix);
	Matrix(size_t n,T v=0);
	Matrix(size_t n,size_t m,T v=0);
	Matrix(std::initializer_list<T> list,size_t n,size_t m);
	~Matrix();

	Matrix& operator=(const Matrix& matrix);
	Matrix& operator=(Matrix&& matrix);
	Matrix& operator=(std::initializer_list<T> list);
	Matrix& operator=(T v);

	T* operator*();
	const T* operator*() const;

	inline T& operator()(size_t i,size_t j);
	inline const T& operator()(size_t i,size_t j) const;
	inline T& operator[](size_t i);
	inline const T& operator[](size_t i) const;

	void clear();

	inline size_t n() const;
	inline size_t m() const;
	inline size_t rows() const;
	inline size_t cols() const;

	void resize(size_t i,T val=0);
	void resize(size_t i,size_t j,T val=0);
	void reshape(size_t i);
	void reshape(size_t i,size_t j);

	Matrix& add(const Matrix& A,const Matrix& B);
	Matrix add(const Matrix& B) const;
	Matrix& substract(const Matrix& A,const Matrix& B);
	Matrix substract(const Matrix& B) const;
	Matrix& multiply(const Matrix& A,const Matrix& B);
	Matrix multiply(const Matrix& B) const;
	Matrix& scalar_multiply(T x,const Matrix& B);
	Matrix scalar_multiply(T x) const;

	T norm() const;
	T distance(const Matrix& A,const Matrix& B) const;
	T norm_max() const;
	T distance_max(const Matrix& A,const Matrix& B) const;

	Matrix transpose() const;
	Matrix& transpose(const Matrix& B);

	void simetrize();

	std::ostream& print(std::ostream& ost,std::string separator=", ") const;
};
template <typename T> Matrix<T>::Matrix(){
}
template <typename T> Matrix<T>::Matrix(const Matrix& matrix): __rows__(matrix.__rows__),__cols__(matrix.__cols__),__data__(matrix.__data__) {
}
template <typename T> Matrix<T>::Matrix(Matrix&& matrix): __rows__(matrix.__rows__),__cols__(matrix.__cols__),__data__(std::move(matrix.__data__)) {
	matrix.__rows__=0;
	matrix.__cols__=0;
}
template <typename T> Matrix<T>::Matrix(size_t n,T v): __rows__(n),__cols__(n),__data__(std::vector<T>(n*n,v)) {
}
template <typename T> Matrix<T>::Matrix(size_t n,size_t m,T v): __rows__(n),__cols__(m),__data__(std::vector<T>(n*m,v)) {
}
template <typename T> Matrix<T>::Matrix(std::initializer_list<T> list,size_t n,size_t m): Matrix(n,m,0) {
	size_t i=0;
	for(auto it=list.begin();it!=list.end();++it)
		__data__[i++]=*it;
}
template <typename T> Matrix<T>::~Matrix() {
	__rows__=0;
	__cols__=0;
	__data__.clear();
}

template <typename T> Matrix<T>& Matrix<T>::operator=(const Matrix& matrix) {
	__rows__=matrix.__rows__;
	__cols__=matrix.__cols__;
	__data__=matrix.__data__;
	return *this;
}
template <typename T> Matrix<T>& Matrix<T>::operator=(Matrix&& matrix) {
	__rows__=matrix.__rows__;
	__cols__=matrix.__cols__;
	__data__=std::move(matrix.__data__);
	matrix.__rows__=0;
	matrix.__cols__=0;
	return *this;
}
template <typename T> Matrix<T>& Matrix<T>::operator=(std::initializer_list<T> list) {
	resize(list.size(),1,0);
	size_t i=0;
	for(auto it=list.begin();it!=list.end();++it)
		__data__[i++]=*it;
	return *this;
}
template <typename T> Matrix<T>& Matrix<T>::operator=(T v) {
	for(size_t i=0;i<__rows__;++i)
		for(size_t j=0;j<__cols__;++j)
			(*this)(i,j)=v;
	return *this;
}

template <typename T> T* Matrix<T>::operator*() {
	return __data__.data();
}
template <typename T> const T* Matrix<T>::operator*() const {
	return __data__.data();
}

template <typename T> inline T& Matrix<T>::operator()(size_t i,size_t j) {
	assert((i<__rows__)&&(j<__cols__));
	return __data__[i*__cols__+j];
}
template <typename T> inline const T& Matrix<T>::operator()(size_t i,size_t j) const {
	assert((i<__rows__)&&(j<__cols__));
	return __data__[i*__cols__+j];
}
template <typename T> inline T& Matrix<T>::operator[](size_t i) {
	assert(i<__data__.size());
	return __data__[i];
}
template <typename T> inline const T& Matrix<T>::operator[](size_t i) const {
	assert(i<__data__.size());
	return __data__[i];
}

template <typename T> inline void Matrix<T>::clear() {
	__rows__=0;
	__cols__=0;
	__data__.clear();
}

template <typename T> inline size_t Matrix<T>::n() const {
	return __rows__;
}
template <typename T> inline size_t Matrix<T>::m() const {
	return __cols__;
}
template <typename T> inline size_t Matrix<T>::rows() const {
	return __rows__;
}
template <typename T> inline size_t Matrix<T>::cols() const {
	return __cols__;
}

template <typename T> void Matrix<T>::resize(size_t i,T val) {
	resize(i,i,val);
}
template <typename T> void Matrix<T>::resize(size_t i,size_t j,T val) {
	__rows__=i;
	__cols__=j;
	__data__.resize(i*j,val);
}
template <typename T> void Matrix<T>::reshape(size_t i,size_t j) {
	assert((i*j)<=(__rows__*__cols__));
	__rows__=i;
	__cols__=j;
}
template <typename T> void Matrix<T>::reshape(size_t i) {
	reshape(i,i);
}

template <typename T> Matrix<T>& Matrix<T>::add(const Matrix& A,const Matrix& B) {
	assert((A.__rows__==B.__rows__)&&(A.__cols__==B.__cols__));
	(*this)=A;
	for(size_t i=0;i<__rows__;++i)
		for(size_t j=0;j<__cols__;++j)
			(*this)(i,j)+=B(i,j);
	return *this;
}
template <typename T> Matrix<T> Matrix<T>::add(const Matrix& B) const {
	Matrix<T> R;
	R.add(*this,B);
	return R;
}
template <typename T> Matrix<T>& Matrix<T>::substract(const Matrix& A,const Matrix& B) {
	assert((A.__rows__==B.__rows__)&&(A.__cols__==B.__cols__));
	(*this)=A;
	for(size_t i=0;i<__rows__;++i)
		for(size_t j=0;j<__cols__;++j)
			(*this)(i,j)-=B(i,j);
	return *this;
}
template <typename T> Matrix<T> Matrix<T>::substract(const Matrix& B) const {
	Matrix<T> R;
	R.substract(*this,B);
	return R;
}
template <typename T> Matrix<T>& Matrix<T>::multiply(const Matrix& A,const Matrix& B) {
	assert(A.__cols__==B.__rows__);
	this->resize(A.__rows__,B.__cols__,0);
	for(size_t i=0;i<__rows__;++i)
		for(size_t j=0;j<__cols__;++j)
			for(size_t k=0;k<A.__cols__;++k)
				(*this)(i,j)+=A(i,k)*B(k,j);
	return *this;
}
template <typename T> Matrix<T> Matrix<T>::multiply(const Matrix& B) const {
	Matrix<T> R;
	R.multiply(*this,B);
	return R;
}
template <typename T> Matrix<T>& Matrix<T>::scalar_multiply(T x,const Matrix& B) {
	(*this)=B;
	for(size_t i=0;i<__rows__;++i)
		for(size_t j=0;j<__cols__;++j)
			(*this)(i,j)*=x;
	return *this;
}
template <typename T> Matrix<T> Matrix<T>::scalar_multiply(T x) const {
	Matrix<T> R;
	R.scalar_multiply(x,*this);
	return R;
}

template <typename T> Matrix<T>& Matrix<T>::transpose(const Matrix& B) {
	this->resize(B.__cols__,B.__rows__,0);
	for(size_t i=0;i<__rows__;++i)
		for(size_t j=0;j<__cols__;++j)
			(*this)(i,j)=B(j,i);
	return *this;
}
template <typename T> Matrix<T> Matrix<T>::transpose() const {
	Matrix<T> R;
	R.transpose(*this);
	return R;
}

template <typename T> void Matrix<T>::simetrize() {
	for(size_t i=0;i<__rows__;++i)
		for(size_t j=0;j<__cols__;++j)
			if(i<j)
			(*this)(j,i)=(*this)(i,j);
}

template <typename T> T Matrix<T>::norm() const {
	T x=0;
	for(size_t i=0;i<__rows__;++i)
		for(size_t j=0;j<__cols__;++j) {
			T y=(*this)(i,j);
			x+=y*y;
		}
	return sqrt(x);
}
template <typename T> T Matrix<T>::distance(const Matrix& A,const Matrix& B) const {
	assert((A.__rows__==B.__rows__)&&(A.__cols__==B.__cols__));
	T x=0;
	for(size_t i=0;i<__rows__;++i)
		for(size_t j=0;j<__cols__;++j) {
			T y=A(i,j)-B(i,j);
			x+=y*y;
		}
	return sqrt(x);
}
template <typename T> T Matrix<T>::norm_max() const {
	T x=0;
	for(size_t i=0;i<__rows__;++i)
		for(size_t j=0;j<__cols__;++j) {
			T y=std::abs((*this)(i,j));
			x=std::max(x,y);
		}
	return x;
}
template <typename T> T Matrix<T>::distance_max(const Matrix& A,const Matrix& B) const {
	assert((A.__rows__==B.__rows__)&&(A.__cols__==B.__cols__));
	T x=0;
	for(size_t i=0;i<__rows__;++i)
		for(size_t j=0;j<__cols__;++j) {
			T y=std::abs(A(i,j)-B(i,j));
			x=std::max(x,y);
		}
	return x;
}

template <typename T>std::ostream& Matrix<T>::print(std::ostream& ost,std::string separator) const {
	for(size_t i=0;i<__rows__;++i) {
		for(size_t j=0;j<__cols__;++j) {
			if(j<__cols__-1)
				ost<<(*this)(i,j)<<separator;
			else
				ost<<(*this)(i,j);
		}
		if(i<__rows__-1)
			ost<<std::endl;
	}
	return ost;
}

template <typename T> inline Matrix<T> operator+(const Matrix<T>& A,const Matrix<T>& B) {
	Matrix<T> R;
	R.add(A,B);
	return R;
}
template <typename T> inline Matrix<T> operator*(const Matrix<T>& A,const Matrix<T>& B) {
	Matrix<T> R;
	R.multiply(A,B);
	return R;
}
template <typename T> inline Matrix<T> operator*(const T& x,const Matrix<T>& B) {
	Matrix<T> R;
	R.scalar_multiply(x,B);
	return R;
}
template <typename T> inline Matrix<T> operator*(const Matrix<T>& B,const T& x) {
	Matrix<T> R;
	R.scalar_multiply(x,B);
	return R;
}

template <typename T> std::ostream& operator<<(std::ostream& ost,const Matrix<T>& A) {
	return A.print(ost,", ");
}
}
#endif
