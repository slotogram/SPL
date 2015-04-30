#ifndef _SPL_MATRIX_
#define _SPL_MATRIX_

///
/// \file  matrix.h
/// \brief ��������� ����� �������.
/// 
/// �������������� ����� ���������� �������� ����������, ������� ��� ������ � ���������.
/// ����������� ������� (2-������, 3-������, 4-������ � �.�.) ������ ��������� ���������� (������).
/// ���� ������� ������� (2�3, 3�5�2 � �.�.) ����� ���������� �������������.
///

#include "common.h" // size_t

#include <stdarg.h>

/// ������ ��� ������������� ������������� ���������� ���������� � �������.
/// �������� ��������� � ������� �������� ������ "name".
#define va2array(T, name, first, D) \
  va_list vl; \
  va_start(vl, first); \
  name[0] = first; \
  for(int i = 1; i < D; i++) { \
   name[i] = va_arg(vl, T); \
  } \
  va_end(vl);

/// ������ ��� ������������� ������������� ���������� ���������� � �������.
/// ������ ������ "name" � �������� � ���� ���������.
#define va2array2(T, name, first, D) \
	T name[D]; \
	va2array(T, name, first, D)



/// ��������������� ������� ������� ��������� ������� �������.
size_t sub2ind(int D, const size_t dim[], const size_t sub[]);

/// ��������������� ������� ������� ��������� ������� �������.
size_t dim2siz(int D, const size_t dim[]);


/// ����� ����������� �������.
/// 
/// ���������� �������� ���������� �� ����������� ������.
/// ��� ������� ������ � ����������� �� ����������� �������.
/// 
/// ����� �� ������ �����-���� ������ �� ��������� ��� ������������ ������.

template<typename T, int D>
class Matrix {
public:
  /// ����������� (��������� ������ ������������).
  /// ����������� ������ ����������� � ������� ���������� �� ���������� (������� ���������� ������).
  /// ���, ��������� ����������� ������ ���������� ������ �������, ������ - ������ �������, � �.�.
  Matrix(T *data, const size_t dims[]): 
    _data(data), _size(dim2siz(D, dims))
  {
    memcpy(_dim, dims, sizeof(size_t) * D);
  }
  /// ����������� (��������� ����������� � ������ ����������).
  Matrix(T *data, size_t first, ...): 
    _data(data)
  {
  	va2array(size_t, _dim, first, D);
    _size = dim2siz(D, _dim);
  }
  
  //@{
  /// ����������� �������� � �����������.
  /// matrix(1,2) ��������� �������, ������� ������������� � ������ ������, ������ �������.
  /// ��� ���� ������� ������ ���� ���������.
  /// 
  /// ��������: ����������� ������� �� �����������! 
  /// (�� ���������� ��������� � ����� ����������� ���������� ������ �������).
  /// ���� ������� �������� D-������, �� � ��� ������� ������ ���� �������� ����� D ���������� ���� size_t.
  /// ���� ����� �������� ������ ����������, �� � ������ ����� ������� ������ ������ D ����������.
  /// ���� ����� �������� ������ ����������, �� �������� ��������� ������ ��������� � ������.
  const T& operator()(size_t first, ...) const { 
    va2array2(size_t, sub, first, D);
    return _data[sub2ind(sub)]; 
  }

  T& operator()(size_t first, ...) { 
    va2array2(size_t, sub, first, D);
    return _data[sub2ind(sub)]; 
  }
  //@}

  //@{
  /// ���������� �������� � �����������
  const T& operator[](size_t i) const { return _data[i]; }
  T& operator[](size_t i) { return _data[i]; }
  //@}

  //@{
  /// �������������� � ���� T *
  operator const T *() const { return _data; }
  operator T *() { return _data; }
  //@}

  /// ���������� ������ ������������ �������.
  const size_t *dims() const { return _dim; }
  /// ���������� �������� ����������� �������.
  size_t dim(int i) const { return _dim[i]; }
  /// ���������� ���������� ��������� �������.
  size_t size() const { return _size; }
  /// ���������� ��������� �� ������ ������ �������.
  T *ptr() { return _data; }

private:
  /// ���������. 
  /// ������ � ��������� �������. ����������� ������� ������� �� ������� � �������
  size_t _dim[D];

  /// ������ ������.
  size_t _size;

  /// ��������� �� ������.
  /// ���������� ������ �������. ������ ����� ������ _dim[0] * _dim[1] * ...
  T *_data;

  /// ������ ��������� �������
  size_t sub2ind(size_t sub[]) const {
    switch(D) { // ��� ����������� (�������� ���������� ������ �������� �����)
    case 1: return sub[0];
    case 2: return _dim[1] * sub[0] + sub[1];
    case 3: return _dim[2] * (_dim[1] * sub[0] + sub[1]) + sub[2];
    default: return ::sub2ind(D, _dim, sub);
    }
  }

};

/// ������� ��� �������� ��������� ��������� �������.

template<typename T>
Matrix<T,2> matrix_ptr(T *data, size_t A, size_t B) {
  return Matrix<T,2>(data, A, B);
}

/// ������� ��� �������� ��������� ���������� �������.

template<typename T>
Matrix<T,3> matrix_ptr(T *data, size_t A, size_t B, size_t C) {
  return Matrix<T,3>(data, A, B, C);
}

/// ������� ��� �������� ��������� ������������� �������.

template<typename T>
Matrix<T,4> matrix_ptr(T *data, size_t A, size_t B, size_t C, size_t D) {
  return Matrix<T,4>(data, A, B, C, D);
}

//@{
///
/// ��������� ������� ������� - reshape.
/// 

template<typename T, int D1, int D2>
Matrix<T,D2> reshape(Matrix<T,D1>& mtx1, const size_t dims2[]) {
  size_t siz2 = dim2siz(D2, dims2);
  if(mtx1.size() != siz2) 
    throw "Matrix sizes must agree on reshape";
  return Matrix<T,D2>(mtx1.ptr(), dims2);
}

template<typename T, int D1, int D2>
Matrix<T,D2> reshape(Matrix<T,D1>& mtx1, const size_t first, ...) {
  va2array2(size_t, dims, first, D2);
  return reshape<T,D1,D2>(mtx1, dims);
}
//@}

///
/// ���������������� ��������� �������.
/// ��� ����������� ���������������� ������� ������ ���� �������� MxN � NxM, 
///  ��� �� ������� ���� � �������� ������� ������ ���� �����.
/// ������� ����� ����� ���� ������ �����, 
///  �������, ����� ����������� �������� "=" ��� �������������� ����� ����.
///

template<typename T1, typename T2>
bool transpose(const Matrix<T1,2>& m, Matrix<T2,2>& t) {
  // ���������, ���� �� ����� � �������� �������.
  if(m.dim(0) > t.dim(1)
  || m.dim(1) > t.dim(0))
    return false;
  // ���������� ����������������
  for(size_t i = 0; i < m.dim(0); i++) {
  for(size_t j = 0; j < m.dim(1); j++) {
    t(j,i) = m(i,j);
  }}
  return true;
}

#endif//_SPL_MATRIX_
