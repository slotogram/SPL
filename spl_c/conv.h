#ifndef _SPL_CONV_
#define _SPL_CONV_

///
/// \file  conv.h
/// \brief ������� ��� ���������� �������.
///

#include "common.h"

namespace spl {

/// ������ ���� ����������� �������.
/// ��� ������ ���� ����������� �������, 
///  ��� ����� ���������� ��������� ���������� ������� ��������, 
///  � ��� ����� ���������� ��������� ���������� �������� ��������.
const int CONV_WIN_SIZ = 8192;

/// ������� ��� ������������ ������ �� ������� 16 ����.
/// ��� ����� ��� ��������� ���������� ��� ������������� SSE � �.�.

template<typename T>
T *SPL_MEMORY_ALIGN(T *p) { 
	return (T *)(((unsigned long long)(p) + 15) & ~15); 
}

//@{
/// ������� ���������� �������, ��������������� ��� ������������� �� ��������
void EXPORT *conv_alloc_low(size_t N);
void EXPORT conv_free(void *m);

template<typename T>
T *conv_alloc(size_t siz) {
  return (T*) conv_alloc_low(siz * sizeof(T));
}
//@}

/// ���������� ������� A = FFT(a).
/// ���������� ��� �������� ���������� ������� ����������� ������� (a * (b,c))
void EXPORT cconv_calc_A(IN double *a, OUT double *Ar, OUT double *Ai);

/// ���������� ������� (B,C) = FFT((b,c)).
/// ���������� ��� �������� ���������� ������� ����������� ������� (a * (b,c))
void EXPORT cconv_calc_BC(IN double *b, IN double *c, OUT double *B, OUT double *C);

/// ������������ ������������� �������.
/// ���������� ��� �������� ���������� ������� ����������� ������� (a * (b,c))
/// ��������� ������������ ���������� FFTW ��������� ��������� IFFT �� ������ �������, 
///  �� ������� ����� ����������� ��������� ������������ ������� �� ������ ������� (CONV_WIN_SIZ).
void EXPORT cconv_normalize(INOUT double *x, size_t n);

/// ���������� ������ ������� �������� ������������ �������.
/// ��� �������� ��������� ������� ������.
void complex_abs_split(size_t N, IN double *Ar, IN double *Ai, OUT double *Am);

/// ������� ����������� �������, ������������ �� �������� ���������.
/// ����� ����������� ������� ������������ ��� ������� �������� ����������
///  �� ��������� ����������� � ����������� (overlap-save), 
///  � ��� ����� ��� ���������� ������� � ���������� ������������� ����������.
void EXPORT cconv(
 IN  double *Ar, IN  double *Ai, 
 IN  double *C,  IN  double *B, 
 OUT double *ab, OUT double *ac
);

} // namespace spl 

#endif//_SPL_CONV_
