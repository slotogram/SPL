///
/// \file  conv.cpp
/// \brief ������ ������� ����������� �������
///
/// ������� ����������� ������� ����������� ����� ������� �������������� ����� (FFT)
///  � �������������� ���������� FFTW.
///

#include "conv.h"
using spl::SPL_MEMORY_ALIGN;

#include <fftw3.h>
#define _USE_MATH_DEFINES // M_PI, etc
#include <math.h>
#include <stdio.h>
#include <malloc.h>

namespace {

#define SPL_FFTW_WISDOM_FILE "spl-fftw-wisdom"
fftw_plan plan_A, plan_BC, plan_abc;

// ����-�������������
struct my_init {
  my_init() {
    using spl::CONV_WIN_SIZ;
    FILE *f;
    
    // try to open wisdom file:
    f = fopen(SPL_FFTW_WISDOM_FILE, "r");
    if(f) {
     char wisdom[10000];
     fread(wisdom, 1, 10000, f);
     fftw_import_wisdom_from_string(wisdom);
     fclose(f); 
    }

    double array[CONV_WIN_SIZ*4 + 2];
    double *ri = SPL_MEMORY_ALIGN(array);
    double *ii = ri + CONV_WIN_SIZ;
    double *ro = ri + 2*CONV_WIN_SIZ;
    double *io = ri + 3*CONV_WIN_SIZ;
    
    fftw_iodim dim, howmany;
    dim.n = CONV_WIN_SIZ; dim.is = dim.os = 1;
	howmany.n = 256; howmany.is = howmany.os = CONV_WIN_SIZ;

    // generate wisdom:
    plan_A = fftw_plan_guru_split_dft_r2c(1, &dim, 0, 0, ri, ro, io, FFTW_MEASURE);
    plan_BC = fftw_plan_guru_split_dft(1, &dim, 0, 0, ri, ii, ro, io, FFTW_MEASURE);
	plan_abc = plan_BC;
//    plan_abc = fftw_plan_guru_split_dft(1, &dim, 1, &howmany, ri, ii, ro, io, FFTW_ESTIMATE);
    
    // try to save wisdom:
	f = fopen(SPL_FFTW_WISDOM_FILE, "w");
	if(f) {
		fftw_export_wisdom(myputc, f);
		fclose(f);
	}
  }

private:
	static void myputc(char x, void *f) {
		putc(x, (FILE *)f);
	}

} _my_init;

}

namespace spl {

//
// �������� � �������
//

void *conv_alloc_low(size_t N) {
  return fftw_malloc(N);
}

void conv_free(void *x) {
  fftw_free(x);
}

//
// ��������� ����������� ��������
//

/// ��������� ����������� �����.
/// (Ar + iAi) * (Br + iBi) = ArBr - AiBi + i(ArBi + AiBr)
/// ��������� in-place ���������
inline void complex_mul_split(
 IN  double& Ar,  IN  double& Ai, 
 IN  double& Br,  IN  double& Bi, 
 OUT double& ABr, OUT double& ABi) 
{
 double r = Ar * Br;
 double i = Ai * Bi;
 ABi = (Ar + Ai) * (Br + Bi) - r - i;
 ABr = r - i;
}

/// ���������� ������ ������� �������� ������������ �������.
void complex_abs_split(size_t N, IN double *Ar, IN double *Ai, OUT double *Am) {
 for(size_t i = N; i--; ) {
  Am[i] = Ar[i] * Ar[i] + Ai[i] * Ai[i];
 }
}


//
// ��������� ����� �������� ���������� ������� ����������� ������� ����� FFT
//

/// ���������� ������� A = FFT(a).
void cconv_calc_A(IN double *a, OUT double *Ar, OUT double *Ai) {
  fftw_execute_split_dft_r2c(plan_A, const_cast<double*>(a), Ar, Ai);
  // ����������� ���������� FFTW: ��-�� ��������� ���������� ������ �������� A
  // ������� ��������� ������ ��������:
  for(int i = CONV_WIN_SIZ/2 + 1; i < CONV_WIN_SIZ; i++) {
	  Ar[i] =  Ar[CONV_WIN_SIZ - i];
	  Ai[i] = -Ai[CONV_WIN_SIZ - i];
  }
}

/// ���������� ������� (B,C) = FFT((b,c)).
void cconv_calc_BC(IN double *b, IN double *c, OUT double *B, OUT double *C) {
  fftw_execute_split_dft(plan_BC, const_cast<double*>(b), const_cast<double*>(c), B, C);
}

/// ���������� ������� ABC: (AB,AC) = A * (B,C).
/// ���������� ��� �������� ���������� ������� ����������� ������� (a * (b,c))
void cconv_calc_ABC(
 IN  double *Ar, IN  double *Ai, 
 IN  double *B,  IN  double *C, 
 OUT double *AB, OUT double *AC) 
{
	for(int i = 0; i < CONV_WIN_SIZ; i++) {
		complex_mul_split(Ar[i], Ai[i], B[i], C[i], AB[i], AC[i]);
	}
}

/// ���������� ������� abc: (ab,ac) = IFFT((AB,AC)).
/// ���������� ��� �������� ���������� ������� ����������� ������� (a * (b,c))
void cconv_calc_abc(IN double *AB, IN double *AC, OUT double *ab, OUT double *ac) {
  // ����� �������� �������� DFT, re & im �������� ������� (FFTW reference)
  fftw_execute_split_dft(plan_abc, const_cast<double*>(AC), const_cast<double*>(AB), ac, ab);
}

/// ������������ ������������� �������.
void cconv_normalize(INOUT double *x, size_t n) {
  for(size_t j = 0; j < n; j++) 
    x[j] /= CONV_WIN_SIZ;
}

/// ������� ����������� �������, ������������ �� �������� ���������.
void cconv(
 IN  double *Ar, IN  double *Ai, 
 IN  double *B,  IN  double *C, 
 OUT double *ab, OUT double *ac) 
{
  // 0. ������� � ������
  double array[2*CONV_WIN_SIZ + 2];
  double *ABCr = SPL_MEMORY_ALIGN(array);
  double *ABCi = ABCr + CONV_WIN_SIZ;

  // 1. ABC = A .* BC
  cconv_calc_ABC(Ar, Ai, B, C, ABCr, ABCi);

  // 2. (ab,ac) = ifft(ABC)
  cconv_calc_abc(ABCr, ABCi, ab, ac);
}


} // namespace spl 
