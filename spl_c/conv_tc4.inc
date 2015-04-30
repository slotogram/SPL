///
/// \file  conv.cpp
/// \brief Модуль быстрой циклической свертки
///
/// Быстрая циклическая свертка реализована через быстрое преобразование Фурье (FFT)
///  с использованием библиотеки FFTW.
///

#define CONV_VECTORIZE

#include "conv.h"
#include "../conv_c/conv.h"
#include <algorithm>
using std::copy;

namespace spl {

void *conv_alloc_low(size_t N) {
  return spl_alloc_low(N);
}

void conv_free(void *x) {
  spl_free(x);
}

// подготовка коэффициентов фильтра
void cconv_precalc_b(IN double *b, OUT double *B) {
 for(int i = 0; i < CONV_WIN_SIZ; i++) {
  B[i] = b[2*i]; // cos
  B[i+CONV_WIN_SIZ] = b[2*i+1]; // sin
 }
}


//
// a[CONV_WIN_SIZ]
// BC[CONV_WIN_SIZ][2]
// abc[CONV_WIN_SIZ]
//

void cconv(IN double *a, IN double *b, OUT double *abc) {

 double ab[CONV_WIN_SIZ];
 double ac[CONV_WIN_SIZ];
 const double *c = b + CONV_WIN_SIZ;
 base<double,CONV_WIN_SIZ>::cconv(a, b, ab);
 base<double,CONV_WIN_SIZ>::cconv(a, c, ac);
 for(int i = 0; i < CONV_WIN_SIZ; i++) {
  abc[i] = ab[i] * ab[i] + ac[i] * ac[i];
 }
 printf(".");

}

} // namespace spl 
