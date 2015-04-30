#ifndef _SPL_CONV_
#define _SPL_CONV_

///
/// \file  conv.h
/// \brief Функции для вычисления свертки.
///

#include "common.h"

namespace spl {

/// Размер окна циклической свертки.
/// Чем больше окно циклической свертки, 
///  тем более эффективны алгоритмы фильтрации длинных сигналов, 
///  и тем менее эффективны алгоритмы фильтрации коротких сигналов.
const int CONV_WIN_SIZ = 8192;

/// Функция для выравнивания памяти по границе 16 байт.
/// Это нужно для ускорения вычислений при использовании SSE и т.п.

template<typename T>
T *SPL_MEMORY_ALIGN(T *p) { 
	return (T *)(((unsigned long long)(p) + 15) & ~15); 
}

//@{
/// Функции управления памятью, предназначенные для использования со сверткой
void EXPORT *conv_alloc_low(size_t N);
void EXPORT conv_free(void *m);

template<typename T>
T *conv_alloc(size_t siz) {
  return (T*) conv_alloc_low(siz * sizeof(T));
}
//@}

/// Вычисление вектора A = FFT(a).
/// Необходимо для быстрого вычисления двойной циклической свертки (a * (b,c))
void EXPORT cconv_calc_A(IN double *a, OUT double *Ar, OUT double *Ai);

/// Вычисление вектора (B,C) = FFT((b,c)).
/// Необходимо для быстрого вычисления двойной циклической свертки (a * (b,c))
void EXPORT cconv_calc_BC(IN double *b, IN double *c, OUT double *B, OUT double *C);

/// Нормализация коэффициентов фильтра.
/// Необходимо для быстрого вычисления двойной циклической свертки (a * (b,c))
/// Поскольку используемая библиотека FFTW домножает результат IFFT на размер массива, 
///  то разумно перед фильтрацией разделить коэффициенты фильтра на размер массива (CONV_WIN_SIZ).
void EXPORT cconv_normalize(INOUT double *x, size_t n);

/// Вычисление модуля каждого элемента комплексного вектора.
/// Для скорости вычисляет квадрат модуля.
void complex_abs_split(size_t N, IN double *Ar, IN double *Ai, OUT double *Am);

/// Двойная циклическая свертка, рассчитанная по быстрому алгоритму.
/// Такая циклическая свертка используется для быстрой цифровой фильтрации
///  по алгоритму пересечения с накоплением (overlap-save), 
///  в том числе при фильтрации сигнала и вычислении одновременной маскировки.
void EXPORT cconv(
 IN  double *Ar, IN  double *Ai, 
 IN  double *C,  IN  double *B, 
 OUT double *ab, OUT double *ac
);

} // namespace spl 

#endif//_SPL_CONV_
