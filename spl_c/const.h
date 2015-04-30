#ifndef _SPL_CONST_
#define _SPL_CONST_

#include "common.h"
#include "spl_types.h"
#include <math.h>

///
/// \file  const.h
/// \brief Константы библиотеки SPL.
///

namespace spl {

// 
// Константы анатомические (e_ - ear)
// 

/// Максимальная длина основной мембраны внутреннего уха, мм
const int e_Xm = 35;

/// Нижняя частота восприятия звука человеком - 20 Гц
const int e_Flow = 20; 

/// Верхняя частота восприятия звука человеком - 20 кГц
const int e_Fhigh = 20000; 

//
// Константы системы фильтров (f_ - filter)
//

//@{
/// Коэффициенты, определяющие зависимость ширины критической полосы 
/// от шкалы частот фильтров.
const double f_alpha = 0.109;
const double f_beta = 69.095;
//@}

/// Коэффициент связи критической полосы и добротности
const double f_B = 1.96;

/// Константа, которая не расшифровывается, и равна 2.4 в модели.
const double f_C24 = 2.4;

// Максимальное число каналов
//const int f_Kmax = 1024;

// Максимальная длина окна фильтрации
//const int f_Imax = 2436;

// Максимальное число каналов в расширенной области 
// в алгоритме одновременной маскировки
//const int f_Mmax = 1085;


//
// Расчетные величины
//

/// Критическая полоса частот.
inline freq_t model_critical_band(freq_t f) { 
	return f_alpha * f + f_beta;
}

/// Коэффициент пропорциональности в 1/мм.
const double f_C = log(model_critical_band(e_Fhigh) / model_critical_band(e_Flow));

/// Координата точки на основной мембране в мм для заданной частоты анализа f
inline double model_f2x(freq_t f) { 
	return e_Xm * log(model_critical_band(e_Fhigh) / model_critical_band(f)) / f_C;
}

/// Частота анализа для заданной координаты точки x на основной мембране
inline freq_t model_x2f(double x) { 
	return (model_critical_band(e_Fhigh) / exp(x / e_Xm * f_C) - f_beta) / f_alpha;
}

/// Добротность фильтра с резонансной частотой f
inline double model_filter_quality(freq_t f) {
	return f_B * f / model_critical_band(f);
}

/// Параметр СКО в окне Гаусса для фильтрации.
/// Fr - резонансная частота фильтра
/// Fs - частота дискретизации сигнала
inline double filter_std(freq_t Fr, freq_t Fs) {
	return f_C24 * M_SQRT1_2 * model_filter_quality(Fr) * Fs / (2 * M_PI * Fr);
}

/// Параметр СКО в окне Гаусса для одновременной маскировки.
/// Fr - резонансная частота фильтра
/// delta - параметр, определяющий ширину маскирующей функции
inline double mask_std(freq_t Fr, double delta) {
//	return Fr / (f_C24 * model_filter_quality(Fr) * delta);
	return model_critical_band(Fr) / (f_C24 * f_B * delta);
}

/// Константа: 1 / sqrt(2*pi)
#define M_1_SQRT2PI 0.39894228040143267793994605993438

/// Формула окна Гаусса.
template<typename T>
inline T gauss_win(T x, T mu, T sigma) {
	T tmp = (x - mu) / sigma;
	return M_1_SQRT2PI / sigma * exp( - 0.5 * tmp*tmp);
}

/// Расчет ширины окна Гаусса.
template<typename T>
inline T gauss_border(T sigma, T err) {
	return ceil( sigma * sqrt( - 2 * log(err) ) );
}

} // namespace spl 

#endif//_SPL_CONST_
