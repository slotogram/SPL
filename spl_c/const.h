#ifndef _SPL_CONST_
#define _SPL_CONST_

#include "common.h"
#include "spl_types.h"
#include <math.h>

///
/// \file  const.h
/// \brief ��������� ���������� SPL.
///

namespace spl {

// 
// ��������� ������������� (e_ - ear)
// 

/// ������������ ����� �������� �������� ����������� ���, ��
const int e_Xm = 35;

/// ������ ������� ���������� ����� ��������� - 20 ��
const int e_Flow = 20; 

/// ������� ������� ���������� ����� ��������� - 20 ���
const int e_Fhigh = 20000; 

//
// ��������� ������� �������� (f_ - filter)
//

//@{
/// ������������, ������������ ����������� ������ ����������� ������ 
/// �� ����� ������ ��������.
const double f_alpha = 0.109;
const double f_beta = 69.095;
//@}

/// ����������� ����� ����������� ������ � �����������
const double f_B = 1.96;

/// ���������, ������� �� ����������������, � ����� 2.4 � ������.
const double f_C24 = 2.4;

// ������������ ����� �������
//const int f_Kmax = 1024;

// ������������ ����� ���� ����������
//const int f_Imax = 2436;

// ������������ ����� ������� � ����������� ������� 
// � ��������� ������������� ����������
//const int f_Mmax = 1085;


//
// ��������� ��������
//

/// ����������� ������ ������.
inline freq_t model_critical_band(freq_t f) { 
	return f_alpha * f + f_beta;
}

/// ����������� ������������������ � 1/��.
const double f_C = log(model_critical_band(e_Fhigh) / model_critical_band(e_Flow));

/// ���������� ����� �� �������� �������� � �� ��� �������� ������� ������� f
inline double model_f2x(freq_t f) { 
	return e_Xm * log(model_critical_band(e_Fhigh) / model_critical_band(f)) / f_C;
}

/// ������� ������� ��� �������� ���������� ����� x �� �������� ��������
inline freq_t model_x2f(double x) { 
	return (model_critical_band(e_Fhigh) / exp(x / e_Xm * f_C) - f_beta) / f_alpha;
}

/// ����������� ������� � ����������� �������� f
inline double model_filter_quality(freq_t f) {
	return f_B * f / model_critical_band(f);
}

/// �������� ��� � ���� ������ ��� ����������.
/// Fr - ����������� ������� �������
/// Fs - ������� ������������� �������
inline double filter_std(freq_t Fr, freq_t Fs) {
	return f_C24 * M_SQRT1_2 * model_filter_quality(Fr) * Fs / (2 * M_PI * Fr);
}

/// �������� ��� � ���� ������ ��� ������������� ����������.
/// Fr - ����������� ������� �������
/// delta - ��������, ������������ ������ ����������� �������
inline double mask_std(freq_t Fr, double delta) {
//	return Fr / (f_C24 * model_filter_quality(Fr) * delta);
	return model_critical_band(Fr) / (f_C24 * f_B * delta);
}

/// ���������: 1 / sqrt(2*pi)
#define M_1_SQRT2PI 0.39894228040143267793994605993438

/// ������� ���� ������.
template<typename T>
inline T gauss_win(T x, T mu, T sigma) {
	T tmp = (x - mu) / sigma;
	return M_1_SQRT2PI / sigma * exp( - 0.5 * tmp*tmp);
}

/// ������ ������ ���� ������.
template<typename T>
inline T gauss_border(T sigma, T err) {
	return ceil( sigma * sqrt( - 2 * log(err) ) );
}

} // namespace spl 

#endif//_SPL_CONST_
