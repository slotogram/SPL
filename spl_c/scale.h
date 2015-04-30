#ifndef _SPL_SCALE
#define _SPL_SCALE

///
/// \file  scale.h
/// \brief ��������� � �������, ����������� � ��������� ������.
///

#include "common.h"
#include "spl_types.h"

namespace spl {

///
/// ���������, ������������ ��� ��������� ����� ������ \ref freq_scale.
///

struct scale_parameters {

	//@{
	/// ������ ������� � ���������� ������ ����������.
	/// ���������� ������� K = b - a + 1;
	int a, b; //@}

	//@{
	/// ������ ���� ��������� ������� ����������.
	int i, j; //@}

	//@{
	/// ��� ��������� �������, ��������������� ��������� ������� ���������� i � j.
	freq_t Fi, Fj; //@}

};

/// ������� ��� ���������� ����� � ����.
bool EXPORT save_scale(
	IN freq_scale& sc, ///< �����, ������� ���������� ���������.
	const char *file   ///< ���� � �����, � ������� ���������� ��������� �����.
);

//
// ������� ��������� �����
// 

/// ������� ����� - ��� ����� �������� ������� ���������� �������� ��������� �������.
typedef double (*scale_func)(freq_t);

/// �������� ������� ����� - ��� �������� ������� ���������� �������.
typedef freq_t (*rscale_func)(double);

/// ������� ��� ��������� ����� �� �������� �������� �����. 
bool EXPORT scale_generate(
	OUT freq_scale& sc,     ///< �����, ������� ����� �������������
	IN scale_parameters& p, ///< ��������� ��������� �����
	IN scale_func f2x, ///< ������ ������� �����
	IN rscale_func x2f ///< �������� ������� �����
);

/// ����� ������������ �����
enum scale_form {
	scale_unknown=0,///< ����������� ����� �����
	scale_linear,  ///< �������� �����
	scale_log,     ///< ��������������� �����
	scale_mel,     ///< ����� ���
	scale_bark,    ///< ����� ����
	scale_model,   ///< ����� ������ �����
};

/// ������� ��� ��������� ����� �� �������� �����.
bool EXPORT scale_generate(
	OUT freq_scale& sc,     ///< �����, ������� ����� �������������
	IN scale_parameters& p, ///< ��������� ��������� �����
	IN scale_form form      ///< ����� �����
);

/// ����� ����������� ������� ��������� �����.
typedef bool (*scale_generate_func)(OUT freq_scale& sc, IN scale_parameters& p);

/// ������� ��� ��������� �������� �����.
bool EXPORT scale_generate_linear(OUT freq_scale& sc, IN scale_parameters& p);

/// ������� ��� ��������� ����� Bark.
bool EXPORT scale_generate_bark(OUT freq_scale& sc, IN scale_parameters& p);

/// ������� ��� ��������� ����� Mel.
bool EXPORT scale_generate_mel(OUT freq_scale& sc, IN scale_parameters& p);

/// ������� ��� ��������� ��������������� ����� ������.
bool EXPORT scale_generate_log(OUT freq_scale& sc, IN scale_parameters& p);

/// ������� ��� ��������� ����� ������ �� ������ ���.
bool EXPORT scale_generate_model(OUT freq_scale& sc, IN scale_parameters& p);

///
/// �� ��, ��� � \ref freq_scale, �� � ����������� ������ (��� C++).
///

struct EXPORT scale_class: freq_scale {

	/// �����������. ���������� ����� �� �������� ����� �����.
	scale_class(const scale_parameters& p, scale_form form);

	/// �����������. ���������� ����� �� �������� �������� �����.
	scale_class(const scale_parameters& p, scale_func f, rscale_func r);

	/// �����������. ��������� ����� �� ��������� �����.
	scale_class(const char *filename);

	/// ���������� ����� � ����. 
	/// �������� ������� save_scale().
	bool save(const char *file) { return save_scale(*this, file); }

	/// ����������. ����������� ������, ���������� ��� �����.
	~scale_class();
private:
	void init(int K);

};

/// ������� ������ ������� � �����. 
/// ��������������, ��� ����� ������������.
int EXPORT scale_find_freq(const freq_scale& sc, freq_t F);

/// ����������� ���������� �����, �� �������� ������������ ��� �����.
const int DEFAULT_SCALE_ESTIMATE_ACCURACY = 5;

/// ������� �������� ���� �����.
bool EXPORT is_scale_of_form(const freq_scale& sc, enum scale_form form, int npoints = DEFAULT_SCALE_ESTIMATE_ACCURACY);

/// ������� ����������� ����� �����.
enum scale_form EXPORT scale_form_estimate(const freq_scale& sc, int npoints = DEFAULT_SCALE_ESTIMATE_ACCURACY);

/// ������ �������� ���� �����.
/// �������� ����������� �� ���� �������� ���� ����� �����.
inline bool is_scale_of_form_accurate(const freq_scale& sc, enum scale_form form) {
	return is_scale_of_form(sc, form, sc.K - 2);
}

/// ������ ����������� ����� �����.
/// �������� ����������� �� ���� �������� ���� ����� �����.
inline enum scale_form scale_form_estimate_accurate(const freq_scale& sc) {
	return scale_form_estimate(sc, sc.K - 2);
}

/// ������������ ����� � ������ �����.
freq_t EXPORT scale_interp_k2f(const freq_scale& sc, scale_form form, double k);

/// ������������ ����� � ������ �������.
double EXPORT scale_interp_f2k(const freq_scale& sc, scale_form form, freq_t f);


} // namespace spl

#endif//_SPL_SCALE
