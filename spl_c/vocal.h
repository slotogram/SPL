#ifndef _SPL_VOCAL_
#define _SPL_VOCAL_

///
/// \file  vocal.h
/// \brief ������� ��� ��������� ������ �� �������� ����������������.
///
/// �������� � ���� ������� ��� ����������� ������� ��������� ���� (���) 
///  � �������������� ����������� �������� �� ��������������/����������������.
///

#include "common.h"
#include "spl_types.h"
#include "mask.h"
#include "io.h"

namespace spl {

/// ������� ���������� (������������ ��� ����������� ������� ��������� ����).

struct mask_templates {

	//@{
	/// ������ � ��������� �����, ��� ������� ���� �������.
	/// ��������� \ref mask_templates �������� ������� ��� ������� k1 ... k2.
	/// ����� �������, � ��������� \ref mask_templates ���������� (k2 - k1 + 1) ��������.
	int k1, k2;

	/// ������ ���������� �������.
	int K;

	/// ���������� ��������� �������.
	int Nt;

	/// ���������� �������.
	/// ������ ������������ ����� ������� ������� (k2 - k1 + 1) x Nt.
	int *T;

};

/// ��������� ��������� �������� ��� ��������� ��� - ������� ��������� ���� (��. \ref mask_templates).

struct pitch_parameters {

	/// ���������� ����������� �������� (��� ����������� ���: 2).
	int Nh;

	//@{ 
	/// ������ � ������� ������� ����������� ���.
	freq_t F1, F2; //@}

};

/// ����������� �������� ���������� ��������� ��������, ��� ����������� ���.
const pitch_parameters DEFAULT_PITCH_PARAMETERS = { 2, 75, 400 };


/// ��������� ����������� �� �������� ����������������.

struct vocal_parameters {

	/// ����������� ������������ ��������������� �������� (� ��������)
	double minV;

	/// ����������� ������������ ����������������� �������� (� ��������)
	double minNV;

};

/// ����������� �������� ���������� ����������� �� �������� ����������������.
const vocal_parameters DEFAULT_VOCAL_PARAMETERS = { 0.030, 0.030 };

/// ��������� ������� ����� ��� ����������� ������� ��������� ����.
bool EXPORT pitch_templates_generate(
	IN freq_scale& sc,        ///< ����� ������.
	OUT mask_templates& tpl,  ///< ������ �����.
	IN mask_parameters& pm,   ///< ��������� ��������� ����������� �������.
	IN pitch_parameters& p    ///< ��������� ��������� ������� �����.
);

/// ����� �������� ����������� ����� � �������������� �������� ������.
struct EXPORT pitch_templates_class: mask_templates {
	/// �����������. �������� ������� pitch_templates_generate
	pitch_templates_class(
		const freq_scale& sc, 
		const mask_parameters& pm, 
		const pitch_parameters& p) 
	{
		if(!pitch_templates_generate(sc, *this, pm, p)) {
			throw "pitch templates generation error";
		}
	}
	/// ����������. ������� ���������� ������.
	~pitch_templates_class() { spl_free(T); }
};

/// ����� � ������ ����� �� �������. 
/// �� ����� -- ����� �����, �� ������ - ������ ������� ��� ��� ������� �������.
/// ���������� ���������� 
size_t EXPORT mask_template_search(
	IN mask_templates& tpl,       ///< ������� �����, �� ������� �������������� �����
	io::istream<mask_t>& in_str,  ///< ������� ����� (�����)
	io::ostream<short>& out_str,  ///< �������� ����� (������ ��������� �������)
	IN int max_diff               ///< ������������ ������� ����� �� �������
);

size_t EXPORT mask_template_fast_search(
	IN mask_templates& tpl,       ///< ������� �����, �� ������� �������������� �����
	io::istream<mask_t>& in_str,  ///< ������� ����� (�����)
	io::ostream<short>& out_str,  ///< �������� ����� (������ ��������� �������)
	IN int max_diff               ///< ������������ ������� ����� �� �������
);

/// ����������� �������� ������������� ���������� ����� ��������������� ����� �� �������.
const int DEFAULT_PITCH_MAX_DIFF = 6;

/// ������� ����������� ������� ��������� ����.
size_t EXPORT pitch_track(
	IN freq_scale& sc,           ///< ����� ������
	io::istream<mask_t>& in_str, ///< ������� ����� (�����)
	io::ostream<short>& out_str, ///< �������� ����� (������ ������� ���)
	IN mask_parameters& pm,      ///< ��������� ��������� ����������� �������
	IN pitch_parameters& p       ///< ��������� ��������� ������� �����
);

/// ������� ����������� ������� �� �������� ����������������
size_t EXPORT vocal_segment(
	io::istream<short>& in_str,  ///< ������� ����� (������ ������� ���)
	io::ostream<short>& out_str, ///< �������� ����� (������ ��������-������ ���������)
	IN freq_t F,                 ///< ������� ������������� �������
	IN vocal_parameters& p       ///< ��������� �����������
);

} // namespace spl

#endif//_SPL_VOCAL_