#ifndef _SPL_MASK_
#define _SPL_MASK_

///
/// \file  mask.h
/// \brief ���� ������, ��������� � �������, ����������� � ����������.
///

#include "common.h"
#include "spl_types.h"
#include "io.h"

namespace spl {

/// ���������, ������������ ��� ��������� ������������� ���������� (\ref mask_filters).

struct mask_parameters {

	/// ��������, ������������ �������� ���������� (������ ���� ��������).
	double ksi;

	/// �����������, ������������ ������ ����������� �������.
	double delta;

	/// �����������, ������������ ��� ����������� �������.
	double rho;

	/// ������� - ��������� ������� ������ ��� ��� (��. mask_filters_generate()).
	bool border_effect;

	/// ������� - ��������� ������� ������ ��� ���.
	bool allow_fast;

};

/// ����������� �������� ���������� ��������� ����������� �������.
const mask_parameters DEFAULT_MASK_PARAMETERS = { 0.001, 1, 1, true, true };

///
/// ���������� ������ ������� �� ����� ������.
/// ���������� ���������� ���������� ��������� �����.
/// 

size_t EXPORT mask_stream(
	IN freq_scale& sc,               ///< ����� ����������� ������
	io::istream<spectrum_t>& in_str, ///< ������� �����  (������)
	io::ostream<mask_t>& out_str,    ///< �������� ����� (�����)
	IN mask_parameters& p            ///< ��������� ����������
);

///
/// ���������� � ������
/// � ���� \a m.Y ������ ������ ���� �������� ����� �������, 
///  � ����� ����������� ���� \a m.K � \a m.N (��. \ref mask).
/// ������� ��������� ������������ ��������� ������ �� ���� �����.
///

size_t EXPORT mask_memory(
	IN freq_scale& sc,    ///< ����� ����������� ������ �������
	IN spectrum& sp,      ///< ������ (�������)
	OUT mask& m,          ///< ����� (��������)
	IN mask_parameters& p ///< ��������� ����������
);

/// ���������� � ������ �� ������������ �����������.

size_t EXPORT mask_memory(
	IN freq_scale& sc,    ///< ����� ����������� ������ �������
	IN spectrum& sp,      ///< ������ (�������)
	OUT mask& m,          ///< ����� (��������)
	IN double ksi         ///< ��������� ����������
);

///
/// ���������� ��������� �����.
///

size_t EXPORT mask_binary_file(
	IN freq_scale& sc,    ///< ����� ����������� ������ �������.
	IN char *in_file,     ///< ������� ���� (������)
	IN char *out_file,    ///< �������� ���� (�����)
	IN bool out_bits,     ///< �������� ������ (��� �������)
	IN mask_parameters& p ///< ��������� ����������
);

/// ���������� ��������� ����� �� ������������ �����������.

size_t EXPORT mask_binary_file(
	IN freq_scale& sc,    ///< ����� ����������� ������ �������.
	IN char *in_file,     ///< ������� ���� (������)
	IN char *out_file,    ///< �������� ���� (�����)
	IN bool out_bits,     ///< �������� ������ (��� �������)
	IN double ksi         ///< ��������� ����������
);

} // namespace spl

#endif//_SPL_MASK_
