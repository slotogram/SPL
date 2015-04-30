#ifndef _SPL_SPECTRUM_
#define _SPL_SPECTRUM_

///
/// \file  spectrum.h
/// \brief ���� ������, ��������� � �������, ����������� � ��������.
///

#include "common.h"
#include "spl_types.h"
#include "io.h"

namespace spl {

///
/// ���������� ������ ������� �� ����� ������.
/// ���������� ���������� ���������� ��������� �������.
///

size_t EXPORT filter_stream(
	IN freq_scale& sc,                ///< ����� ������
	io::istream<signal_t>& in_str,    ///< ������� ����� (������)
	io::ostream<spectrum_t>& out_str, ///< �������� ����� (������)
	IN freq_t F,                      ///< ������� �������������
	IN double ksi                     ///< ��������, ������������ �������� ���������� (������ ����)
);

///
/// ���������� � ����������� ������.
/// ���������� ���������� ���������� ��������� �������.
///
/// � ���� \a sp.Y ������ ������ ���� �������� ����� �������, 
///  � ����� ����������� ���� \a sp.K � \a sp.N (��. \ref spectrum).
/// ������� ��������� ������������ ��������� ������ �� ���� �����.
///

size_t EXPORT filter_memory(
	IN freq_scale& sc, ///< ����� ����������� ������ ��������
	IN signal& s,      ///< ������� ������ (� ������)
	OUT spectrum& sp,  ///< �������� ������ (� ������)
	IN double ksi      ///< ��������, ������������ �������� ���������� (������ ����)
);

///
/// ���������� ��������� ����� ��� �������.
/// ���������� ���������� ���������� ��������� �������.
///

size_t EXPORT filter_binary_file(
	IN freq_scale& sc, ///< ����� ����������� ������ ��������
	IN char *in_file,  ///< ������� ������ (��� �����)
	IN char *out_file, ///< �������� ������ (��� �����)
	IN freq_t F,       ///< ������� �������������
	IN double ksi      ///< ��������, ������������ �������� ���������� (������ ����)
);

/// 
/// ���������� WAV-�����.
/// ���������� ���������� ���������� ��������� �������.
/// 

size_t EXPORT filter_wave_file(
	IN freq_scale& sc, ///< ����� ����������� ������ ��������
	IN char *in_file,  ///< ������� ������ (��� �����)
	IN char *out_file, ///< �������� ������ (��� �����)
	IN double ksi      ///< ��������, ������������ �������� ���������� (������ ����)
);


/// 
/// ���������� WAV-����� � ������.
/// ���������� ���������� ���������� ��������� �������.
/// 
size_t EXPORT filter_wmemory(
	IN freq_scale& sc, ///< ����� ����������� ������ ��������
	IN char *in_file,  ///< ������� ������ (��� �����)
	OUT spectrum& sp,  ///< �������� ������ (� ������)
	IN double ksi,     ///< ��������, ������������ �������� ���������� (������ ����)
	IN int k			///< ������
	);
} // namespace spl 

#endif//_SPL_SPECTRUM_
