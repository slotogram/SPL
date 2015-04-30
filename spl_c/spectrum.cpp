///
/// \file  spectrum.cpp
/// \brief ������ ��������� ������� ��������� ������� ("����������").
///
/// ��������� ������� ����������� ����� ������� ������������� ��������. 
/// ������ �������� ���������� ��� ��������� ���������� ������ �������� � ���������� ��� ���������� (�������).
///

#include "matrix.h"
#include "const.h"
#include "spl_types.h"
#include "signal.h"
#include "spectrum.h"
#include "spectrum_helper.h"
#include "conv.h"

#include "iomem.h"
#include "iofile.h"
#include "iowave.h"
using io::istream;
using io::ostream;
using io::imstream;
using io::omstream;
using io::ifstream;
using io::ofstream;

#define _USE_MATH_DEFINES // M_PI, etc
#include <math.h>
#include <algorithm>


namespace spl {

/// 
/// ��������� ������������� ���������� �� ����� ����������� ������ ��������.
/// � ��������� �������� \a f ��� ������ ���� �������� ����� ��� ��������� ������������� ����������.
/// 

bool spectrum_filters_generate(
	IN freq_scale& s,        ///< ����� ����������� ������ ��������
	OUT spectrum_filters& f, ///< ������������ ����������
	IN filter_parameters& p  ///< ��������� ��������� ��������
) {
	int n, j, k; // ���������
	int K = s.K;
	freq_t F = p.F;

	double array[CONV_WIN_SIZ * 2 + 2];
	double *Hc = SPL_MEMORY_ALIGN(array);
	double *Hs = Hc + CONV_WIN_SIZ;

	if(s.K != f.K) return false;

	// ��������� Ws - ������ ��������� ���� ��������
	// ���������� ��� ������������ �� ����������� ��� ������� ������
	int Ws = 0;
	for(k = 0; k < K; k++) {
		const double std = filter_std(s.Fr[k], F);
		int Wsk = (int) gauss_border(std, p.ksi);
		if(Wsk > Ws) Ws = Wsk; // �� ����� ������ ��� ���������, ��� ���� �� �����
	}
	f.Ws = 2*Ws+1;
	
	// ������� - ��� �������� ������� � ������������� ����������
	Matrix<double, 3> H = matrix_ptr(f.H, 2, f.K, CONV_WIN_SIZ);

	// ��������� ���������� ������������ �������
	// ��� ������� ���� Ws
	for(k = 0; k < K; k++) {
		const double std = filter_std(s.Fr[k], F);
		const double Wf = 2 * M_PI * s.Fr[k] / F;
		for(n = -Ws, j = 0; n <= Ws; n++, j++) {
			double norm = gauss_win<double>(n, 0.0, std); // ���� ������ ��������� ���������
			Hc[j] = norm * cos(Wf * n); // �������������� �����
			Hs[j] = norm * sin(Wf * n); // ������ �����
		}
		// ��������� ���������� ������������ ������
		std::fill(Hc + j, Hc + CONV_WIN_SIZ, 0);
		std::fill(Hs + j, Hs + CONV_WIN_SIZ, 0);
		// ��������������� ���������� ������� BC
		cconv_calc_BC(Hc, Hs, &H(0,k,0), &H(1,k,0));
	}
	// ������������ ������������� ����������:
	cconv_normalize(H, H.size());
	return true;
}

spectrum_filters_class::spectrum_filters_class(
	const freq_scale& s, 
	const filter_parameters& p)
{
	K = s.K;
	H = conv_alloc<double>(K * 2 * CONV_WIN_SIZ);
	if(H == 0) 
		throw "Can't allocate memory for spectrum filters coefficients";
	if(!spectrum_filters_generate(s, *this, p)) 
		throw "Error while generating spectrum filters";
}

spectrum_filters_class::~spectrum_filters_class() {
	conv_free(H);
}

bool spectrum_filters_class::save(const char *file) { 
	return io::array_to_file(H, K * CONV_WIN_SIZ * 2, file); 
}

///
/// ���������� ������ � �������� ������� ��������.
/// ���������� ���������� ���������� �� ����� ���������.
/// 
/// �������� �������� ����� �������� ��� ����������, 
///  �������������� ���� � �����, ��� ����������� ������.
/// ����� ��������� ���������� �������� ���������� ���������� ������� 
///  (������, �����, ��������, ������� � �.�.),
///  �������������� ������� ��������� �������� � ������� 
///  � �������� ������ �������.
/// ��� ������������ �������� ����������, ��� ����������� �� ������������ ����� �������.
///
/// ������ ������� ���������� ����������� ���������� ������� ����� FFT.
///

size_t filter_stream(
	IN spectrum_filters& f,       ///< ������������ ����������
	istream<signal_t>& in_str,    ///< ������� �����  (������)
	ostream<spectrum_t>& out_str  ///< �������� ����� (������)
) {
	int K = f.K;
	size_t Ws = f.Ws - 1; // ����� ����� �� 1 ������, ��� ���� - ��������� �� ��������
	size_t Os = CONV_WIN_SIZ - Ws;
	spectrum_t *out_buf = NULL;

	// ������������ ���������� �������� � ������ �������
	iwstream_extend<signal_t> in_wrap(in_str, Ws/2);

	// ������� �������� - �������� ��������
	size_t written = 0;

	double array[5*CONV_WIN_SIZ + 2];

	// ����� �������� ������� - ������� �� ���� ������:
	// CONV_WIN_SIZ = Ws + Os, 
	// ��� Ws - Window size - ������ ��������� ���� �������, 
	//     CONV_WIN_SIZ - ������ ����������� ����������� �������
	//     Os - Output size - ������ ��������� ������ �������
	// ����� ������������ ��� ������� ����� �������
	double *conv_in_buf = SPL_MEMORY_ALIGN(array);

	// �������� ����� �������
	// ����� ����� �� ��������� ��� � ������� ����� (2*Ws+1) + Os
	// ������ �������� ����� - ��������� Os ��������� - ���� �� �����
	double *tmp_buf1 = conv_in_buf + CONV_WIN_SIZ;
	double *tmp_buf2 = tmp_buf1 + CONV_WIN_SIZ;
	double *tmp_buf3 = tmp_buf2 + CONV_WIN_SIZ;
	double *tmp_buf4 = tmp_buf3 + CONV_WIN_SIZ;

	// ����� ��������� �������
	// ������� ������� K x Os - ��� ��������� ���������� �������
	// � Os x K - ��� ������ ������
	out_buf = spl_alloc<spectrum_t>(2 * K * Os);

	Matrix<spectrum_t,2> out_mtx1 = matrix_ptr(out_buf, K, Os);
	Matrix<spectrum_t,2> out_mtx2 = matrix_ptr(out_buf+K*Os, Os, K);

	// ������� - ��� �������� ������� � ������������� ����������
	Matrix<double, 3> H = matrix_ptr(f.H, 2, f.K, CONV_WIN_SIZ);

	if(!out_buf) 
		goto end; 

	//
	// ���������� �������� ������
	//

	// ������� ��������� Ws ��������� �������� ������
	std::fill(conv_in_buf+CONV_WIN_SIZ-Ws, conv_in_buf+CONV_WIN_SIZ-Ws/2, 0);

	// ������������ ���������� �������� � ������ �������
	in_wrap.read(conv_in_buf+CONV_WIN_SIZ-Ws/2, Ws/2);

	// �������� ���� ����������
	while(!in_wrap.eos() && !out_str.eos()) {

		// �������� ��������� Ws ��������� ������� � ������
		std::copy(conv_in_buf+CONV_WIN_SIZ-Ws, conv_in_buf+CONV_WIN_SIZ, conv_in_buf);

		// ������ Os ����� �������� �������
		// ���� ����� ����� �������������� ��������� ��� ������� ������
		// �������� rOs - real output size - ���������� ��������� ���������
		// � ����� ������ rOs == Os, ������� ����� ���� ������ � ����� �������
		size_t rOs = in_wrap.read(conv_in_buf + Ws, Os);

		// ���� ������� ������, ��� Os ���������,
		//  �� ����� ��������� ����� �����
		// ������� ��������� ������� �������, 
		// � �������� ������� �� ������� �������� ������� (K x rOs)
		if(rOs != Os) {
			std::fill(conv_in_buf + Ws + rOs, conv_in_buf + CONV_WIN_SIZ, 0);
			out_mtx1 = matrix_ptr(out_buf, K, rOs);
			out_mtx2 = matrix_ptr(out_buf+K*rOs, rOs, K);
		}

		spectrum_t *out_ptr = out_mtx1;

		cconv_calc_A(conv_in_buf, tmp_buf1, tmp_buf2);

		// ���� �� �������
		for(int k = 0; k < K; k++) {

			// �������
			cconv(tmp_buf1, tmp_buf2, &H(0,k,0), &H(1,k,0), tmp_buf3, tmp_buf4);

			// ���������� ������ ����������� �����
			complex_abs_split(rOs, tmp_buf3 + Ws, tmp_buf4 + Ws, tmp_buf4 + Ws);

			// ����� ����� � �������� ������� - ������ �� ��������� [Ws, Ws+rOs]
			out_ptr = std::copy(tmp_buf4 + Ws, tmp_buf4 + Ws + rOs, out_ptr);

		}

		// ������������� �������� �������
		// �� K x rOs � rOs x K
		transpose(out_mtx1, out_mtx2);

		// ������� ������� out_mtx2
		written += out_str.write(out_mtx2, out_mtx2.size());

	}

end:
	// ������� ���������� ������
	spl_free(out_buf);

	return written;

}

///
/// ���������� ������ ������� �� ����� ������.
/// 

size_t filter_stream(
	IN freq_scale& sc, ///< ����� ������ 
	istream<signal_t>& in_str, ///< ������� ����� (������)
	ostream<spectrum_t>& out_str, ///< �������� ����� (������)
	IN freq_t F, ///< ������� �������������
	IN double ksi) ///< ��������, ������������ �������� ���������� (������ ����)
{
	// ��������� �������� ��������.
	filter_parameters p;
	p.F = F;
	p.ksi = ksi;

	try {

		// �������
		spectrum_filters_class f(sc, p);

		// ���������� ����������
		return filter_stream(f, in_str, out_str);

	} catch(const char *e) {
		fprintf(stderr, "%s", e);
		return 0;
	}

}

///
/// ���������� � ����������� ������.
///

size_t filter_memory(
	IN freq_scale& sc, ///< ����� ����������� ������ ��������
	IN signal& s,      ///< ������� ������ (� ������)
	OUT spectrum& sp,  ///< �������� ������ (� ������)
	IN double ksi)     ///< ��������, ������������ �������� ���������� (������ ����)
{
	// ��������� ������ ���������� �������
	if(sp.K != sc.K || sp.N != s.N) return 0;

	// ������ ������� �������������
	sp.F = s.F;
	
	// �������������� ������� ����� �� ������
	imstream<signal_t> in_str(s.X, s.N);

	// �������������� �������� ����� � ������
	omstream<spectrum_t> out_str(sp.Y, sp.N * sp.K);

	// ���������� ����������
	return filter_stream(sc, in_str, out_str, s.F, ksi);

}

///
/// ���������� �� wav � ����������� ������.
///

size_t filter_wmemory(
	IN freq_scale& sc, ///< ����� ����������� ������ ��������
	IN char *in_file,  ///< ������� ������ (��� �����)
	OUT spectrum& sp,  ///< �������� ������ (� ������)
	IN double ksi,     ///< ��������, ������������ �������� ���������� (������ ����)
	IN int k)			///< ������
{
	
	// �������������� ������� ����� �� �����
	iwstream in_str(in_file);
	
	//��������� � ������
	
	spl::signal sig = in_str.getInMemory();

	//�������� ������ �� ������
	//spectrum* spec = new spectrum;
	sp.F = sig.F;
	sp.N = sig.N;
	sp.K = k;
	sp.Y = new spectrum_t[sp.N*k];


	// ���������� ����������
	size_t ret = filter_memory(sc,sig,sp,ksi);
	
	delete [] sig.X;
	//delete &sig;

	return ret;

}


///
/// ���������� ������� ��� ��������� �����.
///

size_t filter_binary_file(
	IN freq_scale& sc, ///< ����� ����������� ������ ��������
	IN char *in_file,  ///< ������� ������ (��� �����)
	IN char *out_file, ///< �������� ������ (��� �����)
	IN freq_t F,       ///< ������� �������������
	IN double ksi      ///< ��������, ������������ �������� ���������� (������ ����)
) {

	// �������������� ������� ����� �� �����
	ifstream<signal_t> in_str(in_file);

	// �������������� �������� ����� � ����
	ofstream<spectrum_t> out_str(out_file);

	// ���������� ����������
	return filter_stream(sc, in_str, out_str, F, ksi);

}

/// ���������� wave-�����.

size_t filter_wave_file(
	IN freq_scale& sc, ///< ����� ����������� ������ ��������
	IN char *in_file,  ///< ������� ������ (��� �����)
	IN char *out_file, ///< �������� ������ (��� �����)
	IN double ksi      ///< ��������, ������������ �������� ���������� (������ ����)
) {

	// �������������� ������� ����� �� �����
	iwstream in_str(in_file);

	// �������������� �������� ����� � ����
	ofstream<spectrum_t> out_str(out_file);

	// ���������� ����������
	return filter_stream(sc, in_str, out_str, in_str.freq(), ksi);

}


} // namespace spl 
