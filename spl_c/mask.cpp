///
/// \file  mask.cpp
/// \brief ������� ��� ������� ������������� ���������� (\ref mask).
///
/// ������ �������� ���������� ��� ��������� �������� ���������� (������� mask_filters_generate())
///  � ���������� ��� ���������� (�������) - ������� mask_stream(), mask_memory(), mask_binary_file().
///

#include "mask.h"
#include "matrix.h"
#include "const.h"
#include "scale.h"
#include "conv.h"
#include "mask_helper.h"
#include <math.h>
#include <stdio.h>

#include <algorithm>

#include "iomem.h"
#include "iofile.h"
#include "iobit.h"
using io::istream;
using io::ostream;
using io::imstream;
using io::omstream;
using io::ifstream;
using io::ofstream;

namespace spl {

/// ������� ������� ���� �������� ����������.
/// ����������� ��� ���������������� ������� - ����� ���������� ������.

static int mask_window_size(const freq_scale& s, scale_form form, const mask_parameters& p) {
	// �������� ���������� �� ����
	int Ws = 0;

	double Nmin, Nmax;
	for(int k = 0; k < s.K; k++) {
		double std = mask_std(s.Fr[k], p.delta);
		double dF = gauss_border(std, p.ksi);
		if(p.border_effect) { 
			Nmin = scale_interp_f2k(s, form, s.Fr[k] - dF); // ���� ������� ������ �����������, 
			Nmax = scale_interp_f2k(s, form, s.Fr[k] + dF);	// �� ������������ ������������� ������
		} else {
			Nmin = scale_find_freq(s, s.Fr[k] - dF); // ���� ������� ������ �� �����������, 
			Nmax = scale_find_freq(s, s.Fr[k] + dF); // �� ������������ ������� ����� �� ��������
		}
		if(Nmax - Nmin + 1 > Ws) Ws = (int) ceil(Nmax - Nmin + 1);
	}
	return Ws % 2 ? Ws : Ws + 1; // ������ ���� ����������� �������� ����� - ����� ��� �������� [-Ws;+Ws]
}

/// ������� ��������� ����������� �������.
static void mask_win(
	const freq_scale& s, // ������ ����������� ������ ��������
	double *H,           // �������� ������ - ������������ ����������� �������
	int k, 
	int Ws, 
	enum scale_form form, 
	const mask_parameters& p
) {
	int K = s.K;
	freq_t *Fr = s.Fr;
	double *H2 = H + Ws;

	// ����������� ���������� (���) ������� ������
	double std = mask_std(Fr[k], p.delta);
	// ��������� ������������
	for(int m = -Ws; m <= Ws; m++) {
		double Hkm;
		if(k + m >= 0 && k + m < K) { // ����� k + m ���� � ����� ������?
			Hkm = gauss_win(Fr[k+m], Fr[k], std); // �� - ����� �� ����� ������
		} else if(p.border_effect) { // ���� ����������� ������� ������, 
			double Fm = scale_interp_k2f(s, form, k+m); // �� ������ �������������
			Hkm = gauss_win(Fm, Fr[k], std);
		} else { // ���� ������� ������ �� �����������,
			Hkm = 0; // �� ������ ������ 0 � ����������� �������
		}
		H2[m] = Hkm;
	}

	// ������������
	double sum = 0;
	for(int m = 2*Ws; m >= 0; m--) // ���� �� 0 �� 2*Ws ������������, �.�. ���-�� ��������� = 2*Ws+1
		sum += H[m];
	sum *= p.rho;
	for(int m = 2*Ws; m >= 0; m--)
		H[m] /= sum;
}


/// 
/// ���������� ������������� ������������� ���������� �� ����� ����������� ������ ��������.
/// � ��������� �������� \a f ��� ������ ���� �������� ����� (Ws*K) ��� ��������� ������������� ����������.
///
/// ��������� ������������� ������������ ������������� ���������� � ������ ��� ��� ����� �������� �������.
/// ��������� ������������ ������������� ���������� �������������� 
///  ��� ������� ������ ����� ������ \a s � ������ ������ �������� �������, 
///  �� ���������� ������������ �� ������� (������ � ���������) ������� ����� \a s 
///  �� �������������� ���������. 
/// ��� ������� ����� ������ ������������ ������������ ����� ������ \a s � ������� ������� scale_interp_k2f().
/// ��� ���� ����� ������������ ���� ��������, �������������� ����������� 
///  ����� ����� (\ref scale_form) � ������� ������� scale_form_estimate().
/// ������� ��� ����������� ������������ ������������� ���������� � ������ �������� �������, 
///  ���������� ������������ ����� ������ ��������� ����� (��. \ref scale_form).
///
/// ��� ������������� ����� ������ ����������� ����� 
///  ����� ������������ ������ ������������ ��� ����� �������� �������.
///
/// ���� �������� ������� ������������ ���������� \a p.border_effect.
///

bool mask_filters_generate(
	const freq_scale& s, 
	mask_filters& f, 
	const mask_parameters& p) 
{
	if(f.K != s.K) return false;

	// ������ ���� ������ ���� ��������
	// �� ���� ������� �������� - �� �������� �� -Ws/2 �� Ws/2
	int Ws = f.Ws/2; // �������� ����� ���� - ��� ��������
	if(f.Ws != 2*Ws+1) return false;

	// ��������� ������������
	for(int k = 0; k < f.K; k++) {
		mask_win(s, f.H + k * f.Ws, k, Ws, f.sc_form, p);
	}

	return true;
}

///
/// ���������� ������������� ������������� ����������, 
///  ���������������� ��� �������� ����������.
///
/// ������� ���������� ������������� ���������� ��������, 
///  ���� ����� ������ ����� ����� scale_model (��. \ref scale_form).
/// � ���� ������ ������������ ������������� ���������� �� ���� ��������� ������� ���������, 
///  ��� ��������� ������ ������������� ���������� � �������� ���������� � 
///  ��������� �� �� ������� � �������.
///
/// � ������ ������� ������������ ��������������� ������������� � �������� ����������:
///  ����������� ����� (������ �) �� ����������� �������.
/// 
/// � ��������� �������� \a f ��� ������ ���� �������� ����� (2*CONV_WIN_SIZ) 
///  ��� ��������� ������������� ����������.
/// 

bool mask_filters_generate_fast(
	const freq_scale& s, 
	mask_filters& f, 
	const mask_parameters& p) 
{
	if(f.K != s.K) return false;
	
	// ������ ���� ������ ���� ��������
	// �� ���� ������� �������� - �� �������� �� -Ws/2 �� Ws/2
	int Ws = f.Ws/2; // �������� ����� ���� - ��� ��������
	if(f.Ws != 2*Ws+1) return false;

	double H[CONV_WIN_SIZ];
	// ��������� ����������� ������� ��� k = K/2
	// ����� ����������� k �� ����� ���� �������
	mask_win(s, H, s.K/2, Ws, scale_model, p); // ��������� ������ f.Ws ����
	// ������ ����������� - �.�. ������� �������� ��� ��� ��� ;)
	std::reverse(H, H + f.Ws);
	// ��������� ������� ������
	std::fill(H + f.Ws, H + CONV_WIN_SIZ, 0.0);
	// �������������� ������� A
	cconv_calc_A(H, f.H, f.H + CONV_WIN_SIZ);
	// ���������� ������� �
	cconv_normalize(f.H, 2 * CONV_WIN_SIZ);

	return true;
}


mask_filters_class::mask_filters_class(
	const freq_scale& s, 
	const mask_parameters& p)
{
	K = s.K;
	sc_form = scale_form_estimate(s);
	Ws = mask_window_size(s, sc_form, p);

	allow_fast = (sc_form == scale_model && p.allow_fast);

	// ��� ��������� ����� - ������� ���������� - ������ ������
	size_t N = allow_fast ? 2 * CONV_WIN_SIZ : K * Ws;
	H = conv_alloc<double>(N);
	if(H == 0) 
		throw "Can't allocate memory for mask filters coefficients";

	typedef bool (*generate_func)(const freq_scale&, mask_filters&, const mask_parameters&);
	// ��� ��������� ����� - ������� ����������
	generate_func generate = allow_fast ? mask_filters_generate_fast : mask_filters_generate;

	if(!generate(s, *this, p)) 
		throw "Error while generating mask filters";
}


///
/// ��������� ������������� ����������.
/// ���������� ���������� ���������� �� ����� ���������.
///
/// �������� �������� ����� �������� ��� ���������� ������������� ����������, 
///  �������������� ���� � ����� ��� ����������� ������.
/// ����� ��������� ���������� �������� ���������� ���������� ������� 
///  (������, �����, ��������, ������� � �.�.),
///  �������������� ������� ��������� �������� � ������� 
///  � �������� ������ �������.
/// ��� ������������ �������� ����������, ��� ����������� �� ������������ ����� �������.
/// 
/// ������� �������� �� �������� ���������, ��� �����������.
///

size_t mask_stream(
	IN mask_filters& f,          ///< ������������ ����������
	istream<spectrum_t>& in_str, ///< ������� �����  (������)
	ostream<mask_t>& out_str     ///< �������� ����� (�����)
) {
	int K = f.K;
	int Ws = f.Ws/2;

	spectrum_t *spec_wide1 = spl_alloc<spectrum_t>(Ws + K + Ws);
	spectrum_t *spec_input = spec_wide1 + Ws;
	spectrum_t *spec_wide2 = spec_input + K;
	if(!spec_wide1) return 0;

	size_t written = 0;

	// �������� ���� ����������
	while(in_str.read(spec_input, K) == K && !out_str.eos()) {

		// ���������� ������� ������ �������
		std::fill(spec_wide1, spec_wide1 + Ws, spec_input[0]);
		std::fill(spec_wide2, spec_wide2 + Ws, spec_input[K-1]);

		// ���� ������� ����������:
		double *H = f.H;
		for(int k = 0; k < K; k++) {
			double sum = 0;
			for(int m = -Ws; m <= Ws; m++) {
				sum += spec_input[k+m] * (*H++);
			}
			written += out_str.put(spec_input[k] > sum);
		}
	}

	spl_free(spec_wide1);

	return written;

}


///
/// ��������� ������������� ���������� - ������� �������.
/// ���������� ������ ��� ��������� ����� ������.
///

size_t mask_stream_fast(
	IN mask_filters& f,          ///< ������������ ����������
	istream<spectrum_t>& in_str, ///< ������� �����  (������)
	ostream<mask_t>& out_str     ///< �������� ����� (�����)
) {
	int K = f.K;
	int Ws = f.Ws/2;
	int Os = CONV_WIN_SIZ - 2*Ws;

	spectrum_t array[CONV_WIN_SIZ * 6 + 2];
	// ��� ������ ��� ������ �������
	spectrum_t *input_buf1 = SPL_MEMORY_ALIGN(array);
	spectrum_t *input_buf2 = input_buf1 + CONV_WIN_SIZ;
	// ������ ��������� ������ - B, C, abcr, abci
	spectrum_t *tmp_buf1 = input_buf2 + CONV_WIN_SIZ;
	spectrum_t *tmp_buf2 = tmp_buf1 + CONV_WIN_SIZ;
	spectrum_t *tmp_buf3 = tmp_buf2 + CONV_WIN_SIZ;
	spectrum_t *tmp_buf4 = tmp_buf3 + CONV_WIN_SIZ;
	// ���� �������� �����
	mask_t out_buf[CONV_WIN_SIZ * 2];

	// i/o wrappers ��� ����������/������� ����� ������:
	istream_spectrum_extend in_wrap(in_str, K, Ws);
	ostream_mask_reduce out_wrap(out_str, K, Ws);

	// �������� ��������
	size_t written = 0;

	// ��������� ����� ������� ������, ����� ����� ������ ����������� � ������ �������
	std::fill(input_buf2 + Os, input_buf2 + Os + Ws, 0);
	in_wrap.read(input_buf2 + Os + Ws, Ws);

	size_t N1 = 0, N2 = 0;

	// �������� ���� ����������
	while(!in_wrap.eos() && !out_wrap.eos() || N2 > Os - Ws) { // ��������� ���� �������������� ������ �� ��������� ��������

		// �������� �� ����� ������� ������ � ������ �������
		std::copy(input_buf2 + Os, input_buf2 + CONV_WIN_SIZ, input_buf1);
		// ������ ������ ����� 
		N1 = in_wrap.read(input_buf1 + 2*Ws, Os);
		// ���� ����������� ������, ��� ����� - ������� ��������� ������
		if(N1 < Os) {
			std::fill(input_buf1 + 2*Ws + N1, input_buf1 + CONV_WIN_SIZ, 0);
		}

		// �������� �� ����� ������� ������ � ������ �������
		std::copy(input_buf1 + Os, input_buf1 + CONV_WIN_SIZ, input_buf2);
		// ������ ������ �����
		N2 = in_wrap.read(input_buf2 + 2*Ws, Os);
		// ���� ����������� ������, ��� ����� - ������� ��������� ������
		if(N2 < Os) {
			std::fill(input_buf2 + 2*Ws + N2, input_buf2 + CONV_WIN_SIZ, 0);
		}

		// ������� B, C
		cconv_calc_BC(input_buf1, input_buf2, tmp_buf1, tmp_buf2);

		// �������
		cconv(f.H, f.H + CONV_WIN_SIZ, tmp_buf1, tmp_buf2, tmp_buf3, tmp_buf4);

		int j = 0;
		// ��������� ��������� ���������� ��� ����� �������:
		for(size_t i = Ws; i < Ws + N1; i++) {
			out_buf[j++] = (input_buf1[i] > tmp_buf3[i + Ws]);
		}
		for(size_t i = Ws; i < Ws + N2; i++) {
			out_buf[j++] = (input_buf2[i] > tmp_buf4[i + Ws]);
		}

		// ������� ���������
		written += out_wrap.write(out_buf, j);
	}
	return written;
}


/// ���������� ������ ������� �� ����� ������.

size_t mask_stream(
	IN freq_scale& sc,           ///< ����� ����������� ������
	istream<spectrum_t>& in_str, ///< ������� �����  (������)
	ostream<mask_t>& out_str,    ///< �������� ����� (�����)
	IN mask_parameters& p        ///< ��������� ����������
) {
	try {
		// �������
		mask_filters_class f(sc, p);

		// ���������� ����������
		// � ����������� �� ����� ����� �������� ������� ��� ������� ������ ����������
		return f.mask_stream(in_str, out_str);

	} catch(const char *e) {
		fprintf(stderr, "%s", e);
		return 0;
	}
}

/// ���������� ������ ������� �� ������������ ����������� ��������� ����������� �������.

size_t mask_stream(
	IN freq_scale& sc,           ///< ����� ����������� ������
	istream<spectrum_t>& in_str, ///< ������� �����  (������)
	ostream<mask_t>& out_str,    ///< �������� ����� (�����)
	IN double ksi                ///< ��������, ������������ �������� ���������� (������ ���� ��������)
) {
	mask_parameters p = DEFAULT_MASK_PARAMETERS;
	p.ksi = ksi;
	return mask_stream(sc, in_str, out_str, p);
}

/// ���������� � ������

size_t mask_memory(
	IN freq_scale& sc, 
	IN spectrum& sp, 
	OUT mask& m, 
	IN mask_parameters& p
) {
	// ��������� ������ ���������� �������
	if(sp.K != m.K || sp.N != m.N) return 0;
	int NK = m.N * m.K;

	// �������������� ������� ����� �� ������
	imstream<spectrum_t> in_str(sp.Y, NK);

	// �������������� �������� ����� � ������
	omstream<mask_t> out_str(m.Z, NK);

	// ���������� ����������
	return mask_stream(sc, in_str, out_str, p);
}

/// ���������� � ������ �� ������������ �����������.

size_t mask_memory(
	IN freq_scale& sc,    ///< ����� ����������� ������ �������
	IN spectrum& sp,      ///< ������ (�������)
	OUT mask& m,          ///< ����� (��������)
	IN double ksi         ///< ��������� ����������
) {
	mask_parameters p = DEFAULT_MASK_PARAMETERS;
	p.ksi = ksi;
	return mask_memory(sc, sp, m, p);
}

/// ���������� ��������� �����.

size_t mask_binary_file(
	IN freq_scale& sc, 
	IN char *in_file, 
	IN char *out_file, 
	IN bool out_bits, 
	IN mask_parameters& p
) {
	// �������������� ������� ����� �� �����
	ifstream<spectrum_t> in_str(in_file);

	// �������� ��������
	size_t r;

	if(out_bits) { // ���� ����
		// �������������� �������� ����� � ����
		ofstream<unsigned char> out_str(out_file);
		io::obitwrap8 out_bit_str(out_str);
		// ���������� ����������
		r = mask_stream(sc, in_str, out_bit_str, p);
		out_bit_str.flush(); // ������ ���������� ����� (�����)
	} else { // ���� �����
		// �������������� �������� ����� � ����
		ofstream<mask_t> out_str(out_file);
		// ���������� ����������
		r = mask_stream(sc, in_str, out_str, p);
	}
	return r;
}

/// ���������� ��������� ����� �� ������������ �����������.

size_t mask_binary_file(
	IN freq_scale& sc,
	IN char *in_file,
	IN char *out_file,
	IN bool out_bits,
	IN double ksi
) {
	mask_parameters p = DEFAULT_MASK_PARAMETERS;
	p.ksi = ksi;
	return mask_binary_file(sc, in_file, out_file, out_bits, p);
}

} // namespace spl
