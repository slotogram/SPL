#ifdef _WIN32
# define NOMINMAX
# include <windows.h>
#endif
#include <time.h>


#include "const.h"
#include "iofile.h"
#include "iomem.h"
#include "iobit.h"
#include "iowave.h"
#include "iomic.h"
#include "scale.h"
#include "spectrum.h"
#include "spectrum_helper.h"
#include "mask.h"
#include "mask_helper.h"
#include "vocal.h"
#include "conv.h"
#include "fftw3.h"
#include <limits>
#include <stdio.h>
#include <math.h>

#include <algorithm>

using namespace spl;
using namespace io;
using std::auto_ptr;

namespace {

// time functions
#define SPL_TEST_USE_WIN32_TIME

#if defined(_WIN32) && defined(SPL_TEST_USE_WIN32_TIME)
static DWORD test_time = 0; 
long long get_time() { return GetTickCount(); }
void tic() { test_time = get_time(); }
DWORD toc() { return get_time() - test_time; }
#else
static long long test_time = 0; 
long long get_time() { time_t t; ctime(&t); return t; }
void tic() { test_time = get_time(); }
long long toc() { return get_time() - test_time; }
#endif

// исходные эталоны
const char *original_std = "original.txt";
const char *signal_wav_std = "testdata/signal.wav";
const char *numbers_std = "testdata/256.bin";

// производные эталоны
const char *signal_std = "testdata/signal.bin";
const char *scale_std = "scale-model-freq.bin";
const char *spectrum_std = "testdata/spectrum.bin";
const char *spectrum_filters_std = "testdata/spectrum-filters.bin";
const char *mask_byte_std = "testdata/sync-mask.bin";
const char *mask_bit_std = "testdata/sync-mask.bit";
const char *pitch_chan_std = "testdata/pitch-chan.bin";
const char *pitch_freq_std = "testdata/pitch-freq.bin";
const char *vocal_chan_std = "testdata/vocal-chan.bin";

// производные тестовые данные
const char *scale_test = "scale-model-freq.bin";
const char *spectrum_test = "spectrum.bin";
const char *spectrum_filters_test = "spectrum-filters.bin";
const char *mask_byte_test = "sync-mask.bin";
const char *mask_bit_test = "sync-mask.bit";
const char *mask_fast_test = "sync-mask.bin";
const char *mask_fast_bit_test = "sync-mask.bit";
const char *pitch_chan_test = "pitch-chan.bin";
const char *pitch_freq_test = "pitch-freq.bin";
const char *vocal_chan_test = "vocal-chan.bin";
const char *numbers_test = "256.bin";

// максимальные ошибки
const signal_t max_error_iwstream = 1E-10;
const freq_t max_error_scale_generate = 1E-10;
const freq_t max_error_scale_interp = 1E-10;
const spectrum_t max_error_filter = 1E+0;
const double max_error_fft_approx = 1E-10;
const double max_error_fftw_split = 5E-10;
const double max_error_cconv_naive = 1E-100;
const double max_error_cconv_fast = 1E-10;
const unsigned char max_error_mask = 0;
const double max_error_pitch_track = 0;
const double max_error_vocal_segment = 0;

// 
// функции для печати
// 

void print_time() { 
	printf("%i ms elapsed\n", toc()); 
}

void print_error(double e) {
	printf("error = %lg\n", e);
}

void print_failed_at(int i) {
	printf(" at %i\n", i);
}

void print_status(bool status) {
	printf(status ? "OK\n" : "FAILED\n");
}


//
// Функции для работы со структурой original - содержащей исходные данные по тестам
//

struct Original {
	double alpha, beta, B, Xm, Fsl, Fsh, Fn, Fv, K, delta, rho, ksi, Fon, Fov, Ng, F, border_effect, minV, minNV;
} original;

struct original_field {
	const char *name;
	double Original::*field;
	double default_value;
} original_fields[] = {
	{ "alpha", &Original::alpha, f_alpha },
	{ "beta", &Original::beta, f_beta },
	{ "B", &Original::B, f_B },
	{ "Xm", &Original::Xm, e_Xm },
	{ "Fsl", &Original::Fsl, e_Flow },
	{ "Fsh", &Original::Fsh, e_Fhigh },
	{ "Fn", &Original::Fn, 0 },
	{ "Fv", &Original::Fv, 0 },
	{ "K", &Original::K, 0 },
	{ "delta", &Original::delta, 0 },
	{ "rho", &Original::rho, 0 },
	{ "ksi", &Original::ksi, 0 },
	{ "Fon", &Original::Fon, 0 },
	{ "Fov", &Original::Fov, 0 },
	{ "Ng", &Original::Ng, 0 },
	{ "F", &Original::F, 0 },
	{ "border_effect", &Original::border_effect, 0 },
	{ "minV", &Original::minV, 0 },
	{ "minNV", &Original::minNV, 0 },
};

void read_original(const char *filename, Original *o) {
	char name[100];
	double value;
	FILE *f = fopen(filename, "r");
	if(f == NULL) {
		fprintf(stderr, "can't read original data from file: %s\n", filename);
		return;
	}
	while(!feof(f)) {
		if(fscanf(f, "%s%lf", name, &value) != 2) break;
		size_t n = strlen(name) - 1;
		if(name[n] == ':') name[n] = 0;
		for(int i = 0; i < sizeof(original_fields)/sizeof(original_field); i++) {
			original_field& f = original_fields[i];
			if(0 == strcmp(name, f.name)) {
				o->*(f.field) = value;
				break;
			}
		}
	}
}

bool check_original(const Original *o) {
	for(int i = 0; i < sizeof(original_fields)/sizeof(original_field); i++) {
		original_field& f = original_fields[i];
		if(f.default_value 
		&& o->*(f.field) != f.default_value) 
		{
			fprintf(stderr, "original field '%s' is not of default value\n", f.name);
			return false;
		}
	}
	return true;
}


// функции для сравнения потоков

template<typename T>
double compare_streams(io::istream<T>& str1, io::istream<T>& str2) {
	double e = 0.0;
	T x1, x2;
	// use &, not && - for both streams to be fetched
	while(str1.get(x1) & str2.get(x2)) {
		e += fabs(double(x1) - double(x2));
	}
	if(str1.eos() != str2.eos()) {
		printf("compare_streams: warning: streams have different sizes\n");
	}
	return e;
}

template<typename T>
double compare_streams(const char *filename1, const char *filename2) {
	io::ifstream<T> str1(filename1), str2(filename2);
	return compare_streams(str1, str2);
}

////////////////////////////////////////////////////////////////
// ТЕСТЫ
////////////////////////////////////////////////////////////////

void test_iwstream() {
	signal_t e = compare_streams(ifstream<signal_t>(signal_std), iwstream(signal_wav_std));
	print_error(e);
	print_status(e <= max_error_iwstream);
}

void test_scale_generate() {

	// generate scale
	scale_parameters p;
	p.a = p.i = 0;
	p.b = p.j = int(original.K - 1);
	p.Fi = original.Fn;
	p.Fj = original.Fv;
	scale_class sc(p, scale_model);
	sc.save(scale_test);

	freq_t e = compare_streams<freq_t>(scale_std, scale_test);

	print_error(e);
	print_status(e <= max_error_scale_generate);
}

void test_scale_form_estimate()	{
	// generate scale
	scale_parameters p;
	p.a = p.i = 0; p.b = p.j = 255; p.Fi = 50; p.Fj = 4000;
	scale_class sc(p, scale_model);

	bool status = true;
	int i;
	enum scale_form form;
	enum scale_form forms[] = { scale_linear, scale_log, scale_mel, scale_bark, scale_model, scale_unknown };
	for(i = 0; forms[i] != scale_unknown; i++) {
		if(!scale_generate(sc, p, forms[i])) { 
			printf("scale_form_estimate: can't generate scale\n");
			continue;
		}
		form = scale_form_estimate(sc);
		if(forms[i] != form) { 
			print_failed_at(i);
			print_error(double(form));
			status = false;
		}
	}

	print_status(status);
}

void test_scale_interp() {

	// базовая шкала
	scale_parameters p;
	p.a = p.i = 0; p.b = p.j = 255; p.Fi = 50; p.Fj = 4000;
	scale_class sc(p, scale_model);

	// расширенная шкала
	const int dk = 50;
	p.a = p.i - dk;
	p.b = p.j + dk; 
	scale_class sc2(p, scale_model);

	bool status;

	// тест интерполяции scale_generate
	int m;
	freq_t e;
	status = true;
	for(m = p.i; m < p.j; m++) {
		e = fabs(sc.Fr[m] - sc2.Fr[m+dk]);
		if(e > max_error_scale_interp) {
			status = false;
			break;
		}
	}
	print_status(status);
	if(!status) {
		print_failed_at(m);
		print_error(e);
	}
		
	// тест интерполяции scale_interp_k2f
	freq_t fm;
	status = true;
	for(m = p.a; m <= p.b; m++) {
		fm = scale_interp_k2f(sc, scale_model, m);
		e = fabs(fm - sc2.Fr[m+dk]);
		if(e > max_error_scale_interp) {
			status = false;
			break;
		}
	}
	print_status(status);
	if(!status) {
		print_failed_at(m);
		print_error(e);
	}

	// тест интерполяции scale_interp_f2k
	double k;
	status = true;
	for(m = p.a; m <= p.b; m++) {
		k = scale_interp_f2k(sc, scale_model, sc2.Fr[m+dk]);
		e = fabs(k - m);
		if(e > max_error_scale_interp) {
			status = false;
			break;
		}
	}
	print_status(status);
	if(!status) {
		print_failed_at(m);
		print_error(e);
	}
}

void test_spectrum_filters() {

	// standard scale
	scale_class sc(scale_std);

	// filter parameters
	filter_parameters p;
	p.F = original.F;
	p.ksi = original.ksi;

	// generate and save scale
	spectrum_filters_class f(sc, p);
	f.save(spectrum_filters_test);
}

void test_filter_binary_file() {
	// standard scale
	scale_class sc(scale_std);

	tic();
	size_t written = filter_binary_file(sc, signal_std, spectrum_test, original.F, original.ksi);
	print_time();

	spectrum_t e = compare_streams<spectrum_t>(spectrum_std, spectrum_test);
	print_error(e);
	print_status(e <= max_error_filter);
}

void test_filter_wave_file() {
	// standard scale
	scale_class sc(scale_std);

	tic();
	size_t written = filter_wave_file(sc, signal_wav_std, spectrum_test, original.ksi);
	
	print_time();

	spectrum_t e = compare_streams<spectrum_t>(spectrum_std, spectrum_test);
	print_error(e);
	print_status(e <= max_error_filter);
}

void test_mask_extend_reduce() {
	int Ws = 12;
	int K = int(original.K);
	{
		io::ifstream<spectrum_t> in(spectrum_std);
		io::ofstream<spectrum_t> out(spectrum_test);
		spl::istream_block_extend<spectrum_t> in1(in, K, Ws);
		spl::ostream_block_reduce<spectrum_t> out1(out, K, Ws);

		tic();
		spectrum_t buf[1500];
		while(!in1.eos()) {
			size_t N = in1.read(buf, 1500);
			size_t M = out1.write(buf, N);
		}
		print_time();
	}

	spectrum_t e = compare_streams<spectrum_t>(spectrum_std, spectrum_test);
	print_error(e);
	print_status(e <= 0.0);
}

void test_mask_bytes(bool allow_fast) {
	// standard scale
	scale_class sc(scale_std);

	// mask parameters
	mask_parameters p;
	p.border_effect = original.border_effect != 0.0;
	p.ksi = original.ksi;
	p.rho = original.rho;
	p.delta = original.delta;
	p.allow_fast = allow_fast;

	// test byte version
	tic();
	size_t written = mask_binary_file(sc, spectrum_std, mask_byte_test, false, p);
	print_time();

	double std = filter_std(sc.Fr[0], original.F);
	int Ws = int(ceil(gauss_border(std, original.ksi)));

	double e = compare_streams<unsigned char>(mask_byte_std, mask_byte_test);
	print_error(e);
	print_status(e <= max_error_filter);
}

void test_mask_bits(bool allow_fast) {
	// standard scale
	scale_class sc(scale_std);

	// mask parameters
	mask_parameters p;
	p.border_effect = original.border_effect != 0.0;
	p.ksi = original.ksi;
	p.rho = original.rho;
	p.delta = original.delta;
	p.allow_fast = allow_fast;

	// test byte version
	tic();
	size_t written = mask_binary_file(sc, spectrum_std, mask_bit_test, true, p);
	print_time();

	double std = filter_std(sc.Fr[0], original.F);
	int Ws = int(ceil(gauss_border(std, original.ksi)));

	double e = compare_streams<unsigned char>(mask_bit_std, mask_bit_test);
	print_error(e);
	print_status(e <= max_error_filter);
}

void test_mask_bytes_naive() {
	test_mask_bytes(false);
}

void test_mask_bytes_fast() {
	test_mask_bytes(true);
}

void test_mask_bits_naive() {
	test_mask_bits(false);
}

void test_mask_bits_fast() {
	test_mask_bits(true);
}

void test_iobit() {

	typedef unsigned short ushort;
	typedef unsigned long ulong;
	typedef io::imstream<ushort> imstream16;
	typedef io::omstream<ushort> omstream16;
	using io::ibitwrap16;
	using io::obitwrap16;

	ushort x[4] = { 0x3210, 0x7654, 0xBA98, 0xFEDC };
	ushort y[4] = { 0, 0, 0, 0 };
	imstream16 input_byte(x, sizeof(x)/sizeof(ushort));
	omstream16 output_byte(y, sizeof(y)/sizeof(ushort));
	ibitwrap16 input_bit(input_byte);
	obitwrap16 output_bit(output_byte);

	printf("testing ibitwrap::get and obitwrap::put: ");
	for(int i = 0; i < 17; i++) {
		for(int j = 0; j < 4; j++) {
			bool b;
			if(input_bit.get(b)) { 
				printf("%i", b); 
				if(output_bit.put(b)) { printf("v"); }
				else { printf("x"); }
			}
			else { printf("X"); }
		}
		printf(" ");
	}
	printf("- %X %X %X %X ", y[0], y[1], y[2], y[3]);
	printf("- %s\n", y[0] == x[0] && y[1] == x[1] && y[2] == x[2] && y[3] == x[3] ? "OK" : "FAILED");

	input_byte.pos(0);
	output_byte.pos(0);
	y[0] = y[1] = y[2] = y[3] = 0;

	printf("testing ibitwrap::read && obitwrap::write: ");
	while(!input_bit.eos() && !output_bit.eos()) {
		bool b[19];
		size_t r, w;
		r = input_bit.read(b, 19); // read by 19 bits
		w = output_bit.write(b, 19); // write by 19 bits
		printf("%i-%i ", r, w);
	}
	printf("- %X %X %X %X ", y[0], y[1], y[2], y[3]);
	printf("- %s\n", y[0] == x[0] && y[1] == x[1] && y[2] == x[2] && y[3] == x[3] ? "OK" : "FAILED");

}

void test_iobuf() {
	typedef short elem_t;

	ifstream<elem_t> input(numbers_std);
	ofstream<elem_t> output(numbers_test);
	iobuf<elem_t> buffer;

	elem_t block[25];
	size_t size1 = 20, 
		size2 = 15,
		size3 = 25,
		read, written;
	while(!input.eos()) {
		read = input.read(block, size1);
		written = buffer.output().write(block, read);
		if(written != read) {
			printf("written != read: %d != %d\n", written, read);
			print_status(false);
			break;
		}
		read = buffer.input().read(block, size2);
		if(read != size2) {
			printf("should have read %d elements, but read %d elements\n", size2, read);
			print_status(false);
		}
		output.write(block, read);
	}
	buffer.output().close();

	// write the remainder
	while(!buffer.input().eos()) {
		read = buffer.input().read(block, size3);
		output.write(block, read);
	}

	input.close();
	output.close();

	double e = compare_streams<elem_t>(numbers_std, numbers_test);
	print_status(e == 0.0);

}

void test_imicstream() {

	typedef short elem_t;
	int sample_rate = 16000;

	iobuf<elem_t> pipe;
	istream<elem_t>& input = pipe.input();
	mic_writer<elem_t> writer(pipe.output(), sample_rate);
	ofstream<elem_t> output("microphone.bin");
	
	elem_t buffer[150];
	size_t bufsize = 150;
	writer.start();
	for(int i = 0; i < 150; i++) {
		size_t written = input.read(buffer, bufsize);
		if(written != bufsize) {
			printf("written %d, fullsize %d \n", written, bufsize);
			print_status(false);
			return;
		}
		output.write(buffer, written);
	}

	writer.stop();
	pipe.output().close();

	print_status(true);
}


void fill_test_array(fftw_complex *in, double *ri, double *ii) {
	for(int i = 0; i < CONV_WIN_SIZ; i++) {
		double I = double(i);
		in[i][0] = ri[i] = sin(I) * sin(2*I);
		in[i][1] = ii[i] = sin(I) * sin(2*I) + 5;
	}
}

double compare_fftw_arrays(int N, IN fftw_complex *out, IN double *ro, IN double *io) {
	double e = 0.0;
	for(int i = 0; i < N; i++) {
		e += fabs((out[i][0] - ro[i]));
		e += fabs((out[i][1] - io[i]));
	}
	return e;
}

void test_fftw_split() {
	double e;

	int K = 1;
	double *array = conv_alloc<double>(CONV_WIN_SIZ * 5 * 8 + 2);
	double *ri = SPL_MEMORY_ALIGN(array);
	double *ii = ri + CONV_WIN_SIZ * K;
	double *ro = ii + CONV_WIN_SIZ * K;
	double *io = ro + CONV_WIN_SIZ * K;
	fftw_complex *in = (fftw_complex *)(ro + CONV_WIN_SIZ * K);
	fftw_complex *out= in + CONV_WIN_SIZ * K;

	

	fftw_iodim dim, howmany;
	dim.n = CONV_WIN_SIZ;
	dim.is = dim.os = 1;
	fftw_plan plan1, plan2;
	howmany.n = K;
	howmany.is = howmany.os = CONV_WIN_SIZ;

	plan1 = fftw_plan_dft_1d(CONV_WIN_SIZ, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
	plan2 = fftw_plan_guru_split_dft(1, &dim, 0, 0, ri, ii, ro, io, FFTW_ESTIMATE);
	
	fill_test_array(in, ri, ii);
	for(int k = 0; k < K; k++) fftw_execute_dft(plan1, in + k * howmany.is, out + k * howmany.os);
	fftw_execute(plan2);
	fftw_destroy_plan(plan1);
	fftw_destroy_plan(plan2);

	e = compare_fftw_arrays(K * howmany.os, out, ro, io);
	printf("fftw forward split vs. interleaved\n");
	print_error(e);
	print_status(e <= max_error_fftw_split);

	plan1 = fftw_plan_dft_1d(CONV_WIN_SIZ, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
	plan2 = fftw_plan_guru_split_dft(1, &dim, 0, 0, ii, ri, io, ro, FFTW_ESTIMATE);
	
	fill_test_array(in, ri, ii);
	for(int k = 0; k < K; k++) fftw_execute_dft(plan1, in + k * howmany.is, out + k * howmany.os);
	fftw_execute(plan2);
	fftw_destroy_plan(plan1);
	fftw_destroy_plan(plan2);

	e = compare_fftw_arrays(K * howmany.os, out, ro, io);
	printf("fftw backward split vs. interleaved\n");
	print_error(e);
	print_status(e <= max_error_fftw_split);

	plan1 = fftw_plan_dft_r2c_1d(CONV_WIN_SIZ, ri, out, FFTW_ESTIMATE);
	plan2 = fftw_plan_guru_split_dft_r2c(1, &dim, 0, 0, ri, ro, io, FFTW_ESTIMATE);

	std::fill(ro, ro + 6 * K * CONV_WIN_SIZ, 0.0);
	fill_test_array(in, ri, ii);
	for(int k = 0; k < K; k++) fftw_execute_dft_r2c(plan1, ri + k * howmany.is, out + k * howmany.os);
	fftw_execute(plan2);
	fftw_destroy_plan(plan1);
	fftw_destroy_plan(plan2);

	e = 0;
	for(int k = 0; k < K; k++) {
		e += compare_fftw_arrays(CONV_WIN_SIZ/2 + 1, out + k*CONV_WIN_SIZ, ro + k*CONV_WIN_SIZ, io + k*CONV_WIN_SIZ);
	}
	printf("fftw real-to-complex split vs. interleaved\n");
	print_error(e);
	print_status(e <= max_error_fftw_split);

	conv_free(array);
	
}

#define CCONV_TEST_SAMPLE_FREQ 16000
#define CCONV_TEST_SIZE 8192

#define CCONV_TEST_WINDOW_SIZE 999
#define CCONV_TEST_WINDOW_FREQ { 50, 100, 200, 400, 1000, 2000, 4000 }
#define CCONV_TEST_WINDOW_STD 0.01
#define CCONV_TEST_REPEAT 1

#define CCONV_TEST_SIGNAL_SIZE CCONV_TEST_SIZE
#define	CCONV_TEST_SIGNAL_PITCH 100
#define CCONV_TEST_SIGNAL_HARM_POWER { 1, 0.7, 0.3 }
#define CCONV_TEST_SIGNAL_HARM_PHASE { 0, M_PI_4, M_PI/3 }

template<typename T>
void cconv_generate_test_window(T *window_re, T *window_im, freq_t freq) {
	for(int i = 0; i < CCONV_TEST_WINDOW_SIZE; i++) {
		T t = T(i - CCONV_TEST_WINDOW_SIZE/2) / CCONV_TEST_SAMPLE_FREQ;
		T power = gauss_win<T>(t, 0, CCONV_TEST_WINDOW_STD);
		T omega_t = 2 * M_PI * freq * t;
		window_re[i] = power * cos(omega_t);
		window_im[i] = power * sin(omega_t);
	}
}

template<typename T>
void cconv_generate_test_signal(T *signal) {
	size_t size = CCONV_TEST_SIGNAL_SIZE;
	T freq = CCONV_TEST_SAMPLE_FREQ;
	T pitch = CCONV_TEST_SIGNAL_PITCH;
	T power[] = CCONV_TEST_SIGNAL_HARM_POWER;
	T phase[] = CCONV_TEST_SIGNAL_HARM_PHASE;

	for(int i = 0; i < size; i++) {
		T t = T(i) / CCONV_TEST_SAMPLE_FREQ;
		T value = 0.0;
		for(int j = 0; j < sizeof(power) / sizeof(T); j++) {
			value += power[j] * sin(2 * M_PI * (1+j) * pitch * t + phase[j]);
		}
		signal[i] = value;
	}
}

/* void cconv_generate_testdata(
	double *signal, 
	double *window_re, 
	double *window_im, 
	double *cconv_re_naive, 
	double *cconv_im_naive,
	double window_freq) 
{
	typedef long double extended;

	extended signal_ext[CCONV_TEST_SIZE];
	extended window_ext[2*CCONV_TEST_SIZE];
	extended cconv_ext[2*CCONV_TEST_SIZE];
	extended 
		*window_re_ext = window_ext,
		*window_im_ext = window_ext + CCONV_TEST_SIZE,
		*cconv_re_ext = cconv_ext,
		*cconv_im_ext = cconv_ext + CCONV_TEST_SIZE;

	// generate inputs
	cconv_generate_test_signal(signal_ext);
	std::fill(window_ext, window_ext + 2*CCONV_TEST_SIZE, 0.0);
	cconv_generate_test_window(window_re_ext, window_im_ext, window_freq);

	// compute cconv with long double
	cconv_naive<extended, CCONV_TEST_SIZE>(signal_ext, CCONV_TEST_SIGNAL_SIZE, window_re_ext, CCONV_TEST_WINDOW_SIZE, cconv_re_ext);
	cconv_naive<extended, CCONV_TEST_SIZE>(signal_ext, CCONV_TEST_SIGNAL_SIZE, window_im_ext, CCONV_TEST_WINDOW_SIZE, cconv_im_ext);

	// downgrade long double to double
	std::copy(signal_ext, signal_ext + CCONV_TEST_SIZE, signal);
	std::copy(window_re_ext, window_re_ext + CCONV_TEST_SIZE, window_re);
	std::copy(window_im_ext, window_im_ext + CCONV_TEST_SIZE, window_im);
	std::copy(cconv_re_ext, cconv_re_ext + CCONV_TEST_SIZE, cconv_re_naive);
	std::copy(cconv_im_ext, cconv_im_ext + CCONV_TEST_SIZE, cconv_im_naive);
}
*/
template<typename T, typename T2>
void cconv_calc_error(const T *cconv_re, const T *cconv_im, const T2 *cconv_re_ideal, const T2 *cconv_im_ideal) {

	// calculate error
	T2 e_re = 0.0, 
		e_im = 0.0,
		max_re = 0.0, min_re = 0.0,
		max_im = 0.0, min_im = 0.0;
	for(int i = 0; i < CCONV_TEST_SIZE; i++) {
		e_re += fabs((cconv_re[i] - cconv_re_ideal[i]) / cconv_re_ideal[i]);
		e_im += fabs((cconv_im[i] - cconv_im_ideal[i]) / cconv_im_ideal[i]);
		min_re = std::min(min_re, fabs(cconv_re[i]));
		min_im = std::min(min_im, fabs(cconv_im[i]));
		max_re = std::max(max_re, fabs(cconv_re[i]));
		max_im = std::max(max_im, fabs(cconv_im[i]));
	}
	e_re = e_re / CCONV_TEST_SIZE;
	e_im = e_im / CCONV_TEST_SIZE;

	// print error
	print_error(e_re);
	print_status(e_re < max_error_cconv_fast);
	print_error(e_im);
	print_status(e_im < max_error_cconv_fast);
	print_error((e_re + e_im)/2);
}

/*
void test_cconv_naive() {

	typedef long double extended;
	typedef double cut_t;

	// generate test data
	extended signal[CCONV_TEST_SIGNAL_SIZE];
	extended window_re[CCONV_TEST_WINDOW_SIZE];
	extended window_im[CCONV_TEST_WINDOW_SIZE];
	cconv_generate_test_signal(signal);
	cconv_generate_test_window(window_re, window_im, 200);
	
	// compute cconv
	extended cconv_re[CCONV_TEST_SIZE];
	extended cconv_im[CCONV_TEST_SIZE];
	tic();
	for(int i = 0; i < CCONV_TEST_REPEAT; i++) {
		cconv_naive<extended, CCONV_TEST_SIZE>(signal, CCONV_TEST_SIGNAL_SIZE, window_re, CCONV_TEST_WINDOW_SIZE, cconv_re);
		cconv_naive<extended, CCONV_TEST_SIZE>(signal, CCONV_TEST_SIGNAL_SIZE, window_im, CCONV_TEST_WINDOW_SIZE, cconv_im);
	}
	print_time();

	// copy from long double to double
	cut_t cconv_re_cut[CCONV_TEST_SIZE];
	cut_t cconv_im_cut[CCONV_TEST_SIZE];
	std::copy(cconv_re, cconv_re + CCONV_TEST_SIZE, cconv_re_cut);
	std::copy(cconv_im, cconv_im + CCONV_TEST_SIZE, cconv_im_cut);

	// calculate error
	extended e = 0.0;
	for(int i = 0; i < CCONV_TEST_SIZE; i++) {
		extended re = cconv_re[i] - extended(cconv_re_cut[i]);
		extended im = cconv_im[i] - extended(cconv_im_cut[i]);
		e += re > 0 ? re : -re;
		e += im > 0 ? im : -im;
	}

	// print error
	print_error(e);
	print_status(e <= max_error_cconv_naive);

	// save computed cconvs
	io::array_to_file(signal, CCONV_TEST_SIGNAL_SIZE, "cconv-naive-signal.bin");
	io::array_to_file(window_re, CCONV_TEST_WINDOW_SIZE, "cconv-naive-window-re.bin");
	io::array_to_file(window_im, CCONV_TEST_WINDOW_SIZE, "cconv-naive-window-im.bin");
	io::array_to_file(cconv_re, CCONV_TEST_SIZE, "cconv-naive-cconv-re.bin");
	io::array_to_file(cconv_im, CCONV_TEST_SIZE, "cconv-naive-cconv-im.bin");
	io::array_to_file(cconv_re_cut, CCONV_TEST_SIZE, "cconv-naive-cconv-re-cut.bin");
	io::array_to_file(cconv_im_cut, CCONV_TEST_SIZE, "cconv-naive-cconv-im-cut.bin");

}


void test_cconv_tc4() {

	double signal[CCONV_TEST_SIZE];
	double window_re[CCONV_TEST_SIZE];
	double window_im[CCONV_TEST_SIZE];
	double cconv_re_naive[CCONV_TEST_SIZE];
	double cconv_im_naive[CCONV_TEST_SIZE];
	freq_t window_freq[] = CCONV_TEST_WINDOW_FREQ;
	int n_freq = sizeof(window_freq) / sizeof(freq_t);

	double cconv_re[CCONV_TEST_SIZE];
	double cconv_im[CCONV_TEST_SIZE];

	unsigned long long time = 0;
	for(int i = 0; i < CCONV_TEST_REPEAT; i++) {
		for(int j = 0; j < n_freq; j++) {

			// generate test data
			cconv_generate_testdata(signal, window_re, window_im, cconv_re_naive, cconv_im_naive, window_freq[j]);

			// compute cconv with toom-cook 4x4
			tic();
			base<double,CCONV_TEST_SIZE>::cconv(signal, window_re, cconv_re);
			base<double,CCONV_TEST_SIZE>::cconv(signal, window_im, cconv_im);
			time += toc();

			// calculate errors
			printf("window freq: %lg\n", double(window_freq[j]));
			cconv_calc_error(cconv_re, cconv_im, cconv_re_naive, cconv_im_naive);

		}
	}

	// print calculation time:
	printf("%lg ms per cconv\n", double(time)/CCONV_TEST_REPEAT/n_freq); 

	// save files
	io::array_to_file(cconv_re, CCONV_TEST_SIZE, "cconv-tc4-cconv-re.bin");
	io::array_to_file(cconv_im, CCONV_TEST_SIZE, "cconv-tc4-cconv-im.bin");
	
}

void test_cconv_fft() {

	// generate test data
	auto_ptr<double> memory(new double[11*CCONV_TEST_SIZE]);
	double 
		*signal = memory.get(),
		*window_re = memory.get() + CCONV_TEST_SIZE,
		*window_im = memory.get() + 2*CCONV_TEST_SIZE,
		*cconv_re_naive = memory.get() + 3*CCONV_TEST_SIZE,
		*cconv_im_naive = memory.get() + 4*CCONV_TEST_SIZE,
		*signal_fft = memory.get() + 5*CCONV_TEST_SIZE,
		*signal_fft_re = signal_fft,
		*signal_fft_im = signal_fft + CCONV_TEST_SIZE,
		*window_fft_re = memory.get() + 7*CCONV_TEST_SIZE,
		*window_fft_im = memory.get() + 8*CCONV_TEST_SIZE,
		*cconv_re = memory.get() + 9*CCONV_TEST_SIZE,
		*cconv_im = memory.get() + 10*CCONV_TEST_SIZE;
	freq_t window_freq[] = CCONV_TEST_WINDOW_FREQ;
	int n_freq = sizeof(window_freq) / sizeof(freq_t);

//	assert(CCONV_TEST_SIZE == CONV_WIN_SIZ);

	// real calculations start
	unsigned long long time = 0;
	for(int i = 0; i < CCONV_TEST_REPEAT; i++) {
		for(int j = 0; j < n_freq; j++) {

			// generate test data
			cconv_generate_testdata(signal, window_re, window_im, cconv_re_naive, cconv_im_naive, window_freq[j]);

			// compute cconv with toom-cook 4x4
			tic();
			cconv_calc_A(signal, signal_fft_re, signal_fft_im);
			cconv_calc_BC(window_re, window_im, window_fft_re, window_fft_im);
			cconv_normalize(signal_fft, 2 * CONV_WIN_SIZ);
			cconv(signal_fft_re, signal_fft_im, window_fft_re, window_fft_im, cconv_re, cconv_im);
			time += toc();

			// calculate errors
			printf("window freq: %lg\n", double(window_freq[j]));
			cconv_calc_error(cconv_re, cconv_im, cconv_re_naive, cconv_im_naive);

		}
	}

	// print calculation time:
	printf("%lg ms per cconv\n", double(time)/CCONV_TEST_REPEAT/n_freq); 

	// save files
	io::array_to_file(cconv_re, CCONV_TEST_SIZE, "cconv-fft-cconv-re.bin");
	io::array_to_file(cconv_im, CCONV_TEST_SIZE, "cconv-fft-cconv-im.bin");
	
}
*/
void test_fft_approx() {

	double array[CONV_WIN_SIZ * 11 + 2];
	double *a = SPL_MEMORY_ALIGN(array);
	double *b = a + CONV_WIN_SIZ;
	double *c = b + CONV_WIN_SIZ;
	double *Ar= c + CONV_WIN_SIZ;
	double *Ai= Ar+ CONV_WIN_SIZ;
	double *B = Ai+ CONV_WIN_SIZ;
	double *C = B + CONV_WIN_SIZ;
	double *ab1 = C + CONV_WIN_SIZ;
	double *ac1 = ab1 + CONV_WIN_SIZ;
	double *ab2 = ac1 + CONV_WIN_SIZ;
	double *ac2 = ab2 + CONV_WIN_SIZ;

	std::fill(ab1, ab1 + 4*CONV_WIN_SIZ, 0.0);

	int Ws = 12;
	freq_t Fr = 200, dF = 3; // какие тут частоты указывать - не имеет особого значения
	double std = mask_std(Fr, 1.0);
	for(int m = -Ws; m <= Ws; m++) {
		a[Ws + m] = gauss_win(Fr + m * dF, Fr, std);
	}
	std::fill(a + 2*Ws + 1, a + CONV_WIN_SIZ, 0.0);

	fill_test_array((fftw_complex*)B, b, c);

	cconv_calc_A(a, Ar, Ai);
	cconv_normalize(Ar, 2 * CONV_WIN_SIZ);
	cconv_calc_BC(b, c, B, C);
	cconv(Ar, Ai, B, C, ab1, ac1);

	for(int i = Ws; i < CONV_WIN_SIZ - Ws; i++) {
		double sumb = 0.0, sumc = 0.0;
		for(int m = -Ws; m <= Ws; m++) {
			sumb += a[Ws + m] * b[i - m];
			sumc += a[Ws + m] * c[i - m];
		}
		ab2[i] = sumb;
		ac2[i] = sumc;
	}

	io::array_to_file(ab1, CONV_WIN_SIZ, "cconv-1-fft.bin");
	io::array_to_file(ac1, CONV_WIN_SIZ, "cconv-2-fft.bin");
	io::array_to_file(ab2, CONV_WIN_SIZ, "cconv-1-naive.bin");
	io::array_to_file(ac2, CONV_WIN_SIZ, "cconv-2-naive.bin");

	double e1 = 0.0, e2 = 0.0;
	for(int i = Ws; i < CONV_WIN_SIZ - Ws; i++) {
		e1 += fabs(ab1[i+Ws] - ab2[i]);
		e2 += fabs(ac1[i+Ws] - ac2[i]);
	}

	print_error(e1);
	print_error(e2);
	print_status(e1 < max_error_fft_approx);

}

void test_pitch_track() {
	
	
}

void find_max_harm(double* mharm, int ws_s, short feat_num, int harm_num, int ptr, bool norm, ofstream<float>& output2,const scale_class& sc,const spl::spectrum& sp,int k1, int k2, const  mask& msk)
{
	int d =harm_num;
	int dd = harm_num+harm_num;
	if (harm_num>20) harm_num +=5; else harm_num+=10;


	double freq=0;
    spectrum_t* tmp = new spectrum_t[harm_num];
	double* avg = new double[harm_num];
	double* last_freq = new double[harm_num];
	double* intensity = new double[harm_num];
	bool* used = new bool[harm_num];
	int* max = new int[harm_num]; // массив номеров каналов с макс по интенсивности гармониками
	int* count = new int[harm_num]; // массив количества использованных семпплов в формантах.

	bool inside_segment = false;
	spectrum_t max_segment = -std::numeric_limits<double>::max();
	int num_segment = 0;
	int ptr2 = ptr;
	for (int k=0; k<harm_num; k++)
	{
		max[k] = 0;
		
		tmp[k] = -std::numeric_limits<double>::max();
		avg[k] = 0;
		count[k] = 0;
		intensity[k] =0;
		last_freq[k] =0;
	}
	if (ws_s == 1) //не усредняем
	{
		for (int j=k1; j<k2; j++)
		{
			
			if (msk.Z[ptr*sp.K+j] != 0)
			{
				inside_segment = true;
				if (max_segment < sp.Y[ptr*sp.K+j])
				{
					max_segment = sp.Y[ptr*sp.K+j];
					num_segment = j;
				}
			}
			else 
			{ 
				if (inside_segment) //закончился наш сегмент не маскируемых частот
				{
					inside_segment = false; /* обнуляем максимумы в сегменте */
					for (int k=0; k<harm_num; k++)
					{				
						if (tmp[k]<max_segment)
						{
							for (int l=harm_num-1; l>k;l--)
							{
								max[l]=max[l-1];
								tmp[l]=tmp[l-1];
							}
							tmp[k] = max_segment;
							max[k] = num_segment;
							break;
						}
					}
					max_segment = -std::numeric_limits<double>::max();
					num_segment = 0;
				}
			}
		}
		delete [] tmp;
		delete [] count;
		harm_num = d; //меняем количество гармоник на изначальное
		for (int k=0; k<harm_num; k++)
		{			
			mharm[k] = sc.Fr[max[k]];
		}
		delete [] max;
		delete [] intensity;
	}
	else //вычисляем значения в окне
		// алгоритм следующий: в первом семпле находим максимумы в не маскируемых сегментах и запоминаем.
		//в последующих сегментах находим X максимумов и складываем их частоты с ближайшими частотами (из предыдщего сепмла), запоминаем.
		//при этом суммируем значения интенсивностей, чтобы затем упорядочить все по убыванию.
	{
		for (int n=0; n<ws_s; n++)	
		{
			for (int k=0; k<harm_num; k++) //обнуляем массивы
			{
				max[k] = 0;
				tmp[k] = -std::numeric_limits<double>::max();
				used[k] = false;
			}

			//1. находим массив максимальных гармоник в данном семпле.

			for (int j=k1; j<k2; j++)
		{
			
			if (msk.Z[ptr2*sp.K+j] != 0)
			{
				inside_segment = true;
				if (max_segment < sp.Y[ptr2*sp.K+j])
				{
					max_segment = sp.Y[ptr2*sp.K+j];
					num_segment = j;
				}
			}
			else 
			{ 
				if (inside_segment) //закончился наш сегмент не маскируемых частот
				{
					
					for (int k=0; k<harm_num; k++)
					{				
						if (tmp[k]<max_segment)
						{
							for (int l=harm_num-1; l>k;l--)
							{
								max[l]=max[l-1];
								tmp[l]=tmp[l-1];
							}
							tmp[k] = max_segment;
							max[k] = num_segment;
							break;
						}
					}
					inside_segment = false; /* обнуляем максимумы в сегменте */
					max_segment = -std::numeric_limits<double>::max();
					num_segment = 0;
				}
			}
		}
			//2. добавляем к наиболее близким.
			for (int k=0; k<harm_num; k++)
			{
				freq = sc.Fr[max[k]];
				//avg[k]+= sc.Fr[max[k]];
				if (n!=0) //не первый семпл, добавляем к наиболее близким
				{
					if (max[k] != 0) //если не нашли максимума, то ничего не делаем
					{
						int close = 0;
					
						while (used[close]) close++; //ищем первую неиспользованную гармонику
					
						for (int l=close+1; l<harm_num;l++) //ищем наиболее близкую частоту, при этом еще не использованную
						{
							if ((!used[l])&&abs(last_freq[l]-freq)<abs(last_freq[close]-freq)) close = l; 
						}
						if ((last_freq[close]!=0)&&(abs(last_freq[close]-freq) > 50)) //если большая разница с ближайшим, то пытаемся найти свободное место
						{
							int l;
							for (l=0; l<harm_num;l++)
							{
								if (last_freq[l] == 0) break;
							}
							if (l!=harm_num) {close = l; }

						}
						avg[close]+= freq;
						count[close]++;
						//if ( -std::numeric_limits<double>::infinity() != tmp[k]) 
						intensity[close] += tmp[k];
						used[close] = true;
						last_freq[close] = freq;
					}
				}
				else //первый семпл, просто запоминаем
				{
					if (max[k] == 0) { intensity[k] = 0; last_freq[k] = 0; avg[k] = 0;} else {
					avg[k]+= freq; count[k]++;
					if ( -std::numeric_limits<double>::infinity() == tmp[k]) 
						intensity[k] = 0; else
					intensity[k] = tmp[k];
					last_freq[k] = freq;}
				}
				
			}
			ptr2++;
			inside_segment = false; /* обнуляем максимумы в сегменте */
			max_segment = -std::numeric_limits<double>::max();
			num_segment = 0;
		}
		//а затем находим средние
			for (int k=0; k<harm_num; k++)
			{
				if (count[k] != 0 ) avg[k] = avg[k]/count[k];
				else avg[k] = 0;
			}
			double swap;
			double srch_max = -std::numeric_limits<double>::max();
			int srch_ind = 0; int srch_len = 0;
			//сортируем по интенсивностям
			for (int m=0;m<harm_num; m++) {
				srch_max = -std::numeric_limits<double>::max(); srch_ind = m;
			
			/*for (int k=m; k<harm_num; k++)//ищем максимальную длину
			{
				if (count[k] > srch_len) {srch_len = count[k];}
			}
			for (int k=m; k<harm_num; k++)
			{
				if ((count[k]==srch_len)&&(intensity[k]>srch_max)) {srch_max=intensity[k]; srch_ind = k;}
			}*/
			for (int k=m; k<harm_num; k++)
			{
				if (intensity[k]>srch_max) {srch_max=intensity[k]; srch_ind = k;}
			}
			swap = avg[m]; 
			avg[m] = avg[srch_ind];
			avg[srch_ind] = swap;
			swap = intensity[m]; 
			intensity[m] = intensity[srch_ind];
			intensity[srch_ind] = swap;
			}
			harm_num = d; //меняем количество гармоник на изначальное
			for (int k=0; k<harm_num; k++)
			{
				mharm[k] = avg[k];
			}
			delete [] tmp;
			delete [] count;
			delete [] max;
			delete [] intensity;
	}
	//а теперь выводим результат

		if (norm) 
	{  
		for (int k=0; k<harm_num; k++)
			{
				mharm[k]= ((mharm[k]-sc.Fr[k1])/(sc.Fr[k2]-sc.Fr[k1]))-0.5;
			}
	}
	
	//выводим
  	for (int k=0; k<harm_num; k++)
	{
		output2.put((float)mharm[k]);
		if (feat_num > 1) output2.put((float)(mharm[k]-mharm[k+d]));
		if (feat_num > 2) output2.put((float)(mharm[k]-mharm[k+dd]));
	}


	//запоминаем дельты
	for (int j = 0; j<harm_num; j++)
	{
		mharm[j+dd] = mharm[j+d];
		mharm[j+d] = mharm[j];
	}

	delete [] avg;
	//delete [] tmp;
	
    delete [] last_freq;
	
	delete [] used;
	
}


void out_channels(double* avg,short feat_num,int ws_s, int ptr, bool norm, ofstream<float>& output2,const scale_class& sc,const spl::spectrum& sp,int k1, int k2, const  mask& msk)
{
	int d = msk.K;
	int dd = msk.K+msk.K;
	int totalK = feat_num*msk.K;
	
	int ptr2 = ptr*msk.K;

	int totalN = ptr+ws_s;
	if (totalN>msk.N) { totalN = msk.N-ptr;}
	for (int j = ptr; j<totalN; j++)
	{
		for (int i =0; i<msk.K; i++)
		{
			avg[i] += sp.Y[ptr2]*msk.Z[ptr2];
			ptr2++;
		}
	}
	for (int i=0; i<msk.K; i++)
	{
		avg[i]=avg[i]/totalN;
	}

	//возможно тут нормализация?

	//выводим результат

	for (int i=0; i<msk.K; i++)
	{
		output2.put(avg[i]);
	}
	if (totalK>msk.K) // если есть дельты
	{
		for (int i=msk.K; i<totalK; i++)
	{
		output2.put(avg[i-msk.K]-avg[i]);
	}
	
	}

	//запоминаем дельты
	for (int j = 0; j<msk.K; j++)
	{
		avg[j+dd] = avg[j+d];
		avg[j+d] = avg[j];
	}
	
}
void windowing(int ws, int wm, short feat_num, int harm_num, bool norm,const spl::spectrum& sp,const  mask& msk,short* channels,const char *output_path, short low_ch, short ch_num,const  scale_class& sc, bool outChannels)
{
	// выходной файл
	ofstream<float> output2(output_path);
	
	int ws_s,wm_s;
	//вывод как нам надо.
	if (ws != 0) {
	 ws_s = ((ws*sp.F)/1000); //окно в семплах
	 wm_s = ((wm*sp.F)/1000); //сдвиг в семплах
	}
	else { ws_s=1; wm_s=1;}
	
		
	int result_n = (msk.N-(ws_s-wm_s))/wm_s; //всего результирующих значений каналов ЧОТ

	//выводим сколько у нас признаков
	int feats =feat_num+(harm_num*feat_num);
	if (outChannels) feats += msk.K*feat_num;
	output2.put(feats);
	//выводим количество блоков на один признак
	if (result_n*wm_s<msk.N) output2.put(result_n+1);
	else output2.put(result_n);
	
	//выводим количество окон в 1 секунде
	output2.put(sp.F/(wm_s));

	int ptr = 0;
	int total_v =0, total_nv =0;
	float current=0, delta = -1, ddelta = -1;
	double* mharm = new double[harm_num*3];
	for (int j=0; j<harm_num*3; j++)
	{
		if (norm) mharm[j] = -0.5; 
		else mharm[j] = 0;
	}
	double* K = new double[msk.K*3];
	for (int j=0; j<msk.K*3; j++)
	{
		K[j]=0;
	}
	int k1 = scale_find_freq(sc, original.Fov);
	int k2 = scale_find_freq(sc,original.Fv);
		
	for (int i=0; i<result_n; i++)
	{
		
	if (harm_num!=0) if (((ws_s==1)&&(channels[ptr] !=-1))||(ws_s!=1))   find_max_harm(mharm,ws_s, feat_num, harm_num,ptr,norm, output2,sc, sp,k1,k2,msk);
		if (ws_s == 1) { if (channels[ptr] != -1) current = sc.Fr[channels[ptr]]; else current = -1; }
		else 
		{
			short not_voiced =0;
			double avg_channel =0;
			
			for (int j=0; j<ws_s; j++)	
			{
				if (channels[ptr+j] == -1) not_voiced++;
				else avg_channel += sc.Fr[channels[ptr+j]];
			}
		
			//если хотя бы половина в окне невокализованных, считаем что целиком окно невокализованное
			if ((not_voiced>=ws_s/2)) { current=-1; total_nv ++;}
			else current = (avg_channel/(ws_s-not_voiced));
		}


		if (norm&&current!=-1) 
		{
			current = (current/200)-0.5;
		}

		//очень странное условие. Вообще не понимаю, почему так?
		//if ((ws_s==1)&&(current !=-1))
		//вообще вопрос, надо ли выводить по разному, когда ЧОТ =0 и когда != 0.Корректно ли выводить 200-0 Гц дельты?
		
		
		//вывод
		//if (current !=-1)
		{
			
			output2.put(current);
			if (feat_num>1) //выводим 1 дельту
			{output2.put(current-delta);}
			if (feat_num>2) //выводим 2 дельту
			{output2.put(current-ddelta);}
		}
		
		ddelta = delta;
		delta = current;

		if (outChannels) out_channels(K,feat_num, ws_s, ptr,norm, output2,sc, sp,k1,k2,msk);
		ptr+=wm_s;
		
	}
	if (ptr<msk.N)	
		{
			if (harm_num!=0) find_max_harm(mharm,msk.N-ptr, feat_num, harm_num,ptr,norm, output2,sc, sp,k1,k2,msk);
			current = -1;
			
			output2.put(current);
			if (feat_num>1) //выводим 1 дельту
			{output2.put(current-delta);}
			if (feat_num>2) //выводим 2 дельту
			{output2.put(current-ddelta);}
			if (outChannels) out_channels(K,feat_num, ws_s, ptr,norm, output2,sc, sp,k1,k2,msk);
		}

	
	delete [] mharm;
	delete [] K;
	output2.close();
}

void get_mtf(int ws, int wm, short feat_num, int harm_num, bool norm,const char* wav_path, const char* ch_path,const char* voc_path,bool chan_out)
{
	{
	// открываем шкалу частот
	scale_class sc(scale_std);

	//делаем фльтрацию из вав в память
	//size_t written = filter_wave_file(sc, signal_wav_std, spectrum_test, original.ksi);
	spl::spectrum sp;
	size_t written = spl::filter_wmemory(sc,wav_path,sp,original.ksi, original.K);

	//iwstream in_str(signal_wav_std);
	// открываем маску
	//ifstream<mask_t> input(mask_byte_std);

	//сделать так, чтобы открыть поток из фильтра.

	// mask parameters
	mask_parameters p;
	p.border_effect = original.border_effect != 0.0;
	p.ksi = original.ksi;
	p.rho = original.rho;
	p.delta = original.delta;
	p.allow_fast = true;
	mask msk;
	msk.K = sp.K;
	msk.N = sp.N;
	msk.Z = new mask_t[msk.K*msk.N];

	written = mask_memory(sc, sp, msk, p);
	

	double std = filter_std(sc.Fr[0], original.F);
	int Ws = int(ceil(gauss_border(std, original.ksi)));
	
	
	//входной поток
	
	imstream<mask_t> input(msk.Z, msk.K*msk.N);

	//выходной поток в память
	
	short* channels = new short[msk.N];

	omstream<short> output(channels,msk.N);

	// параметры генерации шаблонов
	pitch_parameters pp;
	pp.Nh = original.Ng;
	pp.F1 = original.Fon;
	pp.F2 = original.Fov;

	// собственно выделение ЧОТ

	written = pitch_track(sc, input, output, p, pp);

	// 1. вычисление k1, k2
	int k1 = scale_find_freq(sc, pp.F1);
	int k2 = scale_find_freq(sc, pp.F2);
	int nk = k2 - k1 + 1;

	
	//делим на окна и выводим
	if (voc_path == nullptr) windowing(ws,wm,feat_num , harm_num, norm, sp, msk, channels, ch_path, k1, nk, sc,chan_out);
	
		
	//printf("%i",total_v);
	//total_nv = result_n-total_v;

	if (voc_path != nullptr) {
	
	//проведем сегментацию

	imstream<short> input_ch(channels, msk.N);
	// выходной файл
	ofstream<short> output_seg(vocal_chan_test);
	short* segmentation = new short[msk.N];
	omstream<short> mem_seg(segmentation, msk.N);
	
	// параметры сегментации
	vocal_parameters voc_p;
	voc_p.minV = original.minV;
	voc_p.minNV = original.minNV;

	// частота дискретизации
	freq_t F = original.F;

	// собственно сегментация
	written = vocal_segment(input_ch, mem_seg, F, voc_p);


	//делим на окна и выводим.
	windowing(ws,wm,feat_num ,harm_num, norm, sp, msk, segmentation, voc_path,k1,nk,sc,chan_out);

	//очищаем память
	delete [] segmentation;}	
	delete [] msk.Z;
	delete [] channels;

	
}
}



void test_vocal_segment() {
	
{
	// открываем шкалу частот
	scale_class sc(scale_std);

	// открываем каналы ЧОТ
	ifstream<short> input(pitch_chan_test);

	// выходной файл
	ofstream<short> output(vocal_chan_test);

	// параметры сегментации
	vocal_parameters p;
	p.minV = original.minV;
	p.minNV = original.minNV;

	// частота дискретизации
	freq_t F = original.F;

	// собственно сегментация
	size_t written = vocal_segment(input, output, F, p);
}

	// сравниваем выход и эталон
	double e = compare_streams<short>(vocal_chan_std, vocal_chan_test);

	print_error(e);
	print_status(e <= max_error_vocal_segment);
}

void test_mask_win_size() {
	int Ks[] = { 128, 147, 168, 194, 222, 256, 294, 337, 388, 445, 512 };
	double Es[] = { 0.1, 0.063096, 0.039811, 0.025119, 0.015849, 0.01, 0.0063096, 0.0039811, 0.0025119, 0.0015849, 0.001 };
	scale_parameters sp;

	mask_parameters mp = DEFAULT_MASK_PARAMETERS;
	mp.allow_fast = false;
	mp.border_effect = true;
	char scale_file_name[30];
	for(int i = 0; i < sizeof(Ks)/sizeof(int); i++) {
		sp.a = sp.i = 1;
		sp.b = sp.j = Ks[i];
		sp.Fi = 50.0;
		sp.Fj = 4000.0;
		scale_class sc(sp, scale_model);
		for(int j = 0; j < sizeof(Es)/sizeof(double); j++) {
			mp.ksi = Es[j];
			mask_filters_class mf(sc, mp);
			printf("K = %d, e = %lf, L = %d\n", Ks[i], Es[j], mf.Ws);
			sprintf(scale_file_name, "scale_K%dE%lf.bin", Ks[i], Es[j]);
			mf.save(scale_file_name);
		}
	}
}

void help();
void test_all();

struct test {
	const char *name;
	typedef void (*test_ptr)();
	test_ptr func;
} tests[] = {
	{ "all", test_all },
	{ "iwstream", test_iwstream },
	{ "scale_generate", test_scale_generate },
	{ "scale_form_estimate", test_scale_form_estimate },
	{ "scale_interp", test_scale_interp },
	{ "filter_binary", test_filter_binary_file },
	{ "filter_wave", test_filter_wave_file },
	{ "spectrum_filters", test_spectrum_filters },
	{ "mask_bytes_naive", test_mask_bytes_naive },
	{ "mask_bytes_fast", test_mask_bytes_fast },
	{ "mask_bits_naive", test_mask_bits_naive },
	{ "mask_bits_fast", test_mask_bits_fast },
	{ "mask_extend_reduce", test_mask_extend_reduce },
	{ "mask_win_size", test_mask_win_size },
	{ "iobit", test_iobit },
	{ "iobuf", test_iobuf },
	{ "imicstream", test_imicstream },
	{ "fftw_split", test_fftw_split },
	{ "fft_approx", test_fft_approx },
/*	{ "cconv_naive", test_cconv_naive },
	{ "cconv_tc4", test_cconv_tc4 },
	{ "cconv_fft", test_cconv_fft },*/
	{ "pitch_track", test_pitch_track },
	{ "vocal_segment", test_vocal_segment },
	{ "help", help },
};

struct helper{ const char *name; 
} help_[] = {{"Speech parameters estimation. Command line help:"},
			{"test.exe .wav_path output_path -ws* -wm* -nw** -X* -d* -dd* -n* -seg* "},
			{"* - not required"},			
			{"-ws int / window size, default 20 ms"},
			{"-wm int / window move, default 10 ms"},
			{"-nw int / do not use windows, samples only"},
			{"-X int / Output most intensive harmonics, not Fundamental frequency"},
			{"-K / Output channels"},
			{"-d / include deltas"},
			{"-dd / include deltas and double deltas"},
			{"-n / normalize"},
			{"-seg  / do speech algo voiced-unvoiced segmentation"},
			{""},


};

void help() {
	//printf("availible tests: \n");
	for(int i = 0; i < sizeof(help_)/sizeof(helper); i++) {
		printf("%s ", help_[i].name);
	}
	printf("\n");
}

void test_all() {
	for(int j = 0; j < sizeof(tests)/sizeof(test); j++) {
		if(tests[j].func == test_all) continue;
		printf("--- Attemting %s test ---\n", tests[j].name);
		tests[j].func(); 
	}
}

} // namespace {


// дебаг строки signal.wav signal.bin -ws 20 -wm 10 -X 10 -dd -n -nw -seg signal.bin
//"E:\temp\123\1 (2).wav" "E:\temp\123\mft2dN\1 (2).mft" -dd -n
int main(int argc, const char *argv[]) {
	
	read_original(original_std, &original);
	if(!check_original(&original)) {
		printf("Can't proceed with tests, as test original is different\n");
		return -1;
	}
	int ws = 20, wm=10; 
	short feat_num = 1;
	bool norm = false;
	const char * seg_out = nullptr;
	int harm_num = 0;
	bool channels = false; // для вывода К каналов

	if (argc <3) { help(); return 0;}

	for(int i = 3; i < argc; i++) {

		
			if(0 == strcmp(argv[i],"-ws")) { ws=atoi(argv[i+1]); i++; continue;}
			if(0 == strcmp(argv[i],"-wm")) { wm=atoi(argv[i+1]); i++; continue;}
			if(0 == strcmp(argv[i],"-n")) {norm = true; continue;}
			if(0 == strcmp(argv[i],"-seg")) {seg_out = argv[2]; continue;}
			if(0 == strcmp(argv[i],"-d")) {feat_num = 2; continue;}
			if(0 == strcmp(argv[i],"-dd")) {feat_num = 3; continue;}
			if(0 == strcmp(argv[i],"-nw")) {ws = 0; continue;}
			if(0 == strcmp(argv[i],"-X")) {harm_num=atoi(argv[i+1]); i++; continue;}
			if(0 == strcmp(argv[i],"-K")) {channels=true; continue;}




			/*	try {
					tests[j].func(); 
				} catch(const char *msg) {
					printf("=== Got error message: '%s' ===\n", msg);
				}
		
				break;
			}
			*/
		
		
	}
    get_mtf(ws,wm,feat_num, harm_num,norm,argv[1], argv[2],seg_out,channels);
	return 0;
}


