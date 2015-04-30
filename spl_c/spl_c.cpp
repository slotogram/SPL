/// 
/// \file  spl_c.cpp
/// \brief Реализация Си-интерфейса библиотеки SPL
/// 

#include "spl_types.h"
#include "scale.h"
#include "spectrum.h"
#include "mask.h"
#include "mask_helper.h"
#include "spl_c.h"
#include "iomem.h"
#include "iofile.h"
#include "iowave.h"
#include "iobit.h"
#include "vocal.h"

using namespace spl;

// disable: conversion from size_t to unsigned - possible loss of data
#pragma warning (disable: 4267)

extern "C" {

define_wrap_func(unsigned, spl_filter_memory_linear, 8, (K, N, F, F1, F2, s, sp, ksi));
define_wrap_func(unsigned, spl_filter_memory, 7, (K, N, F, sc, s, sp, ksi));
define_wrap_func(unsigned, spl_filter_binary_file, 6, (K, F, sc, in_file, out_file, ksi));
define_wrap_func(unsigned, spl_filter_wave_file, 5, (K, sc, in_file, out_file, ksi));

define_wrap_func(unsigned, spl_load_signal, 4, (file, N, buf, F));

define_wrap_func(unsigned, spl_mask_function_generate, 4, (K, Ws, sc, H));

define_wrap_func(unsigned, spl_mask_memory, 6, (K, N, sc, sp, m, ksi));
define_wrap_func(unsigned, spl_mask_binary_file, 6, (K, sc, in_file, out_file, out_bits, ksi));
define_wrap_func(unsigned, spl_mask_memory_ext, 10, (K, N, sc, sp, m, ksi, delta, rho, border_effect, allow_fast));
define_wrap_func(unsigned, spl_mask_binary_file_ext, 10, (K, sc, in_file, out_file, out_bits, ksi, delta, rho, border_effect, allow_fast));

define_wrap_func(unsigned, spl_pitch_memory, 6, (K, N, sc, m, p, ksi));
define_wrap_func(unsigned, spl_pitch_binary_file, 7, (K, N, sc, in_file, out_file, in_bits, ksi));
define_wrap_func(unsigned, spl_pitch_memory_ext, 12, (K, N, sc, m, p, ksi, delta, rho, border_effect, F1, F2, Nh));
define_wrap_func(unsigned, spl_pitch_binary_file_ext, 13, (K, N, sc, in_file, out_file, in_bits, ksi, delta, rho, border_effect, F1, F2, Nh));

/// Загрузка wav-файла.

unsigned spl_load_signal(char *file, int N, signal_t *buf, freq_t *F) {
	try {
		iwstream in(file);
		*F = in.freq();
		return in.read(buf, N);
	} catch(const char *) {
		return 0;
	} catch(...) {
		return 0;
	}
}

/// Генерация маскирующей функции 

unsigned spl_mask_function_generate(int K, int Ws, freq_t *sc, double *H) {
	mask_filters mf;
	mf.K = K;
	mf.Ws = Ws;
	mf.H = H;
	mask_parameters mp = DEFAULT_MASK_PARAMETERS;
	mp.allow_fast = false;
	mp.border_effect = false;
	freq_scale scale;
	scale.K = K;
	scale.Fr = sc;

	mask_filters_generate(scale, mf, mp);
	return K*Ws;
}

/// Генерация шкалы заданной формы по заданным параметрам.

#define define_scale_generate(name) \
void spl_scale_generate_##name(int K, freq_t F1, freq_t F2, freq_t *sc) { \
	scale_parameters p; p.a = p.i = 0; p.b = p.j = K-1; p.Fi = F1; p.Fj = F2; \
	freq_scale s; s.K = K; s.Fr = sc; \
	scale_generate_##name(s,p); } \
define_wrap_func(void, spl_scale_generate_##name, 4, (K, F1, F2, sc));

define_scale_generate(linear);
define_scale_generate(log);
define_scale_generate(mel);
define_scale_generate(bark);
define_scale_generate(model);

/// Фильтрация по заданным характеристикам шкалы.

unsigned spl_filter_memory_linear(
	int K,             ///< Количество каналов фильтрации
	int N,             ///< Количество отсчетов в сигнале
	freq_t F,          ///< Частота дискретизации сигнала
	freq_t F1,         ///< Нижняя частота фильтрации
	freq_t F2,         ///< Верхняя частота фильтрации
	const signal_t *s, ///< Массив сигнала (N элементов)
	spectrum_t *sp,    ///< Выходной массив спектра (K x N элементов)
	double ksi         ///< Точность вычислений
) {
	scale_parameters p; p.a = p.i = 0; p.b = p.j = K-1; p.Fi = F1; p.Fj = F2;
	signal S; S.N = N; S.F = F; S.X = const_cast<signal_t *>(s);
	spectrum Sp; Sp.K = K; Sp.N = N; Sp.Y = sp;
	try {
		scale_class Sc(p, scale_linear);
		return filter_memory(Sc, S, Sp, ksi);
	} catch(...) {
		return 0;
	}
}

/// Фильтрация в оперативной памяти.

unsigned spl_filter_memory(
	int K,             ///< Количество каналов фильтрации
	int N,             ///< Количество отсчетов в сигнале
	freq_t F,          ///< Частота дискретизации сигнала
	const freq_t *sc,  ///< Массив шкалы (K элементов)
	const signal_t *s, ///< Массив сигнала (N элементов)
	spectrum_t *sp,    ///< Выходной массив спектра (K x N элементов)
	double ksi         ///< Точность вычислений
) {
	freq_scale Sc; Sc.K = K; Sc.Fr = const_cast<freq_t *>(sc);
	signal S; S.N = N; S.F = F; S.X = const_cast<signal_t *>(s);
	spectrum Sp; Sp.K = K; Sp.N = N; Sp.Y = sp;
	return filter_memory(Sc, S, Sp, ksi);
}

/// Фильтрация файла.

unsigned spl_filter_binary_file(
	int K,                ///< Количество каналов фильтрации
	freq_t F,             ///< Частота дискретизации сигнала
	const freq_t *sc,     ///< Массив шкалы
	const char *in_file,  ///< Входной сигнал (имя файла)
	const char *out_file, ///< Выходной спектр (имя файла)
	double ksi            ///< Точность вычислений
) {
	freq_scale Sc; Sc.K = K; Sc.Fr = const_cast<freq_t *>(sc);
	return filter_binary_file(Sc, in_file, out_file, F, ksi);
}

/// Фильтрация wav-файла.

unsigned spl_filter_wave_file(
	int K,                ///< Количество каналов фильтрации
	const freq_t *sc,     ///< Массив шкалы
	const char *in_file,  ///< Входной сигнал (имя файла)
	const char *out_file, ///< Выходной спектр (имя файла)
	double ksi            ///< Точность вычислений
) {
	freq_scale Sc; Sc.K = K; Sc.Fr = const_cast<freq_t *>(sc);
	return filter_wave_file(Sc, in_file, out_file, ksi);
}

/// Маскировка в оперативной памяти.

unsigned spl_mask_memory(int K, int N, 
	const freq_t *sc, const spectrum_t *sp, mask_t *m, 
	double ksi)
{
	freq_scale Sc; Sc.K = K; Sc.Fr = const_cast<freq_t *>(sc);
	spectrum Sp; Sp.K = K; Sp.N = N; Sp.Y = const_cast<spectrum_t *>(sp);
	mask M; M.K = K; M.N = N; M.Z = m;
	return mask_memory(Sc, Sp, M, ksi);
}

unsigned spl_mask_memory_ext(int K, int N, 
	const freq_t *sc, const spectrum_t *sp, mask_t *m, 
	double ksi, double delta, double rho, bool border_effect, bool allow_fast)
{
	mask_parameters p; p.ksi = ksi; p.delta = delta; p.rho = rho; p.border_effect = border_effect; p.allow_fast = allow_fast;
	freq_scale Sc; Sc.K = K; Sc.Fr = const_cast<freq_t *>(sc);
	spectrum Sp; Sp.K = K; Sp.N = N; Sp.Y = const_cast<spectrum_t *>(sp);
	mask M; M.K = K; M.N = N; M.Z = m;
	return mask_memory(Sc, Sp, M, p);
}

/// Маскировка бинарного файла.

unsigned spl_mask_binary_file(int K, 
	const freq_t *sc, const char *in_file, const char *out_file, 
	bool out_bits, double ksi)
{
	freq_scale Sc; Sc.K = K; Sc.Fr = const_cast<freq_t *>(sc);
	return mask_binary_file(Sc, in_file, out_file, out_bits, ksi);
}

unsigned spl_mask_binary_file_ext(int K, 
	const freq_t *sc, const char *in_file, const char *out_file, 
	bool out_bits, double ksi, double delta, double rho, bool border_effect, bool allow_fast)
{
	mask_parameters p; p.ksi = ksi; p.delta = delta; p.rho = rho; p.border_effect = border_effect; p.allow_fast = allow_fast;
	freq_scale Sc; Sc.K = K; Sc.Fr = const_cast<freq_t *>(sc);
	return mask_binary_file(Sc, in_file, out_file, out_bits, p);
}

/// Отслеживание ЧОТ в памяти

unsigned spl_pitch_memory(int K, int N, 
	const freq_t *sc, const mask_t *m, short *p, 
	double ksi)
{
	freq_scale Sc; Sc.K = K; Sc.Fr = const_cast<freq_t *>(sc);
	mask_parameters pm = DEFAULT_MASK_PARAMETERS; pm.ksi = ksi;
	pitch_parameters pp = DEFAULT_PITCH_PARAMETERS;
	
	// инициализируем входной поток из памяти
	io::imstream<mask_t> in_str(m, N * K);

	// инициализируем выходной поток в память
	io::omstream<short> out_str(p, N);

	return pitch_track(Sc, in_str, out_str, pm, pp);
}

unsigned spl_pitch_memory_ext(int K, int N, 
	const freq_t *sc, const mask_t *m, short *p, 
	double ksi, double delta, double rho, bool border_effect, double F1, double F2, int Nh)
{
	freq_scale Sc; Sc.K = K; Sc.Fr = const_cast<freq_t *>(sc);
	mask_parameters pm; pm.ksi = ksi; pm.delta = delta; pm.rho = rho; pm.border_effect = border_effect;
	pitch_parameters pp; pp.Nh = Nh; pp.F1 = F1; pp.F2 = F2;

	// инициализируем входной поток из памяти
	io::imstream<mask_t> in_str(m, N * K);

	// инициализируем выходной поток в память
	io::omstream<short> out_str(p, N);

	return pitch_track(Sc, in_str, out_str, pm, pp);
}

/// Отслеживание ЧОТ в файлах.

unsigned spl_pitch_binary_file(int K, int N, 
	const freq_t *sc, const char *in_file, const char *out_file, 
	bool input_in_bits, double ksi)
{
	freq_scale Sc; Sc.K = K; Sc.Fr = const_cast<freq_t *>(sc);
	mask_parameters pm = DEFAULT_MASK_PARAMETERS; pm.ksi = ksi;
	pitch_parameters pp = DEFAULT_PITCH_PARAMETERS;

	// инициализируем входной поток из памяти
	io::ifstream<unsigned char> input(in_file); // real input
	io::iwrapelem<mask_t, unsigned char> in_bytes(input); // byte wrapper
	io::ibitwrap8 in_bits(input); // bit wrapper
	// выбираем вход:
	io::istream<mask_t> *in_str;
	if(input_in_bits) {
		in_str = &in_bits;
	} else {
		in_str = &in_bytes;
	}

	// инициализируем выходной поток в память
	io::ofstream<short> out_str(out_file);

	return pitch_track(Sc, *in_str, out_str, pm, pp);
}

unsigned spl_pitch_binary_file_ext(int K, int N, 
	const freq_t *sc, const char *in_file, const char *out_file, bool input_in_bits, 
	double ksi, double delta, double rho, bool border_effect, double F1, double F2, int Nh)
{
	freq_scale Sc; Sc.K = K; Sc.Fr = const_cast<freq_t *>(sc);
	mask_parameters pm; pm.ksi = ksi; pm.delta = delta; pm.rho = rho; pm.border_effect = border_effect;
	pitch_parameters pp; pp.Nh = Nh; pp.F1 = F1; pp.F2 = F2;

	// инициализируем входной поток из памяти
	io::ifstream<unsigned char> input(in_file); // real input
	io::iwrapelem<mask_t, unsigned char> in_bytes(input); // byte wrapper
	io::ibitwrap8 in_bits(input); // bit wrapper
	// выбираем вход:
	io::istream<mask_t> *in_str;
	if(input_in_bits) {
		in_str = &in_bits;
	} else {
		in_str = &in_bytes;
	}

	// инициализируем выходной поток в память
	io::ofstream<short> out_str(out_file);

	return pitch_track(Sc, *in_str, out_str, pm, pp);
}

} // extern "C"

