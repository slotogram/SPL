///
/// \file  mask.cpp
/// \brief Функции для расчета одновременной маскировки (\ref mask).
///
/// Модуль содержит функционал для генерации фильтров маскировки (функция mask_filters_generate())
///  и собственно для маскировки (свертки) - функции mask_stream(), mask_memory(), mask_binary_file().
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

/// Функция расчета окна фильтров маскировки.
/// Применяется для предварительного расчета - перед выделением памяти.

static int mask_window_size(const freq_scale& s, scale_form form, const mask_parameters& p) {
	// выбираем наибольшее из окон
	int Ws = 0;

	double Nmin, Nmax;
	for(int k = 0; k < s.K; k++) {
		double std = mask_std(s.Fr[k], p.delta);
		double dF = gauss_border(std, p.ksi);
		if(p.border_effect) { 
			Nmin = scale_interp_f2k(s, form, s.Fr[k] - dF); // если краевой эффект учитывается, 
			Nmax = scale_interp_f2k(s, form, s.Fr[k] + dF);	// то используется экстраполяция частот
		} else {
			Nmin = scale_find_freq(s, s.Fr[k] - dF); // если краевой эффект не учитывается, 
			Nmax = scale_find_freq(s, s.Fr[k] + dF); // то используется простой поиск по частотам
		}
		if(Nmax - Nmin + 1 > Ws) Ws = (int) ceil(Nmax - Nmin + 1);
	}
	return Ws % 2 ? Ws : Ws + 1; // должно быть обязательно нечетное число - чтобы был интервал [-Ws;+Ws]
}

/// Функция генерации маскирующей функции.
static void mask_win(
	const freq_scale& s, // массив резонансных частот фильтров
	double *H,           // выходной массив - коэффициенты маскирующей функции
	int k, 
	int Ws, 
	enum scale_form form, 
	const mask_parameters& p
) {
	int K = s.K;
	freq_t *Fr = s.Fr;
	double *H2 = H + Ws;

	// стандартное отклонение (СКО) функции Гаусса
	double std = mask_std(Fr[k], p.delta);
	// посчитать коэффициенты
	for(int m = -Ws; m <= Ws; m++) {
		double Hkm;
		if(k + m >= 0 && k + m < K) { // точка k + m есть в шкале частот?
			Hkm = gauss_win(Fr[k+m], Fr[k], std); // да - берем из шкалы частот
		} else if(p.border_effect) { // если учитывается краевой эффект, 
			double Fm = scale_interp_k2f(s, form, k+m); // то делаем экстраполяцию
			Hkm = gauss_win(Fm, Fr[k], std);
		} else { // если краевой эффект не учитывается,
			Hkm = 0; // то просто ставим 0 в маскирующую функцию
		}
		H2[m] = Hkm;
	}

	// нормирование
	double sum = 0;
	for(int m = 2*Ws; m >= 0; m--) // цикл от 0 до 2*Ws ВКЛЮЧИТЕЛЬНО, т.к. кол-во элементов = 2*Ws+1
		sum += H[m];
	sum *= p.rho;
	for(int m = 2*Ws; m >= 0; m--)
		H[m] /= sum;
}


/// 
/// Вычисление коэффициентов одновременной маскировки по шкале резонансных частот фильтров.
/// В структуре фильтров \a f уже должно быть выделено место (Ws*K) для помещения коэффициентов фильтрации.
///
/// Позволяет сгенерировать коэффициенты одновременной маскировки с учетом или без учета краевого эффекта.
/// Поскольку коэффициенты одновременной маскировки рассчитываются 
///  для каждого канала шкалы частот \a s с учетом частот соседних каналов, 
///  то рассчитать коэффициенты на крайних (первых и последних) каналах шкалы \a s 
///  не представляется возможным. 
/// Для расчета таких частот производится интерполяция шкалы частот \a s с помощью функции scale_interp_k2f().
/// Для того чтобы интерполяция была возможна, предварительно угадывается 
///  форма шкалы (\ref scale_form) с помощью функции scale_form_estimate().
/// Поэтому для возможности использовать одновременную маскировку с учетом краевого эффекта, 
///  необходимо использовать шкалы частот известной формы (см. \ref scale_form).
///
/// При использовании шкалы частот неизвестной формы 
///  можно генерировать только коэффициенты без учета краевого эффекта.
///
/// Учет краевого эффекта регулируется параметром \a p.border_effect.
///

bool mask_filters_generate(
	const freq_scale& s, 
	mask_filters& f, 
	const mask_parameters& p) 
{
	if(f.K != s.K) return false;

	// размер окна должен быть нечетным
	// на этом основан алгоритм - он движется от -Ws/2 до Ws/2
	int Ws = f.Ws/2; // половина длины окна - для удобства
	if(f.Ws != 2*Ws+1) return false;

	// посчитать коэффициенты
	for(int k = 0; k < f.K; k++) {
		mask_win(s, f.H + k * f.Ws, k, Ws, f.sc_form, p);
	}

	return true;
}

///
/// Вычисление коэффициентов одновременной маскировки, 
///  оптимизированных для быстрого вычисления.
///
/// Быстрое вычисление одновременной маскировки возможно, 
///  если шкала частот имеет форму scale_model (см. \ref scale_form).
/// В этом случае коэффициенты одновременной маскировки на всех частотных каналах одинаковы, 
///  что позволяет свести одновременную маскировку к цифровой фильтрации и 
///  вычислять ее по теореме о свертке.
///
/// В данной функции производится предварительное приготовление к быстрому вычислению:
///  вычисляется Фурье (вектор А) от маскирующей функции.
/// 
/// В структуре фильтров \a f уже должно быть выделено место (2*CONV_WIN_SIZ) 
///  для помещения коэффициентов фильтрации.
/// 

bool mask_filters_generate_fast(
	const freq_scale& s, 
	mask_filters& f, 
	const mask_parameters& p) 
{
	if(f.K != s.K) return false;
	
	// размер окна должен быть нечетным
	// на этом основан алгоритм - он движется от -Ws/2 до Ws/2
	int Ws = f.Ws/2; // половина длины окна - для удобства
	if(f.Ws != 2*Ws+1) return false;

	double H[CONV_WIN_SIZ];
	// вычисляет маскирующую функцию для k = K/2
	// выбор конкретного k на самом деле неважен
	mask_win(s, H, s.K/2, Ws, scale_model, p); // заполняет первые f.Ws байт
	// меняем направление - т.к. свертка поменяет его еще раз ;)
	std::reverse(H, H + f.Ws);
	// заполняем остаток нулями
	std::fill(H + f.Ws, H + CONV_WIN_SIZ, 0.0);
	// предвычисление вектора A
	cconv_calc_A(H, f.H, f.H + CONV_WIN_SIZ);
	// нормировка вектора А
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

	// для модельной шкалы - быстрое вычисление - другой размер
	size_t N = allow_fast ? 2 * CONV_WIN_SIZ : K * Ws;
	H = conv_alloc<double>(N);
	if(H == 0) 
		throw "Can't allocate memory for mask filters coefficients";

	typedef bool (*generate_func)(const freq_scale&, mask_filters&, const mask_parameters&);
	// для модельной шкалы - быстрое вычисление
	generate_func generate = allow_fast ? mask_filters_generate_fast : mask_filters_generate;

	if(!generate(s, *this, p)) 
		throw "Error while generating mask filters";
}


///
/// Потоковая одновременная маскировка.
/// Возвращает количество записанных на выход элементов.
///
/// Является наиболее общей функцией для вычисления одновременной маскировки, 
///  воспринимающей вход и выход как абстрактные потоки.
/// Более детальная реализация выбирает конкретную реализацию потоков 
///  (память, файлы, микрофон, колонки и т.п.),
///  инициализирует входные структуры фильтров и потоков 
///  и вызывает данную функцию.
/// Это обеспечивает единство вычислений, вне зависимости от используемых типов потоков.
/// 
/// Функция работает по наивному алгоритму, без оптимизаций.
///

size_t mask_stream(
	IN mask_filters& f,          ///< Коэффициенты маскировки
	istream<spectrum_t>& in_str, ///< Входной поток  (спектр)
	ostream<mask_t>& out_str     ///< Выходной поток (маска)
) {
	int K = f.K;
	int Ws = f.Ws/2;

	spectrum_t *spec_wide1 = spl_alloc<spectrum_t>(Ws + K + Ws);
	spectrum_t *spec_input = spec_wide1 + Ws;
	spectrum_t *spec_wide2 = spec_input + K;
	if(!spec_wide1) return 0;

	size_t written = 0;

	// основной цикл маскировки
	while(in_str.read(spec_input, K) == K && !out_str.eos()) {

		// расширение области частот спектра
		std::fill(spec_wide1, spec_wide1 + Ws, spec_input[0]);
		std::fill(spec_wide2, spec_wide2 + Ws, spec_input[K-1]);

		// цикл расчета маскировки:
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
/// Потоковая одновременная маскировка - быстрый вариант.
/// Вызывается только для модельной шкалы частот.
///

size_t mask_stream_fast(
	IN mask_filters& f,          ///< Коэффициенты маскировки
	istream<spectrum_t>& in_str, ///< Входной поток  (спектр)
	ostream<mask_t>& out_str     ///< Выходной поток (маска)
) {
	int K = f.K;
	int Ws = f.Ws/2;
	int Os = CONV_WIN_SIZ - 2*Ws;

	spectrum_t array[CONV_WIN_SIZ * 6 + 2];
	// Два буфера для чтения спектра
	spectrum_t *input_buf1 = SPL_MEMORY_ALIGN(array);
	spectrum_t *input_buf2 = input_buf1 + CONV_WIN_SIZ;
	// Четыре временных буфера - B, C, abcr, abci
	spectrum_t *tmp_buf1 = input_buf2 + CONV_WIN_SIZ;
	spectrum_t *tmp_buf2 = tmp_buf1 + CONV_WIN_SIZ;
	spectrum_t *tmp_buf3 = tmp_buf2 + CONV_WIN_SIZ;
	spectrum_t *tmp_buf4 = tmp_buf3 + CONV_WIN_SIZ;
	// один выходной буфер
	mask_t out_buf[CONV_WIN_SIZ * 2];

	// i/o wrappers для расширения/сужения шкалы частот:
	istream_spectrum_extend in_wrap(in_str, K, Ws);
	ostream_mask_reduce out_wrap(out_str, K, Ws);

	// выходная величина
	size_t written = 0;

	// заполняем конец второго буфера, чтобы потом оттуда перенеслось в начало первого
	std::fill(input_buf2 + Os, input_buf2 + Os + Ws, 0);
	in_wrap.read(input_buf2 + Os + Ws, Ws);

	size_t N1 = 0, N2 = 0;

	// основной цикл маскировки
	while(!in_wrap.eos() && !out_wrap.eos() || N2 > Os - Ws) { // последний блок обрабатывается только на следующей итерации

		// копируем из конца второго буфера в начало первого
		std::copy(input_buf2 + Os, input_buf2 + CONV_WIN_SIZ, input_buf1);
		// читаем первый буфер 
		N1 = in_wrap.read(input_buf1 + 2*Ws, Os);
		// если прочиталось меньше, чем нужно - остаток заполняем нулями
		if(N1 < Os) {
			std::fill(input_buf1 + 2*Ws + N1, input_buf1 + CONV_WIN_SIZ, 0);
		}

		// копируем из конца первого буфера в начало второго
		std::copy(input_buf1 + Os, input_buf1 + CONV_WIN_SIZ, input_buf2);
		// читаем второй буфер
		N2 = in_wrap.read(input_buf2 + 2*Ws, Os);
		// если прочиталось меньше, чем нужно - остаток заполняем нулями
		if(N2 < Os) {
			std::fill(input_buf2 + 2*Ws + N2, input_buf2 + CONV_WIN_SIZ, 0);
		}

		// считаем B, C
		cconv_calc_BC(input_buf1, input_buf2, tmp_buf1, tmp_buf2);

		// свертка
		cconv(f.H, f.H + CONV_WIN_SIZ, tmp_buf1, tmp_buf2, tmp_buf3, tmp_buf4);

		int j = 0;
		// вычисляем результат маскировки для обоих буферов:
		for(size_t i = Ws; i < Ws + N1; i++) {
			out_buf[j++] = (input_buf1[i] > tmp_buf3[i + Ws]);
		}
		for(size_t i = Ws; i < Ws + N2; i++) {
			out_buf[j++] = (input_buf2[i] > tmp_buf4[i + Ws]);
		}

		// выводим результат
		written += out_wrap.write(out_buf, j);
	}
	return written;
}


/// Маскировка потока спектра по шкале частот.

size_t mask_stream(
	IN freq_scale& sc,           ///< Шкала резонансных частот
	istream<spectrum_t>& in_str, ///< Входной поток  (спектр)
	ostream<mask_t>& out_str,    ///< Выходной поток (маска)
	IN mask_parameters& p        ///< Параметры маскировки
) {
	try {
		// фильтры
		mask_filters_class f(sc, p);

		// собственно маскировка
		// в зависимости от формы шкалы выбирает быструю или обычную версию маскировки
		return f.mask_stream(in_str, out_str);

	} catch(const char *e) {
		fprintf(stderr, "%s", e);
		return 0;
	}
}

/// Маскировка потока спектра со стандартными параметрами генерации маскирующей функции.

size_t mask_stream(
	IN freq_scale& sc,           ///< Шкала резонансных частот
	istream<spectrum_t>& in_str, ///< Входной поток  (спектр)
	ostream<mask_t>& out_str,    ///< Выходной поток (маска)
	IN double ksi                ///< Величина, определяющая точность вычислений (размер окна фильтров)
) {
	mask_parameters p = DEFAULT_MASK_PARAMETERS;
	p.ksi = ksi;
	return mask_stream(sc, in_str, out_str, p);
}

/// Маскировка в памяти

size_t mask_memory(
	IN freq_scale& sc, 
	IN spectrum& sp, 
	OUT mask& m, 
	IN mask_parameters& p
) {
	// проверяем размер выделенной области
	if(sp.K != m.K || sp.N != m.N) return 0;
	int NK = m.N * m.K;

	// инициализируем входной поток из памяти
	imstream<spectrum_t> in_str(sp.Y, NK);

	// инициализируем выходной поток в память
	omstream<mask_t> out_str(m.Z, NK);

	// собственно маскировка
	return mask_stream(sc, in_str, out_str, p);
}

/// Маскировка в памяти со стандартными параметрами.

size_t mask_memory(
	IN freq_scale& sc,    ///< Шкала резонансных частот спектра
	IN spectrum& sp,      ///< Спектр (входной)
	OUT mask& m,          ///< Маска (выходная)
	IN double ksi         ///< Параметры маскировки
) {
	mask_parameters p = DEFAULT_MASK_PARAMETERS;
	p.ksi = ksi;
	return mask_memory(sc, sp, m, p);
}

/// Маскировка бинарного файла.

size_t mask_binary_file(
	IN freq_scale& sc, 
	IN char *in_file, 
	IN char *out_file, 
	IN bool out_bits, 
	IN mask_parameters& p
) {
	// инициализируем входной поток из файла
	ifstream<spectrum_t> in_str(in_file);

	// выходное значение
	size_t r;

	if(out_bits) { // если биты
		// инициализируем выходной поток в файл
		ofstream<unsigned char> out_str(out_file);
		io::obitwrap8 out_bit_str(out_str);
		// собственно фильтрация
		r = mask_stream(sc, in_str, out_bit_str, p);
		out_bit_str.flush(); // запись последнего блока (байта)
	} else { // если байты
		// инициализируем выходной поток в файл
		ofstream<mask_t> out_str(out_file);
		// собственно фильтрация
		r = mask_stream(sc, in_str, out_str, p);
	}
	return r;
}

/// Маскировка бинарного файла со стандартными параметрами.

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
