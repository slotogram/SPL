///
/// \file  spectrum.cpp
/// \brief Модуль выделения спектра звукового сигнала ("фильтрации").
///
/// Выделение спектра реализовано через систему нерекурсивных фильтров. 
/// Модуль содержит функционал для генерации параметров данных фильтров и собственно для фильтрации (свертки).
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
/// Генерация коэффициентов фильтрации по шкале резонансных частот фильтров.
/// В структуре фильтров \a f уже должно быть выделено место для помещения коэффициентов фильтрации.
/// 

bool spectrum_filters_generate(
	IN freq_scale& s,        ///< шкала резонансных частот фильтров
	OUT spectrum_filters& f, ///< коэффициенты фильтрации
	IN filter_parameters& p  ///< параметры генерации фильтров
) {
	int n, j, k; // итераторы
	int K = s.K;
	freq_t F = p.F;

	double array[CONV_WIN_SIZ * 2 + 2];
	double *Hc = SPL_MEMORY_ALIGN(array);
	double *Hs = Hc + CONV_WIN_SIZ;

	if(s.K != f.K) return false;

	// вычисляем Ws - размер реального окна фильтров
	// выбирается как максимальное из вычисленных для каждого канала
	int Ws = 0;
	for(k = 0; k < K; k++) {
		const double std = filter_std(s.Fr[k], F);
		int Wsk = (int) gauss_border(std, p.ksi);
		if(Wsk > Ws) Ws = Wsk; // на каком канале это произошло, нам даже не важно
	}
	f.Ws = 2*Ws+1;
	
	// матрица - для удобного доступа к коэффициентам фильтрации
	Matrix<double, 3> H = matrix_ptr(f.H, 2, f.K, CONV_WIN_SIZ);

	// вычисляем собственно коэффициенты фильтра
	// для размера окна Ws
	for(k = 0; k < K; k++) {
		const double std = filter_std(s.Fr[k], F);
		const double Wf = 2 * M_PI * s.Fr[k] / F;
		for(n = -Ws, j = 0; n <= Ws; n++, j++) {
			double norm = gauss_win<double>(n, 0.0, std); // окно Гаусса нормирует синусоиды
			Hc[j] = norm * cos(Wf * n); // действительная часть
			Hs[j] = norm * sin(Wf * n); // мнимая часть
		}
		// заполняем оставшиеся коэффициенты нулями
		std::fill(Hc + j, Hc + CONV_WIN_SIZ, 0);
		std::fill(Hs + j, Hs + CONV_WIN_SIZ, 0);
		// предварительное вычисление вектора BC
		cconv_calc_BC(Hc, Hs, &H(0,k,0), &H(1,k,0));
	}
	// нормализация коэффициентов фильтрации:
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
/// Фильтрация потока с заданным набором фильтров.
/// Возвращает количество записанных на выход элементов.
/// 
/// Является наиболее общей функцией для фильтрации, 
///  воспринимающей вход и выход, как абстрактные потоки.
/// Более детальная реализация выбирает конкретную реализацию потоков 
///  (память, файлы, микрофон, колонки и т.п.),
///  инициализирует входные структуры фильтров и потоков 
///  и вызывает данную функцию.
/// Это обеспечивает единство вычислений, вне зависимости от используемых типов потоков.
///
/// Данная функция использует оптимизацию вычисления свертки через FFT.
///

size_t filter_stream(
	IN spectrum_filters& f,       ///< Коэффициенты фильтрации
	istream<signal_t>& in_str,    ///< Входной поток  (сигнал)
	ostream<spectrum_t>& out_str  ///< Выходной поток (спектр)
) {
	int K = f.K;
	size_t Ws = f.Ws - 1; // можно брать на 1 меньше, чем окно - результат не меняется
	size_t Os = CONV_WIN_SIZ - Ws;
	spectrum_t *out_buf = NULL;

	// обеспечиваем отсутствие смещения в начале сигнала
	iwstream_extend<signal_t> in_wrap(in_str, Ws/2);

	// сколько записано - выходная величина
	size_t written = 0;

	double array[5*CONV_WIN_SIZ + 2];

	// буфер входного сигнала - состоит из двух частей:
	// CONV_WIN_SIZ = Ws + Os, 
	// где Ws - Window size - размер реального окна фильтра, 
	//     CONV_WIN_SIZ - размер вычисляемой циклической свертки
	//     Os - Output size - размер полезного выхода свертки
	// также используется как входной буфер свертки
	double *conv_in_buf = SPL_MEMORY_ALIGN(array);

	// выходной буфер свертки
	// имеет такую же структуру как и входной буфер (2*Ws+1) + Os
	// только полезный выход - последние Os элементов - идут на выход
	double *tmp_buf1 = conv_in_buf + CONV_WIN_SIZ;
	double *tmp_buf2 = tmp_buf1 + CONV_WIN_SIZ;
	double *tmp_buf3 = tmp_buf2 + CONV_WIN_SIZ;
	double *tmp_buf4 = tmp_buf3 + CONV_WIN_SIZ;

	// буфер выходного сигнала
	// матрицы размера K x Os - для помещения результата свертки
	// и Os x K - для вывода наружу
	out_buf = spl_alloc<spectrum_t>(2 * K * Os);

	Matrix<spectrum_t,2> out_mtx1 = matrix_ptr(out_buf, K, Os);
	Matrix<spectrum_t,2> out_mtx2 = matrix_ptr(out_buf+K*Os, Os, K);

	// матрица - для удобного доступа к коэффициентам фильтрации
	Matrix<double, 3> H = matrix_ptr(f.H, 2, f.K, CONV_WIN_SIZ);

	if(!out_buf) 
		goto end; 

	//
	// подготовка структур данных
	//

	// очищаем последние Ws элементов входного буфера
	std::fill(conv_in_buf+CONV_WIN_SIZ-Ws, conv_in_buf+CONV_WIN_SIZ-Ws/2, 0);

	// обеспечиваем отсутствие смещения в начале сигнала
	in_wrap.read(conv_in_buf+CONV_WIN_SIZ-Ws/2, Ws/2);

	// основной цикл фильтрации
	while(!in_wrap.eos() && !out_str.eos()) {

		// копируем последние Ws элементов сигнала в начало
		std::copy(conv_in_buf+CONV_WIN_SIZ-Ws, conv_in_buf+CONV_WIN_SIZ, conv_in_buf);

		// вводим Os новых отсчетов сигнала
		// этот буфер будет использоваться неизменно для каждого канала
		// получаем rOs - real output size - количество считанных элементов
		// в общем случае rOs == Os, отличия могут быть только в конце сигнала
		size_t rOs = in_wrap.read(conv_in_buf + Ws, Os);

		// если считано меньше, чем Os элементов,
		//  то будет последний виток цикла
		// очищаем последние отсчеты сигнала, 
		// и заменяем матрицы на матрицы меньшего размера (K x rOs)
		if(rOs != Os) {
			std::fill(conv_in_buf + Ws + rOs, conv_in_buf + CONV_WIN_SIZ, 0);
			out_mtx1 = matrix_ptr(out_buf, K, rOs);
			out_mtx2 = matrix_ptr(out_buf+K*rOs, rOs, K);
		}

		spectrum_t *out_ptr = out_mtx1;

		cconv_calc_A(conv_in_buf, tmp_buf1, tmp_buf2);

		// цикл по каналам
		for(int k = 0; k < K; k++) {

			// свертка
			cconv(tmp_buf1, tmp_buf2, &H(0,k,0), &H(1,k,0), tmp_buf3, tmp_buf4);

			// вычисление модуля комплексных чисел
			complex_abs_split(rOs, tmp_buf3 + Ws, tmp_buf4 + Ws, tmp_buf4 + Ws);

			// пишем выход в выходную матрицу - только из интервала [Ws, Ws+rOs]
			out_ptr = std::copy(tmp_buf4 + Ws, tmp_buf4 + Ws + rOs, out_ptr);

		}

		// транспонируем выходную матрицу
		// из K x rOs в rOs x K
		transpose(out_mtx1, out_mtx2);

		// выводим матрицу out_mtx2
		written += out_str.write(out_mtx2, out_mtx2.size());

	}

end:
	// удаляем выделенные буферы
	spl_free(out_buf);

	return written;

}

///
/// Фильтрация потока сигнала по шкале частот.
/// 

size_t filter_stream(
	IN freq_scale& sc, ///< Шкала частот 
	istream<signal_t>& in_str, ///< Входной поток (сигнал)
	ostream<spectrum_t>& out_str, ///< Выходной поток (спектр)
	IN freq_t F, ///< Частота дискретизации
	IN double ksi) ///< Величина, определяющая точность вычислений (размер окна)
{
	// Параметры создания фильтров.
	filter_parameters p;
	p.F = F;
	p.ksi = ksi;

	try {

		// фильтры
		spectrum_filters_class f(sc, p);

		// собственно фильтрация
		return filter_stream(f, in_str, out_str);

	} catch(const char *e) {
		fprintf(stderr, "%s", e);
		return 0;
	}

}

///
/// Фильтрация в оперативной памяти.
///

size_t filter_memory(
	IN freq_scale& sc, ///< Шкала резонансных частот фильтров
	IN signal& s,      ///< Входной сигнал (в памяти)
	OUT spectrum& sp,  ///< Выходной спектр (в памяти)
	IN double ksi)     ///< Величина, определяющая точность вычислений (размер окна)
{
	// проверяем размер выделенной области
	if(sp.K != sc.K || sp.N != s.N) return 0;

	// ставим частоту дискретизации
	sp.F = s.F;
	
	// инициализируем входной поток из памяти
	imstream<signal_t> in_str(s.X, s.N);

	// инициализируем выходной поток в память
	omstream<spectrum_t> out_str(sp.Y, sp.N * sp.K);

	// собственно фильтрация
	return filter_stream(sc, in_str, out_str, s.F, ksi);

}

///
/// Фильтрация из wav в оперативную память.
///

size_t filter_wmemory(
	IN freq_scale& sc, ///< Шкала резонансных частот фильтров
	IN char *in_file,  ///< Входной сигнал (имя файла)
	OUT spectrum& sp,  ///< Выходной спектр (в памяти)
	IN double ksi,     ///< Величина, определяющая точность вычислений (размер окна)
	IN int k)			///< Каналы
{
	
	// инициализируем входной поток из файла
	iwstream in_str(in_file);
	
	//переводим в память
	
	spl::signal sig = in_str.getInMemory();

	//выделяем память на спектр
	//spectrum* spec = new spectrum;
	sp.F = sig.F;
	sp.N = sig.N;
	sp.K = k;
	sp.Y = new spectrum_t[sp.N*k];


	// собственно фильтрация
	size_t ret = filter_memory(sc,sig,sp,ksi);
	
	delete [] sig.X;
	//delete &sig;

	return ret;

}


///
/// Фильтрация сигнала как бинарного файла.
///

size_t filter_binary_file(
	IN freq_scale& sc, ///< Шкала резонансных частот фильтров
	IN char *in_file,  ///< Входной сигнал (имя файла)
	IN char *out_file, ///< Выходной спектр (имя файла)
	IN freq_t F,       ///< Частота дискретизации
	IN double ksi      ///< Величина, определяющая точность вычислений (размер окна)
) {

	// инициализируем входной поток из файла
	ifstream<signal_t> in_str(in_file);

	// инициализируем выходной поток в файл
	ofstream<spectrum_t> out_str(out_file);

	// собственно фильтрация
	return filter_stream(sc, in_str, out_str, F, ksi);

}

/// Фильтрация wave-файла.

size_t filter_wave_file(
	IN freq_scale& sc, ///< Шкала резонансных частот фильтров
	IN char *in_file,  ///< Входной сигнал (имя файла)
	IN char *out_file, ///< Выходной спектр (имя файла)
	IN double ksi      ///< Величина, определяющая точность вычислений (размер окна)
) {

	// инициализируем входной поток из файла
	iwstream in_str(in_file);

	// инициализируем выходной поток в файл
	ofstream<spectrum_t> out_str(out_file);

	// собственно фильтрация
	return filter_stream(sc, in_str, out_str, in_str.freq(), ksi);

}


} // namespace spl 
