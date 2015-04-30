#ifndef _SPL_MASK_HELPER_
#define _SPL_MASK_HELPER_

#include "iofile.h"
#include "mask.h"
#include "conv.h"
#include "iowrap.h"
#include "spl_types.h"
#include <algorithm>

namespace spl {

//
// Forward declarations
//

/// Коэффициенты фильтров маскировки
struct mask_filters;

/// Вычисление коэффициентов одновременной маскировки по шкале резонансных частот фильтров.
bool mask_filters_generate(const freq_scale& s, mask_filters& f, const mask_parameters& p);

/// Вычисление коэффициентов одновременной маскировки, оптимизированных для быстрого вычисления.
bool mask_filters_generate_fast(const freq_scale& s, mask_filters& f, const mask_parameters& p);

/// Потоковая одновременная маскировка.
size_t EXPORT mask_stream(IN mask_filters& f, io::istream<spectrum_t>& in_str, io::ostream<mask_t>& out_str);

/// Потоковая одновременная маскировка - быстрый вариант.
size_t EXPORT mask_stream_fast(IN mask_filters& f, io::istream<spectrum_t>& in_str, io::ostream<mask_t>& out_str);



// Класс input-wrapper для расширения шкалы частот.
template<typename T>
class istream_block_extend;

// Класс output-wrapper для сужения шкалы частот.
template<typename T> 
class ostream_block_reduce;

/// Препроцессинг для расширения шкалы частот в спектре.
typedef istream_block_extend<spectrum_t> istream_spectrum_extend;

/// Постпроцессинг для сужения шкалы частот в маске.
typedef ostream_block_reduce<mask_t> ostream_mask_reduce;


/// Коэффициенты фильтров маскировки

struct mask_filters {
	
	/// Количество каналов в спектре.
	int K; 

	/// Размер окон фильтров.
	int Ws; 

	/// Массив коэффициентов фильтров маскировки.
	/// Матрица размерности (K, Ws).
	double *H;

	/// Форма шкалы, по которой ведется экстраполяция шкалы частот.
	scale_form sc_form;
};

///
/// то же, что mask_filters, 
/// только с управлением памятью (для С++)
///

struct EXPORT mask_filters_class: mask_filters 
{
	/// Конструктор.
	/// Делает вызов функции mask_filters_generate(), 
	///  предварительно выделив место для фильтра.
	/// При разрушении место освобождается.
	mask_filters_class(
		const freq_scale& s, 
		const mask_parameters& p);
	~mask_filters_class() { conv_free(H); }

	/// Сохранить коэффициенты в файл.
	bool save(const char *file) { return io::array_to_file(H, K * Ws, file); }

	/// Потоковая маскировка в зависимости от формы шкалы.
	/// Для модельной шкалы запускает быструю процедуру вычисления маски.
	/// Для остальных форм шкалы запускает обычную, неоптимизированную версию.
	size_t mask_stream(io::istream<spectrum_t>& in_str, io::ostream<mask_t>& out_str) {
		if(allow_fast) 
			return mask_stream_fast(*this, in_str, out_str);
		else 
			return spl::mask_stream(*this, in_str, out_str);
	}
private:
	bool allow_fast;
};



///
/// Препроцессинг потока: расширение блоков данных.
///
/// Это позволяет реализовывать потоковую одновременную маскировку
///  за счет расширения области спектра в обе стороны на половину размера окна.
///

template<typename T>
class istream_block_extend: public io::iwrap<T> 
{
public:
	/// Конструктор.
	/// Принимает оборачиваемый поток, размер блока и размер расширений.
	istream_block_extend(io::istream<T>& str, int _K, int _W):
	  io::iwrap<T>(str), K(_K), W(_W), _pos(0), _maxpos(0) { 
		// ставя _maxpos = _pos = 0, мы заставляем загрузить первый блок
		_buf = spl_alloc<T>(W + K + W);
	}

	/// Деструктор.
	/// Освобождает выделенную память
	~istream_block_extend() { spl_free(_buf); }

	//@{
	/// Стандартные функции abstract_stream.
	virtual bool eos() const { return _pos >= _maxpos && _understream->eos(); }
	virtual size_t pos() const { return _understream->pos() / K * (W+K+W) + _pos; }
	virtual size_t pos(size_t newsize) { return 0; } // пока не реализован за ненадобностью
	//@}

	/// Получение элемента.
	virtual bool get(T& x) {
		if(_pos >= _maxpos) {
			if(!read_block()) return false;
		}
		x = _buf[_pos++];
		return true;
	}
//*
	/// Получение строки элементов.
	/// Реализация функции istream::read()
	virtual size_t read(T *buf, size_t count) {
		size_t i = 0;
		// пока есть куда читать
		while(i < count) {
			// если элементы кончились - грузим следующий блок
			if(_pos >= _maxpos) {
				if(!read_block()) break; // если элементы и там кончились - выходим
			}
			// копируем минимальное из: оставшегося в _buf и оставшегося в buf
			size_t N = std::min(_maxpos - _pos, count - i);
			size_t r = std::copy(_buf + _pos, _buf + _pos + N, buf + i) - buf - i;
			_pos += r; i += r;
		}
		return i; // возвращаем количество прочтенных символов
	}
//*/
private:

	bool read_block() {
		// читаем блок: K элементов - запоминаем, сколько прочиталось
		_maxpos = _understream->read(_buf + W, K);
		// если ничего не прочиталось - конец файла
		if(_maxpos == 0) return false;
		// прочиталось _maxpos элементов, делаем расширение области:
		std::fill(_buf, _buf + W, _buf[W]);
		T *buf2 = _buf + W + _maxpos;
		std::fill(buf2, buf2 + W, buf2[-1]);
		_pos = 0;
		_maxpos = W + _maxpos + W;
		return true;
	}

	/// Размер блока.
	const int K;
	/// Размер расширения.
	const int W;
	/// Внутренний буфер с расширенной областью.
	T *_buf;
	/// Позиция в буфере.
	size_t _pos;
	/// Максимальная позиция в буфере.
	size_t _maxpos;
};

///
/// Постпроцессинг потока: отбрасывание отрезков.
/// Используется в быстром алгоритме одновременной маскировки
///  для перехода от расширенной шкалы частот к обычной.
///

template<typename T>
class ostream_block_reduce: public io::owrap<T>
{
public:
	/// Конструктор.
	/// Принимает выходной поток, размер блока и размер расширений.
	ostream_block_reduce(io::ostream<T>& str, size_t _K, size_t _W):
	  io::owrap<T>(str), K(_K), W(_W), rel_pos(0) {}

	/// Вывести элемент.
	bool put(const T& x) {
		bool r = true;
		if(rel_pos >= W && rel_pos < W + K) {
			r = _understream->put(x);
		}
		rel_pos = (rel_pos + 1) % (W + K + W);
		return r;
	}
//*
	/// Вывести строку элементов.
	size_t write(const T *buf, size_t count) {
		size_t i = 0, j = 0;
		while(i < count) {
			// если позиция не там, где нужно
			size_t N;
			if(rel_pos < W) {
				N = std::min(W - rel_pos, count - i);
			} else if(rel_pos >= K+W) {
				N = std::min(W+K+W+W - rel_pos, count - i);
			} else {
				N = 0;
			}
			i += N;
			rel_pos = (rel_pos + N) % (W + K + W);
			if(i >= count) break;
			N = std::min(W + K - rel_pos, count - i);
			size_t M = _understream->write(buf + i, N);
			i += M; j += M;
			// если вывелось меньше, чем нужно - выходим
			if(M < N) break;
			rel_pos += M;
		}
		return j; // вернуть число записанных на выход элементов.
	}
//*/
private:
	size_t rel_pos;
	/// Размер блока.
	const size_t K;
	/// Размер расширения.
	const size_t W;
};


} // namespace spl

#endif//_SPL_MASK_HELPER_