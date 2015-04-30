#ifndef _SPL_SPECTRUM_HELPER_
#define _SPL_SPECTRUM_HELPER_

#include "spectrum.h"
#include "iowrap.h"
#include <algorithm>

namespace spl {

///
/// Параметры, используемые для генерации фильтров (\ref spectrum_filters).
///

struct filter_parameters {
	
	/// Частота дискретизации сигнала.
	freq_t F;

	/// Величина, определяющая точность вычислений (размер окна фильтров).
	double ksi;

};

/// 
/// Коэффициенты нерекурсивного фильтра, используемого для получения спектра.
/// 

struct spectrum_filters {

	/// Количество каналов фильтрации.
	int K;
	
	/// Размер окон фильтрации.
	/// Размер окна фильтрации определяется при генерации фильтра и зависит от коэффициента точности ksi и частоты дискретизации.
	/// Чем больше окно фильтрации, тем медленнее работает фильтр, но тем точнее его вычисления.
	size_t Ws;

	/// Массив косинусных и синусных коэффициентов фильтрации.
	/// Матрица размерности (K, CONV_WIN_SIZ, 2)
	/// Синусные и косинусные коэффициенты фильтрации непосредственно используются для фильтрации.
	/// Поскольку синусные и косинусные коэффициенты отстоят по фазе на pi/2, 
	/// то будут учтены все точки в окне фильтрации, 
	/// несмотря на то что отдельные коэффициенты могут обращаться в ноль.
	double *H;

};

/// 
/// Генерация коэффициентов фильтрации по шкале резонансных частот фильтров.
/// В структуре фильтров \a f уже должно быть выделено место для помещения коэффициентов фильтрации.
/// 

bool spectrum_filters_generate(
	IN freq_scale& s,        ///< шкала резонансных частот фильтров
	OUT spectrum_filters& f, ///< коэффициенты фильтрации
	IN filter_parameters& p  ///< параметры генерации фильтров
);


///
/// то же, что spectrum_filters, 
/// только с управлением памятью (для С++)
///

struct EXPORT spectrum_filters_class: spectrum_filters 
{
	/// Конструктор.
	/// Делает вызов функции filters_generate(), 
	///  предварительно выделив место для фильтра.
	/// При разрушении место освобождается.
	spectrum_filters_class(
		const freq_scale& s, 
		const filter_parameters& p);
	/// Сохранить коэффициенты фильтров в файл.
	bool save(const char *file);
	~spectrum_filters_class();
};

///
/// Обертка, продлевающая поток сигнала на \a N нулевых отсчетов
///

template<typename T>
class iwstream_extend: public io::iwrap<T> {
public:

	iwstream_extend(io::istream<T>& str, size_t N):
	  io::iwrap<T>(str), _N(N), _pos(0) {}

	//@{
	/// Стандартные функции abstract_stream.
	virtual bool eos() const { return _pos >= _N; }
	virtual size_t pos() const { return _understream->pos() + _pos; }
	virtual size_t pos(size_t newpos) { return 0; } // пока не реализован за ненадобностью
	//@}

	/// Получение элемента 
	virtual bool get(T& x) {
		if(_understream->eos()) {
			if(!eos()) {
				x = 0;
				_pos++;
			} else {
				return false;
			}
		} else {
			_understream->get(x);
		}
		return true;
	}

	/// Получение строки элементов
	virtual size_t read(T *buf, size_t count) {
		size_t r1 = _understream->read(buf, count);
		size_t r2 = 0;
		if(r1 < count) {
			r2 = std::min(count-r1, _N-_pos);
			std::fill(buf + r1, buf + r1 + r2, 0);
			_pos += r2;
		}
		return r1 + r2;
	}

private:
	size_t _pos, _N;

};

} // namespace spl

#endif//_SPL_SPECTRUM_HELPER_