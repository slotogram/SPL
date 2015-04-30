#ifndef _IO_WAVE_
#define _IO_WAVE_

#include "iofile.h"
#include "spl_types.h"
#include <limits>

namespace spl {

/// 
/// Поток wave-файла.
/// Открывает wave-файл и представляет к нему интерфейс io::istream<signal_t>.
/// Дополнительно имеет функцию freq() для получения частоты дискретизации.
///

class EXPORT iwstream:
	public io::istream<signal_t>
{
public:
	/// Конструктор. Принимает путь к wav-файлу. 
	iwstream(const char *file);
	/// Получить очередной отсчет.
	/// В многоканальных сигналах (например, стерео) данная функция получает сумму по каналам за отсчет.
	virtual bool get(signal_t& x);
	/// Текущая позиция в потоке. Фактически означает время.
	virtual size_t pos() const;
	/// Установить позицию в потоке.
	virtual size_t pos(size_t N);
	/// Пропустить \a N отсчетов.
	virtual size_t skip(size_t N);
	/// Проверить конец потока.
	virtual bool eos() const;

	/// Получить частоту дискретизации.
	freq_t freq() const { return (freq_t)F; }
	//Загрузить сигнал в память.
	signal getInMemory();

	/// Закрыть поток.
	virtual void close() { _f.close(); }

private:
	typedef unsigned char byte;
	typedef unsigned short word;
	typedef unsigned long dword;
	io::ifstream<byte> _f;

	template<typename T>
	bool get_element(T& x) {
		return _f.read((byte*)&x, sizeof(T)) > 0;
	}

	// Функция, получающая один суммарный сэмпл данных по всем каналам.
	// Например для стерео - сумма по двум каналам.
	template<typename T>
	bool get_slice(T& x);
	
	// WAVE FILE FORMAT
	size_t _base; ///< база секции отсчетов
	dword N; ///< размер секции отсчетов
	word M; ///< количество чисел во временном отсчете (каналы - моно, стерео)
	word B; ///< количество байт в числе
	dword F; ///< частота дискретизации
};

template<typename T>
bool iwstream::get_slice(T& x) {
	x = 0; 
	T y;
	for(int i = 0; i < M; i++) {
		if(!get_element(y)) return false;
		x += y;
	}
	return true;
}

template<typename T>
void convert_sample_to_float(const T& x, signal_t& y) {
	T m = std::numeric_limits<T>::min();
	T M = std::numeric_limits<T>::max();
	y = ( 2*signal_t(x) - (signal_t(M) + 1 + m) ) / (signal_t(M) + 1 - m);
}

template<typename T>
void convert_float_to_sample(const signal_t& x, T& y) {
	// TODO implement this
}

} // namespace spl

#endif//_IO_WAVE_
