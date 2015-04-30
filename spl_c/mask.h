#ifndef _SPL_MASK_
#define _SPL_MASK_

///
/// \file  mask.h
/// \brief Типы данных, структуры и функции, относящиеся к маскировке.
///

#include "common.h"
#include "spl_types.h"
#include "io.h"

namespace spl {

/// Параметры, используемые для генерации коэффициентов маскировки (\ref mask_filters).

struct mask_parameters {

	/// Величина, определяющая точность вычислений (размер окна фильтров).
	double ksi;

	/// Коэффициент, определяющий ширину маскирующей функции.
	double delta;

	/// Коэффициент, определяющий вес маскирующей функции.
	double rho;

	/// Триггер - учитывать краевой эффект или нет (см. mask_filters_generate()).
	bool border_effect;

	/// Триггер - позволять быструю версию или нет.
	bool allow_fast;

};

/// Стандартные значения параметров генерации маскирующей функции.
const mask_parameters DEFAULT_MASK_PARAMETERS = { 0.001, 1, 1, true, true };

///
/// Маскировка потока спектра по шкале частот.
/// Возвращает количество записанных элементов маски.
/// 

size_t EXPORT mask_stream(
	IN freq_scale& sc,               ///< Шкала резонансных частот
	io::istream<spectrum_t>& in_str, ///< Входной поток  (спектр)
	io::ostream<mask_t>& out_str,    ///< Выходной поток (маска)
	IN mask_parameters& p            ///< Параметры маскировки
);

///
/// Маскировка в памяти
/// В поле \a m.Y память должна быть выделена перед вызовом, 
///  а также установлены поля \a m.K и \a m.N (см. \ref mask).
/// Функция проверяет правильность выделения памяти по этим полям.
///

size_t EXPORT mask_memory(
	IN freq_scale& sc,    ///< Шкала резонансных частот спектра
	IN spectrum& sp,      ///< Спектр (входной)
	OUT mask& m,          ///< Маска (выходная)
	IN mask_parameters& p ///< Параметры маскировки
);

/// Маскировка в памяти со стандартными параметрами.

size_t EXPORT mask_memory(
	IN freq_scale& sc,    ///< Шкала резонансных частот спектра
	IN spectrum& sp,      ///< Спектр (входной)
	OUT mask& m,          ///< Маска (выходная)
	IN double ksi         ///< Параметры маскировки
);

///
/// Маскировка бинарного файла.
///

size_t EXPORT mask_binary_file(
	IN freq_scale& sc,    ///< Шкала резонансных частот спектра.
	IN char *in_file,     ///< Входной файл (спектр)
	IN char *out_file,    ///< Выходной файл (маска)
	IN bool out_bits,     ///< Выводить битами (или байтами)
	IN mask_parameters& p ///< Параметры маскировки
);

/// Маскировка бинарного файла со стандартными параметрами.

size_t EXPORT mask_binary_file(
	IN freq_scale& sc,    ///< Шкала резонансных частот спектра.
	IN char *in_file,     ///< Входной файл (спектр)
	IN char *out_file,    ///< Выходной файл (маска)
	IN bool out_bits,     ///< Выводить битами (или байтами)
	IN double ksi         ///< Параметры маскировки
);

} // namespace spl

#endif//_SPL_MASK_
