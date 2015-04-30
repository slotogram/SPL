#ifndef _SPL_VOCAL_
#define _SPL_VOCAL_

///
/// \file  vocal.h
/// \brief Функции для обработки звуков по признаку вокализованности.
///
/// Включает в себя функции для определения частоты основного тона (ЧОТ) 
///  и автоматической сегментации сигналов на вокализованные/невокализованные.
///

#include "common.h"
#include "spl_types.h"
#include "mask.h"
#include "io.h"

namespace spl {

/// Шаблоны маскировки (используются для определения частоты основного тона).

struct mask_templates {

	//@{
	/// Первый и последний канал, для которых есть шаблоны.
	/// Структура \ref mask_templates содержит шаблоны для каналов k1 ... k2.
	/// Таким образом, в структуре \ref mask_templates содержится (k2 - k1 + 1) шаблонов.
	int k1, k2;

	/// Полное количество каналов.
	int K;

	/// Количество элементов шаблона.
	int Nt;

	/// Собственно шаблоны.
	/// Шаблон представляет собой матрицу размера (k2 - k1 + 1) x Nt.
	int *T;

};

/// Параметры генерации шаблонов для выделения ЧОТ - частоты основного тона (см. \ref mask_templates).

struct pitch_parameters {

	/// Количество учитываемых гармоник (для определения ЧОТ: 2).
	int Nh;

	//@{ 
	/// Нижняя и верхняя границы определения ЧОТ.
	freq_t F1, F2; //@}

};

/// Стандартные значения параметров генерации шаблонов, для определения ЧОТ.
const pitch_parameters DEFAULT_PITCH_PARAMETERS = { 2, 75, 400 };


/// Параметры сегментации по признаку вокализованности.

struct vocal_parameters {

	/// минимальная длительность вокализованного сегмента (в секундах)
	double minV;

	/// минимальная длительность невокализованного сегмента (в секундах)
	double minNV;

};

/// Стандартные значения параметров сегментации по признаку вокализованности.
const vocal_parameters DEFAULT_VOCAL_PARAMETERS = { 0.030, 0.030 };

/// Генерация шаблона маски для определения частоты основного тона.
bool EXPORT pitch_templates_generate(
	IN freq_scale& sc,        ///< Шкала частот.
	OUT mask_templates& tpl,  ///< Шаблон маски.
	IN mask_parameters& pm,   ///< Параметры генерации маскирующей функции.
	IN pitch_parameters& p    ///< Параметры генерации шаблона маски.
);

/// Класс шаблонов определения маски с автоматической очисткой памяти.
struct EXPORT pitch_templates_class: mask_templates {
	/// Конструктор. Вызывает функцию pitch_templates_generate
	pitch_templates_class(
		const freq_scale& sc, 
		const mask_parameters& pm, 
		const pitch_parameters& p) 
	{
		if(!pitch_templates_generate(sc, *this, pm, p)) {
			throw "pitch templates generation error";
		}
	}
	/// Деструктор. Очищает занимаемую память.
	~pitch_templates_class() { spl_free(T); }
};

/// Поиск в потоке маски по шаблону. 
/// На входе -- поток маски, на выходе - номера каналов ЧОТ для каждого отсчета.
/// Возвращает количество 
size_t EXPORT mask_template_search(
	IN mask_templates& tpl,       ///< шаблоны маски, по которым осуществляется поиск
	io::istream<mask_t>& in_str,  ///< входной поток (маска)
	io::ostream<short>& out_str,  ///< выходной поток (номера найденных каналов)
	IN int max_diff               ///< максимальное отличие маски от шаблона
);

size_t EXPORT mask_template_fast_search(
	IN mask_templates& tpl,       ///< шаблоны маски, по которым осуществляется поиск
	io::istream<mask_t>& in_str,  ///< входной поток (маска)
	io::ostream<short>& out_str,  ///< выходной поток (номера найденных каналов)
	IN int max_diff               ///< максимальное отличие маски от шаблона
);

/// Стандартное значение максимального отклонения маски вокализованного звука от шаблона.
const int DEFAULT_PITCH_MAX_DIFF = 6;

/// Функция определения частоты основного тона.
size_t EXPORT pitch_track(
	IN freq_scale& sc,           ///< шкала частот
	io::istream<mask_t>& in_str, ///< входной поток (маска)
	io::ostream<short>& out_str, ///< выходной поток (номера каналов ЧОТ)
	IN mask_parameters& pm,      ///< параметры генерации маскирующей функции
	IN pitch_parameters& p       ///< параметры генерации шаблона маски
);

/// Функция сегментации сигнала по признаку вокализованности
size_t EXPORT vocal_segment(
	io::istream<short>& in_str,  ///< входной поток (номера каналов ЧОТ)
	io::ostream<short>& out_str, ///< выходной поток (номера отсчетов-границ сегментов)
	IN freq_t F,                 ///< частота дискретизации сигнала
	IN vocal_parameters& p       ///< параметры сегментации
);

} // namespace spl

#endif//_SPL_VOCAL_