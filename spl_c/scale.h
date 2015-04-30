#ifndef _SPL_SCALE
#define _SPL_SCALE

///
/// \file  scale.h
/// \brief Структуры и функции, относящиеся к частотным шкалам.
///

#include "common.h"
#include "spl_types.h"

namespace spl {

///
/// Параметры, используемые для генерации шкалы частот \ref freq_scale.
///

struct scale_parameters {

	//@{
	/// Номера первого и последнего канала фильтрации.
	/// Количество каналов K = b - a + 1;
	int a, b; //@}

	//@{
	/// Номера двух расчетных каналов фильтрации.
	int i, j; //@}

	//@{
	/// Две расчетные частоты, соответствующие расчетным каналам фильтрации i и j.
	freq_t Fi, Fj; //@}

};

/// Функция для сохранения шкалы в файл.
bool EXPORT save_scale(
	IN freq_scale& sc, ///< Шкала, которую необходимо сохранить.
	const char *file   ///< Путь к файлу, в котором необходимо сохранить шкалу.
);

//
// Функции генерации шкалы
// 

/// Функция шкалы - для любой заданной частоты возвращает значение некоторой метрики.
typedef double (*scale_func)(freq_t);

/// Обратная функция шкалы - для значения метрики возвращает частоту.
typedef freq_t (*rscale_func)(double);

/// Функция для генерации шкалы по заданным функциям шкалы. 
bool EXPORT scale_generate(
	OUT freq_scale& sc,     ///< Шкала, которая будет сгенерирована
	IN scale_parameters& p, ///< Параметры генерации шкалы
	IN scale_func f2x, ///< Прямая функция шкалы
	IN rscale_func x2f ///< Обратная функция шкалы
);

/// Форма генерируемой шкалы
enum scale_form {
	scale_unknown=0,///< неизвестная форма шкалы
	scale_linear,  ///< линейная шкала
	scale_log,     ///< логарифмическая шкала
	scale_mel,     ///< шкала мел
	scale_bark,    ///< шкала барк
	scale_model,   ///< шкала модели слуха
};

/// Функция для генерации шкалы по заданной форме.
bool EXPORT scale_generate(
	OUT freq_scale& sc,     ///< Шкала, которая будет сгенерирована
	IN scale_parameters& p, ///< Параметры генерации шкалы
	IN scale_form form      ///< Форма шкалы
);

/// Общее определение функции генерации шкалы.
typedef bool (*scale_generate_func)(OUT freq_scale& sc, IN scale_parameters& p);

/// Функция для генерации линейной шкалы.
bool EXPORT scale_generate_linear(OUT freq_scale& sc, IN scale_parameters& p);

/// Функция для генерации шкалы Bark.
bool EXPORT scale_generate_bark(OUT freq_scale& sc, IN scale_parameters& p);

/// Функция для генерации шкалы Mel.
bool EXPORT scale_generate_mel(OUT freq_scale& sc, IN scale_parameters& p);

/// Функция для генерации логарифмической шкалы частот.
bool EXPORT scale_generate_log(OUT freq_scale& sc, IN scale_parameters& p);

/// Функция для генерации шкалы частот из модели уха.
bool EXPORT scale_generate_model(OUT freq_scale& sc, IN scale_parameters& p);

///
/// То же, что и \ref freq_scale, но с управлением памяти (для C++).
///

struct EXPORT scale_class: freq_scale {

	/// Конструктор. Генерирует шкалу по заданной форме шкалы.
	scale_class(const scale_parameters& p, scale_form form);

	/// Конструктор. Генерирует шкалу по заданным функциям шкалы.
	scale_class(const scale_parameters& p, scale_func f, rscale_func r);

	/// Конструктор. Загружает шкалу из заданного файла.
	scale_class(const char *filename);

	/// Сохранение шкалы в файл. 
	/// Вызывает функцию save_scale().
	bool save(const char *file) { return save_scale(*this, file); }

	/// Деструктор. Освобождает память, выделенную под шкалу.
	~scale_class();
private:
	void init(int K);

};

/// Функция поиска частоты в шкале. 
/// Предполагается, что шкала возрастающая.
int EXPORT scale_find_freq(const freq_scale& sc, freq_t F);

/// Стандартное количество точек, по которому определяется тип шкалы.
const int DEFAULT_SCALE_ESTIMATE_ACCURACY = 5;

/// Функция проверки типа шкалы.
bool EXPORT is_scale_of_form(const freq_scale& sc, enum scale_form form, int npoints = DEFAULT_SCALE_ESTIMATE_ACCURACY);

/// Функция определения формы шкалы.
enum scale_form EXPORT scale_form_estimate(const freq_scale& sc, int npoints = DEFAULT_SCALE_ESTIMATE_ACCURACY);

/// Точная проверка типа шкалы.
/// Точность достигается за счет проверки всех точек шкалы.
inline bool is_scale_of_form_accurate(const freq_scale& sc, enum scale_form form) {
	return is_scale_of_form(sc, form, sc.K - 2);
}

/// Точное определение формы шкалы.
/// Точность достигается за счет проверки всех точек шкалы.
inline enum scale_form scale_form_estimate_accurate(const freq_scale& sc) {
	return scale_form_estimate(sc, sc.K - 2);
}

/// Интерполяция шкалы в другой точке.
freq_t EXPORT scale_interp_k2f(const freq_scale& sc, scale_form form, double k);

/// Интерполяция шкалы в другой частоте.
double EXPORT scale_interp_f2k(const freq_scale& sc, scale_form form, freq_t f);


} // namespace spl

#endif//_SPL_SCALE
