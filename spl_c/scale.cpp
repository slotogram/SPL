#include "scale.h"
#include "iofile.h"
#include "const.h"
#include <math.h>

namespace spl {

///
/// Функция для сохранения шкалы в файл
///

bool save_scale(
	IN freq_scale& sc, 
	IN char *file) 
{
	return io::array_to_file(sc.Fr, sc.K, file);
}


// 
// Функции класса scale_class.
//

void scale_class::init(int K) {
	Fr = spl_alloc<freq_t>(K);
	if(!Fr) 
		throw "Can't allocate place for scale";
	this->K = K;
}

scale_class::scale_class(const scale_parameters& p, scale_form form) {
	init(p.b - p.a + 1);
	if(!scale_generate(*this, p, form)) 
		throw "Can't generate scale";
}

scale_class::scale_class(const scale_parameters& p, scale_func f, rscale_func r)
{
	init(p.b - p.a + 1);
	if(!scale_generate(*this, p, f, r)) 
		throw "Can't generate scale";
}

scale_class::scale_class(const char *filename) {
	io::ifstream<freq_t> in(filename);
	size_t K = in.size();
	this->K = (int)K;
	if(K <= 0) throw "Can't load scale";
	Fr = spl_alloc<freq_t>(K);
	if(!Fr) 
		throw "Can't allocate place for scale";
	if(in.read(Fr, K) != K)
		throw "Can't read scale from file";
}
	
scale_class::~scale_class() { spl_free(Fr); }

//
// Различные функции для генерации шкал.
// 

// Универсальная функция генерации шкалы.
// Фактически любая задача генерации шкалы проходят через эту функцию
// (в том числе интерполяция и экстраполяция шкал).
bool scale_generate(
	OUT freq_scale& sc,
	IN scale_parameters& p,
	IN scale_func f2x,
	IN rscale_func x2f)
{
	if(sc.K != p.b - p.a + 1) return false;

	// расчет шага
	double Xi = f2x(p.Fi);
	double Xj = f2x(p.Fj);
	double dx = (Xi - Xj) / (p.i - p.j);

	// собственно вычисление частот	
	for(int n = p.a, k = 0; n <= p.b; n++, k++) {
		double Xn = Xi + dx * (n - p.i);
		sc.Fr[k] = x2f(Xn);
	}
	return true;
}

// Функции стандартных шкал

//@{
/// Прямая и обратная функция линейной шкалы частот.
double lin_scale(freq_t x) { return (double) x; }
freq_t rlin_scale(double x) { return (freq_t) x; }
//@}

//@{
/// Прямая и обратная функция логарифмической шкалы частот.
double log_scale(freq_t x) { return (double) log(x); }
freq_t rlog_scale(double x) { return (freq_t) exp(x); }
//@}

//@{
/// Прямая и обратная функция мел-шкалы частот.
/// Используется апроксимация B(f) = 1125 * ln(1 + f/700).
double mel(freq_t f) { return 1125 * log(1 + double(f)/700); }
freq_t rmel(double B) { return freq_t(700 * (exp(B/1125) - 1)); }
//@}

//@{
/// Прямая и обратная функция барк-шкалы частот.
/// Используется апроксимация B(f) = 7* arcsinh(f/650)
/// и arcsinh(x) = ln(x + sqrt(x^2 + 1))
double bark(freq_t f) { return 7 * log(f/650 + sqrt(1+f*f/422500)); }
freq_t rbark(double B) { return freq_t(650 * sinh(B/7)); }
//@}


namespace {

struct _scale_funcs {
 scale_func f2x;
 rscale_func x2f;
}
scale_funcs[] = {
 { 0, 0 },
 { lin_scale,  rlin_scale  },
 { log_scale,  rlog_scale  },
 { mel,        rmel        },
 { bark,       rbark       },
 { model_f2x,  model_x2f   },
 { 0, 0 }
};

} // namespace {

// Функция, генерирующая шкалу по ее форме
bool scale_generate(
	OUT freq_scale& sc,     ///< Шкала, которая будет сгенерирована
	IN scale_parameters& p, ///< Параметры генерации шкалы
	IN scale_form form)     ///< Форма шкалы
{
	return scale_generate(sc, p, scale_funcs[form].f2x, scale_funcs[form].x2f);
}

//@ 
/// Функции для генерации стандартных шкал.

bool scale_generate_linear(OUT freq_scale& sc, IN scale_parameters& p) { return scale_generate(sc, p, scale_linear); }
bool scale_generate_model (OUT freq_scale& sc, IN scale_parameters& p) { return scale_generate(sc, p, scale_model ); }
bool scale_generate_log   (OUT freq_scale& sc, IN scale_parameters& p) { return scale_generate(sc, p, scale_log   ); }
bool scale_generate_mel   (OUT freq_scale& sc, IN scale_parameters& p) { return scale_generate(sc, p, scale_mel   ); }
bool scale_generate_bark  (OUT freq_scale& sc, IN scale_parameters& p) { return scale_generate(sc, p, scale_bark  ); }
//@}

namespace {

// Линейный поиск частоты в шкале.
int scale_find_freq_linear(const freq_scale& sc, freq_t F) {
	for(int k = 1; k < sc.K; k++) {
		if(sc.Fr[k] >= F) {
			return (sc.Fr[k] - F < F - sc.Fr[k - 1]) ? k : k - 1;
		}
	}
	return sc.K - 1;
}

// Двоичный поиск частоты в шкале.
int scale_find_freq_binary(const freq_scale& sc, freq_t F) {
	int a = 0, b = sc.K-1;
	if(F < sc.Fr[a]) return a;
	if(F > sc.Fr[b]) return b;
	while(b - a > 1) {
		int c = (b + a) / 2;
		if(F < sc.Fr[c]) b = c;
		else if(F > sc.Fr[c]) a = c;
		else return c;
	}
	return (F - sc.Fr[a] < sc.Fr[b] - F) ? a : b;
}

} // namespace {

/// Функция поиска частоты в шкале. 
int scale_find_freq(const freq_scale& sc, freq_t F) {
	return scale_find_freq_binary(sc, F);
}

/// Максимальная ошибка функции scale_form_estimate
const double SCALE_ESTIMATE_ERROR = 1E-10;

///
/// Функция определения формы шкалы.
/// 
/// Пытается определить форму заданной шкалы частот \a sc.
/// Для этого для всех известных типов шкал вызывается функция \ref is_scale_of_form.
/// Если не подошло ни одной шкалы, то возвращается значение scale_unknown (см. \ref scale_form).
///
/// Для точного определения формы шкалы должно выполняться равенство \a npoints = \a sc.K - 2. 
///  То есть должны проверяться все внутренние точки.
///
enum scale_form scale_form_estimate(const freq_scale& sc, int npoints) {

	// перебираем шкалы
	for(int i = 1; scale_funcs[i].f2x != NULL; i++) {
		enum scale_form form = static_cast<scale_form>(i);
		if(is_scale_of_form(sc, form, npoints)) {
			return form;
		}
	}

	// ни одна из известных шкал не подошла
	return scale_unknown;
}

///
/// Проверка типа шкалы.
///
/// Функция проверяет, совпадает ли форма заданной шкалы с заданной формой.
/// Для этого производится расчет новой шкалы в нескольких точках 
///  (количество точек задается аргументом \a npoints).
/// Если суммарное отличие частот в новой шкале от частот в старой шкале
///  не превышает константы \ref SCALE_ESTIMATE_ERROR, то шкала считается подходящей.
///
/// Для точной проверки формы шкалы должно выполняться равенство \a npoints = \a sc.K - 2. 
///  То есть должны проверяться все внутренние точки.
///
bool is_scale_of_form(const freq_scale& sc, enum scale_form form, int npoints) {
	const freq_t *Fr = sc.Fr;
	int K = sc.K;

	// функции шкалы:
	scale_func f2x = scale_funcs[form].f2x;
	rscale_func x2f = scale_funcs[form].x2f;

	// расчет шага:
	double x1 = f2x(Fr[0]);
	double x2 = f2x(Fr[K-1]);
	double dx = (x2 - x1) / (K - 1);

	// вычисление ошибки интерполяции
	double e = 0.0;
	for(int j = 1; j < npoints; j++) {
		int k = j * K / (npoints + 1);
		e += fabs(Fr[k] - x2f(x1 + k * dx));
	}

	// если суммарная ошибка по точкам меньше порога - то мы обнаружили шкалу
	return e < SCALE_ESTIMATE_ERROR;
}

//
// Функции интерполяции шкалы в других точках (частотах)
//

/// Интерполяция шкалы в другой точке.
freq_t scale_interp_k2f(const freq_scale& sc, scale_form form, double k) {
	int K = sc.K;
	double x1 = scale_funcs[form].f2x(sc.Fr[0]);
	double x2 = scale_funcs[form].f2x(sc.Fr[K-1]);
	double dx = (x2 - x1) / (K - 1);

	return scale_funcs[form].x2f(x1 + k * dx);
}

/// Интерполяция шкалы в другой частоте.
double scale_interp_f2k(const freq_scale& sc, scale_form form, freq_t f) {
	int K = sc.K;
	double x1 = scale_funcs[form].f2x(sc.Fr[0]);
	double x2 = scale_funcs[form].f2x(sc.Fr[K-1]);
	double dx = (x2 - x1) / (K - 1);
	double x = scale_funcs[form].f2x(f);
	
	return (x - x1) / dx;
}


} // namespace spl
