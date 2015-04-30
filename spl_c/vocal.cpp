///
/// \file  vocal.cpp
/// \brief Функции для обработки звуков по признаку вокализованности.
///

#include "vocal.h"
#include "scale.h"
#include "const.h"
#include "mask.h"
#include "matrix.h"
#include "iomem.h"
#include "iofile.h"

#define _TEMPLATE_SEARCH_FAST_

//#define __MPN(x) ___gmpn_##x
#include <mpir.h>
#include <algorithm>
#include <queue>
using std::fill_n;

namespace spl {

/// Генерация шаблона маски для определения частоты основного тона.
/// Память в аргументе \a tpl.T выделяется функцией.
bool pitch_templates_generate(
	IN freq_scale& sc, 
	OUT mask_templates& tpl,
	IN mask_parameters& pm, 
	IN pitch_parameters& p
) {

	bool r;

	// 1. вычисление k1, k2
	int k1 = scale_find_freq(sc, p.F1);
	int k2 = scale_find_freq(sc, p.F2);
	int nk = k2 - k1 + 1;

	// 2. Формирование тестовых сигналов

	freq_t Fs = 12000; // частота тестового сигнала в принципе не имеет значения
	double Sm = floor(double(Fs) / p.F1); // длительность тестовых сигналов
	int s = (int) floor(Sm / 2); // точка, в которой рассчитывается интенсивность

	// выделяем место: 
	spectrum_t *I1 = spl_alloc<spectrum_t>(nk * sc.K);
	Matrix<spectrum_t, 2> I = matrix_ptr(I1, nk, sc.K);

	// собственно вычисления:
	for(int kt = k1; kt <= k2; kt++) {
		for(int k = 0; k < sc.K; k++) {
			double yc = 0.0, ys = 0.0;
			for(int n = 1; n <= p.Nh; n++) {
				double tmp = f_C24 / 2 * model_filter_quality(sc.Fr[k]) * (1 - n * sc.Fr[kt]/sc.Fr[k]);
				double H = exp( - tmp * tmp );
				yc += H * cos(2 * M_PI * n * sc.Fr[kt] / Fs * s);
				ys += H * sin(2 * M_PI * n * sc.Fr[kt] / Fs * s);
			}
			I(kt - k1, k) = yc * yc + ys * ys;
		}
	}

//	io::array_to_file<spectrum_t>(I1, nk * sc.K, "pitch-test-spec.bin");

	// 3. Одновременная маскировка
	
	// выделяем место:
	mask_t *M1 = spl_alloc<mask_t>(nk * sc.K);
	Matrix<mask_t, 2> M = matrix_ptr(M1, nk, sc.K);

	spectrum sp;
	sp.F = Fs; sp.K = sc.K; sp.N = nk; sp.Y = I1;
	mask m;
	m.K = sc.K; m.N = nk; m.Z = M1;

	if(mask_memory(sc, sp, m, pm) != nk * sc.K) { 
		// обработка ошибки
		r = false;
		goto out;
	}

//	io::array_to_file<mask_t>(M1, nk * sc.K, "pitch-test-mask.bin");

	// 4. Вычисление шаблона

	// выделение памяти:
	tpl.k1 = k1;
	tpl.k2 = k2;
	tpl.K = sc.K;
	tpl.Nt = 2 * p.Nh + 2;
	tpl.T = spl_alloc<int>(nk * tpl.Nt);

	// границы единичных областей маскировки (0 -> 1, 1 -> 0 - для каждой учитываемой гармоники)
	int *borders = tpl.T; 

	// для каждого канала ЧОТ:
	for(int kt = 0; kt < nk; kt++) {
		int k;
		
		// вначале смотрим на первую гармонику (ЧОТ)
		// там должна быть единица
		if(M(kt, k1 + kt) != 1) {
			spl_free(tpl.T);
			r = false;
			goto out;
		}

		// ищем первую границу слева:
		for(k = k1 + kt - 1; M(kt, k) == 1; k--);
		borders[1] = k + 1;

		// определяем границы единичных областей маскировки
		int n = 2;
		for(int k = k1 + kt + 1; k < sc.K; k++) {
			if(M(kt, k) == M(kt, k-1)) 
				continue;
			borders[n++] = k;
			if(n > 2 * p.Nh) 
				break;
		}

		// добавляем границы:
		borders[0] = borders[1] - (borders[3] - borders[2] + 2) / 3;
		borders[n] = borders[n-1] + (borders[n-2] - borders[n-3] + 2) / 3;
		borders[0] = borders[0] < 0 ? 0 : borders[0];
		borders[n] = borders[n] >= sc.K ? sc.K - 1 : borders[n];

		borders += tpl.Nt;

	}

	r = true;

out:

	spl_free(I1);
	spl_free(M1);

	return r;
}


/// Поиск в потоке маски по шаблону. 
/// Наиболее общая функция поиска в маске по шаблону - использует абстрактные потоки.
/// На выход пишутся номера наиболее подходящих шаблонов или -1, если превышен порог \a max_diff.
/// Возвращает количество записанных на выход элементов.
size_t mask_template_search(
	IN mask_templates& tpl,       ///< шаблоны маски, по которым осуществляется поиск
	io::istream<mask_t>& in_str,  ///< входной поток (маска)
	io::ostream<short>& out_str,  ///< выходной поток (номера найденных каналов)
	IN int max_diff               ///< максимальное отличие маски от шаблона
) {
	int K = tpl.K, 
		k1 = tpl.k1, 
		k2 = tpl.k2, 
		nk = k2 - k1 + 1;
	int Nt = tpl.Nt;

	// матрица шаблонов
	Matrix<int, 2> T = matrix_ptr(tpl.T, nk, tpl.Nt);

	// входной буфер
	mask_t *M = spl_alloc<mask_t>(K);

	// количество записанных элементов - возвращаемое значение
	size_t N = 0;

	// основной цикл - проходимся по отсчетам
	while(in_str.read(M, K) == K) {

		int dmin = max_diff; // минимальная найденная разница
		int kmin = -1; // номер канала, где была найдена минимальная разница
		// если меньше, чем max_diff, не найдется, то будет выведено -1

		// цикл выбора оптимального шаблона:
		for(int i = 0; i < nk; i++) { // цикл по шаблонам
			int di = 0; // количество расхождений с шаблоном
			mask_t mj = 0; // правильное значение маски на интервале - начинается с нуля
			for(int j = 1; j < Nt; j++, mj ^= 1) { // цикл по элементам шаблона (номера каналов)
				for(int k = T(i,j-1); k < T(i,j); k++) { // цикл по каналам
					// если значение маски не совпадает с правильным - увеличивается разница
					di += (M[k] != mj);
				}
				// если хуже, чем уже было найдено - следующий шаблон
				if(di >= dmin) break;
			}
			// если нашли лучше - запоминаем
			if(di < dmin) {
				dmin = di;
				kmin = k1 + i;
			}
		}

		// выводим номер канала:
		N += out_str.put(kmin);
	}

	return N;
}

struct limb_template {
	mp_limb_t mask;
	mp_limb_t value;
};

mp_limb_t * shift1(mp_limb_t *s1p, mp_size_t s1n, unsigned bit_num) {

	fill_n(s1p, s1n, 0);

	const unsigned limb_bits = sizeof(mp_limb_t) * 8;
	unsigned index = bit_num / limb_bits;
	unsigned bit = bit_num % limb_bits;

	s1p[index] = s1p[index] | (mp_limb_t(1) << bit);
	return s1p;
}

inline 
mp_limb_t get_part(mp_limb_t *s1p, mp_size_t s1n, int num_part_bits, int num_part) {

	if(num_part_bits == 16) {
		unsigned short *x = (unsigned short *) s1p;
		return x[num_part];
	}

	if(num_part_bits == 8) {
		unsigned char *x = (unsigned char *) s1p;
		return x[num_part];
	}

	if(num_part_bits == 4) {
		unsigned char *x = (unsigned char *) s1p;
		unsigned char two_parts = x[num_part / 2];
		if(num_part % 1 == 0) {
			return two_parts & 0x0F;
		} else {
			return two_parts >> 4;
		}
	}

	throw "such num_part_bits value is not supported";
}

#define CEIL_MODULUS(x,y) ( ( (x) + (y) - 1 ) / (y) )

size_t mask_template_fast_search_tpl(
	IN mask_templates& tpl,       ///< шаблоны маски, по которым осуществляется поиск
	io::istream<mask_t>& in_str,  ///< входной поток (маска)
	io::ostream<short>& out_str,  ///< выходной поток (номера найденных каналов)
	IN int num_part_bits,         ///< количество бит в куске
	IN int max_diff               ///< максимальное отличие маски от шаблона
) {
	// тип для хранения количества отличий (diff count) между сэмплом и маской
	// для этого достаточно одного байта
	typedef unsigned char diff_t;

	// параметры шаблонов
	const int k1 = tpl.k1, // канал первого шаблона
		k2 = tpl.k2,       // канал последнего шаблона
		num_templates = k2 - k1 + 1;  // количество шаблонов

	// количество каналов = количество бит в сэмпле
	const int K = tpl.K;
	const int num_sample_bits = tpl.K;

	// количество (неполных) байт в сэмпле
	const int num_sample_bytes = CEIL_MODULUS(num_sample_bits, 8);

	// количество (неполных) чисел в сэмпле
	const int num_sample_limbs = CEIL_MODULUS(num_sample_bytes, sizeof(mp_limb_t));

	// количество кусков в сэмпле
	const int num_sample_parts = CEIL_MODULUS(num_sample_bits, num_part_bits);

	// количество вариантов куска - 2^num_bits
	const int num_part_variants = 1 << num_part_bits;

	// количество чисел в строке с отличиями
	const int num_diff_limbs = CEIL_MODULUS(num_templates * sizeof(diff_t), sizeof(mp_limb_t));

	// количество отличий в одном числе (их должно быть ровное количество)
	const int num_limb_diffs = sizeof(mp_limb_t) / sizeof(diff_t);

	// diff count в строке (вместе с неиспользуемыми)
	const int num_diffs_all = num_diff_limbs * num_limb_diffs;

	// инициализация таблиц поиска
	// выделяем место под таблицы
	mp_limb_t *mem_tables = spl_alloc<mp_limb_t>(num_sample_parts * num_part_variants * num_diff_limbs
		+ 2 * num_templates * num_sample_limbs + num_sample_parts * 2);
	mp_limb_t *buffer = mem_tables;

	// место под шаблоны
	Matrix<mp_limb_t,2> tpl_masks = matrix_ptr(buffer, num_templates, num_sample_limbs);
	buffer += tpl_masks.size();
	Matrix<mp_limb_t,2> tpl_values = matrix_ptr(buffer, num_templates, num_sample_limbs);
	buffer += tpl_values.size();

	// ссылка для быстрого сложения
	Matrix<mp_limb_t,3> sum_tables = matrix_ptr(buffer, num_sample_parts, num_part_variants, num_diff_limbs);
	buffer += sum_tables.size();

	// ссылка для доступа к diff count
	Matrix<diff_t,3> diff_tables = matrix_ptr((diff_t *)sum_tables.ptr(), num_sample_parts, num_part_variants, num_diffs_all);

	// индексы ненулевых элементов для быстрого суммирования
	Matrix<mp_limb_t,2> sum_indices = matrix_ptr(buffer, num_sample_parts, 2);

	// входная матрица шаблонов
	Matrix<int, 2> tpl_indices = matrix_ptr(tpl.T, num_templates, tpl.Nt);

	// предварительно вычисляем шаблоны
	int Nt = tpl.Nt;
	for(int k = 0; k < num_templates; k++) {

		mp_limb_t *mask = &tpl_masks(k, 0);
		mp_limb_t *value = &tpl_values(k, 0);

		// mask = (1 << indices(k,end)) - (1 << indices(k,0))
		mpn_sub_n(mask, 
			shift1(mask, num_sample_limbs, tpl_indices(k, Nt-1)), 
			shift1(buffer, num_sample_limbs, tpl_indices(k, 0)), 
			num_sample_limbs);

		// value заполняется по аналогии
		fill_n(value, num_sample_limbs, 0);
		bool add = true;
		for(int i = Nt - 2; i >= 1; i--) {
			shift1(buffer, num_sample_limbs, tpl_indices(k, i));
			if(add) {
				mpn_add_n(value, value, buffer, num_sample_limbs);
			} else {
				mpn_sub_n(value, value, buffer, num_sample_limbs);
			}
			add = !add;
		}

	}

	// заполняем таблицы разностей
	// заполняем индексы суммирования
	for(int i = 0; i < num_sample_parts; i++) {
		mp_limb_t k1 = 0, k2 = 0;
		bool found = false;
		for(int k = 0; k < num_templates; k++) {
			mp_limb_t x = get_part(&tpl_values(k, 0), num_sample_limbs, num_part_bits, i);
			mp_limb_t y = get_part(&tpl_masks(k, 0), num_sample_limbs, num_part_bits, i);

			// таблица разностей:
			for(mp_limb_t j = 0; j < num_part_variants; j++) {
				// нам нужно расстояние хэмминга между i-той частью k-того шаблона и j-тым вариантом i-той части
				mp_limb_t z = j & y;
				diff_tables(i, j, k) = (diff_t) mpn_hamdist(&x, &z, 1);
			}

			// индексы суммирования
			if(!found) {
				if(y != 0) {
					k1 = k2 = k;
					found = true;
				}
			} else {
				if(y != 0) {
					k2 = k;
				}
			}
		}
		sum_indices(i, 0) = found ? k1 / num_limb_diffs : 0;
		sum_indices(i, 1) = found ? k2 / num_limb_diffs + 1: 0;
	}

#ifdef _DEBUG
	io::array_to_file<mp_limb_t>(tpl_masks, tpl_masks.size(), "pitch-templates-masks.bin");
	io::array_to_file<mp_limb_t>(tpl_values, tpl_values.size(), "pitch-templates-values.bin");
	io::array_to_file<diff_t>(diff_tables, diff_tables.size(), "pitch-diff-tables.bin");
#endif
	// процедура сравнения с шаблонами
	buffer = mem_tables;
	mp_limb_t *diff_limbs = buffer; buffer += num_diff_limbs;
	mp_limb_t *input_limbs = buffer; buffer += num_sample_limbs;
	mp_limb_t *input2_limbs = buffer; buffer += num_sample_limbs;
	mask_t *input = (mask_t *) buffer;
	diff_t *diffs = (diff_t *) diff_limbs;
	
	const int limb_bits = sizeof(mp_limb_t) * 8;
	bool first = true;
	while(in_str.read(input, K) == K) {

		// конвертируем в биты
		fill_n(input_limbs, num_sample_limbs, 0);
		for(size_t i = 0; i < K; i++) {
			size_t n_byte = i / limb_bits;
			size_t n_bit  = i % limb_bits;
			if(input[i]) {
				input_limbs[n_byte] |= mp_limb_t(1) << n_bit;
			}
		}

		// считаем разницы для всех шаблонов
		if(first) {
			// если все происходит в первый раз,
			// нужно посчитать разницу полностью
			fill_n(diff_limbs, num_diff_limbs, 0);
			for(int i = 0; i < num_sample_parts; i++) {
				mp_limb_t part_value = get_part(input_limbs, num_sample_limbs, num_part_bits, i);
				mp_limb_t *tpl_diffs = &sum_tables(i, part_value, 0);
				mp_limb_t k1 = sum_indices(i, 0);
				mp_limb_t k2 = sum_indices(i, 1);
				if(k2 > k1) {
					mpn_add_n(diff_limbs + k1, diff_limbs + k1, tpl_diffs + k1, k2 - k1);
				}
			}

		} else {
			// если все происходит не в первый раз,
			// нужно обновлять только если часть изменилась
			for(int i = 0; i < num_sample_parts; i++) {
				mp_limb_t current_value = get_part(input_limbs, num_sample_limbs, num_part_bits, i);
				mp_limb_t previous_value = get_part(input2_limbs, num_sample_limbs, num_part_bits, i);
				mp_limb_t k1 = sum_indices(i, 0);
				mp_limb_t k2 = sum_indices(i, 1);
				if(current_value != previous_value && k2 > k1) {
					// вычитаем предыдущие значения
					{
						const mp_limb_t part_value = previous_value;
						mp_limb_t *tpl_diffs = &sum_tables(i, part_value, 0);
						mpn_sub_n(diff_limbs + k1, diff_limbs + k1, tpl_diffs + k1, k2 - k1);
					}
					// добавляем текущие значения
					{
						const mp_limb_t part_value = current_value;
						mp_limb_t *tpl_diffs = &sum_tables(i, part_value, 0);
						mpn_add_n(diff_limbs + k1, diff_limbs + k1, tpl_diffs + k1, k2 - k1);
					}
				}
			}
		}

		// ищем наименьшее отличие от шаблона
		diff_t min_index = -1, min_diff = -1;
		for(int k = 0; k < num_templates; k++) {
			if(diffs[k] < min_diff) {
				min_index = k;
				min_diff = diffs[k];
			}
		}

		// если нас устраивает это отличие, мы выводим номер канала, иначе - 0
		if(min_diff < max_diff) {
			int k = k1 + min_index;
			out_str.put(k);
		} else {
			out_str.put(-1);
		}

		// save old input
		mp_limb_t *tmp = input_limbs;
		input_limbs = input2_limbs;
		input2_limbs = tmp;
		first = false;
	
	}

out:
	spl_free(mem_tables);

	return 0;
}

size_t mask_template_fast_search(
	IN mask_templates& tpl,       ///< шаблоны маски, по которым осуществляется поиск
	io::istream<mask_t>& in_str,  ///< входной поток (маска)
	io::ostream<short>& out_str,  ///< выходной поток (номера найденных каналов)
	IN int max_diff               ///< максимальное отличие маски от шаблона
) {
	return mask_template_fast_search_tpl(tpl, in_str, out_str, 8, max_diff);
}


// Функция определения частоты основного тона
size_t pitch_track(
	IN freq_scale& sc, 
	io::istream<mask_t>& in_str,
	io::ostream<short>& out_str,
	IN mask_parameters& pm, 
	IN pitch_parameters& p
) {
	try {
		// генерация шаблонов для поиска ЧОТ:
		pitch_templates_class tpl(sc, pm, p);

		// вызов функции поиска в маске по шаблону:
#ifdef _TEMPLATE_SEARCH_FAST_
		return mask_template_fast_search(tpl, in_str, out_str, DEFAULT_PITCH_MAX_DIFF);
#else
		return mask_template_search(tpl, in_str, out_str, DEFAULT_PITCH_MAX_DIFF);
#endif
	} catch(...) {
		return 0;
	}
}

inline void keep(std::queue<short>& q, io::ostream<short>& o, short c) { q.push(c); }
inline void reset(std::queue<short>& q, io::ostream<short>& o, short c) {
	while(!q.empty()) q.pop(); 
}
inline void flush(std::queue<short>& q, io::ostream<short>& o, short c) {
	while(!q.empty()) {
		o.put(q.front());
		q.pop();
	}
	o.put(c);
}

// функция сегментации сигнала по признаку вокализованности
size_t vocal_segment(
	io::istream<short>& in_str,
	io::ostream<short>& out_str,
	IN freq_t F,
	IN vocal_parameters& p
) {
	// минимальная длительность вокализованного участка в отсчетах:
	int minV = int(p.minV * F);
	// минимальная длительность невокализованного участка в отсчетах:
	int minNV = int(p.minNV * F);

	bool VocP = false, Voc;
	short kffp = 0, kff;
	int Tnv = 0, Tkv = 0, T = 0;

	size_t written = 0;

	std::queue<short> q;

	enum { out_v, out_nv, keep } action = keep;

	while(in_str.get(kff)) {
		kff = kff + 1; // приводим к шкале 1:K
		Voc = kff != 0;

		// посреди вокализованного сегмента
		if(VocP && Voc) { 
			// должен быть плавный переход между каналами
			if( kff - kffp < 2 && kffp - kff < 2 ) {
				// достигли минимальной длительности вокализованного
				if(T - Tnv == minV 
				// и невокализованный сегмент достаточно большой длительности
				&& Tnv - Tkv > minNV
				){
					// фиксируем границу Tnv

					// нужно сбросить все, что запомнили
					action = out_v;
				}
			// если переход неплавный
			} else {
				// то звук невокализованный
				Voc = false;
			}
		}

		// был невокализованный, стал вокализованный
		if(!VocP && Voc) {
			Tnv = T;
			action = keep;
		}

		// был вокализованный, стал невокализованный
		if(VocP && !Voc) {
			// вокализованный сегмент достаточной длительности
			if(T - Tnv > minV) {
				Tkv = T;
				action = keep;
			// невокализованный сегмент достаточной длительности
			} else if(T - Tkv > minNV) {
				// невокализованный сегмент недостаточной длительности
				if(Tnv - Tkv <= minNV) {
					// фиксируем границу - Tkv
					action = out_nv;
				}
			}
		}

		// посреди невокализованного участка
		if(!VocP && !Voc) {
			// достигли минимальной длительности невокализованного участка
			if(T - Tkv == minNV) {
				// фиксируем границу - Tkv
				action = out_nv;
			}
		}

		switch(action) {
		case out_v: 
			while(!q.empty()) {
				if(out_str.put(q.front() - 1))
					written++;
				q.pop();
			}
			if(out_str.put(kff - 1)) 
				written++;
			break;
		case out_nv:
			while(!q.empty()) {
				if(out_str.put(-1))
					written++;
				q.pop();
			}
			if(out_str.put(-1)) 
				written++;
			break;
		case keep:
			q.push(kff);
			break;
		}

		VocP = Voc;
		kffp = kff;
		T++;
	}

	//дополнил Рахманенко. Очищаем очередь, если у нас фикс. файл, а не поток.
	while(!q.empty()) {
				if(out_str.put(-1))
					written++;
				q.pop();
			}
			
	return written;
}

} // namespace spl
