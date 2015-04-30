#include "matrix.h"

///
/// \file  matrix.cpp
/// \brief Реализация функций, обслуживающих класс матрицы, 
///        не зависящих от шаблонных параметров.
///

/// Расчет линейного индекса.
/// Линейный индекс рассчитывается по формуле: 
/// i = dim[N] * ( dim[N-1] * ( ... * dim[3] ( dim[2] * sub[1] + sub[2] ) + sub[3] ) + ... ) + sub[N-1]) + sub[N]

size_t sub2ind(int D, const size_t dim[], const size_t sub[]) {
	size_t index = 0;
	for(int d = 0; d < D; d++) {
// забыто для скорости
//		if(sub[d] >= dim[d] || sub[d] < 0) 
//			throw "Matrix dimension exceeded";
		index = index * dim[d] + sub[d];
	}
	return index;
}

/// Расчет размера многомерной матрицы.
/// Размер расчитывается по формуле: dim[1] * dim[2] * ... * dim[N]
/// Если какая-либо размерность неположительная - выбрасывает исключение.

size_t dim2siz(int D, const size_t dim[]) {
	size_t siz = 1;
	for(int d = 0; d < D; d++) {
		if(dim[d] <= 0) 
			throw "Non-positive dimension received";
		siz *= dim[d];
	}
	return siz;
}
