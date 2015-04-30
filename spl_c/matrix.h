#ifndef _SPL_MATRIX_
#define _SPL_MATRIX_

///
/// \file  matrix.h
/// \brief Шаблонный класс матрицы.
/// 
/// Представленный класс фактически является указателем, удобным для работы с матрицами.
/// Размерность матрицы (2-мерная, 3-мерная, 4-мерная и т.д.) задана шаблонным параметром (жестко).
/// Сами размеры матрицы (2х3, 3х5х2 и т.д.) могут задаваться произвольными.
///

#include "common.h" // size_t

#include <stdarg.h>

/// Макрос для использования произвольного количества аргументов в функции.
/// Копирует аргументы в заранее заданный массив "name".
#define va2array(T, name, first, D) \
  va_list vl; \
  va_start(vl, first); \
  name[0] = first; \
  for(int i = 1; i < D; i++) { \
   name[i] = va_arg(vl, T); \
  } \
  va_end(vl);

/// Макрос для использования произвольного количества аргументов в функции.
/// Задает массив "name" и копирует в него аргументы.
#define va2array2(T, name, first, D) \
	T name[D]; \
	va2array(T, name, first, D)



/// Вспомогательная функция расчета линейного индекса матрицы.
size_t sub2ind(int D, const size_t dim[], const size_t sub[]);

/// Вспомогательная функция расчета линейного размера матрицы.
size_t dim2siz(int D, const size_t dim[]);


/// Класс многомерной матрицы.
/// 
/// Фактически является указателем на многомерный массив.
/// Для удобной работы с указателями на многомерные массивы.
/// 
/// Класс не делает какой-либо работы по выделению или освобождению памяти.

template<typename T, int D>
class Matrix {
public:
  /// Конструктор (принимает массив размерностей).
  /// Размерности должны приниматься в порядке уменьшения их значимости (размера адресуемой памяти).
  /// Так, последняя размерность должна адресовать ячейку матрицы, вторая - строки матрицы, и т.д.
  Matrix(T *data, const size_t dims[]): 
    _data(data), _size(dim2siz(D, dims))
  {
    memcpy(_dim, dims, sizeof(size_t) * D);
  }
  /// Конструктор (принимает размерности в списке аргументов).
  Matrix(T *data, size_t first, ...): 
    _data(data)
  {
  	va2array(size_t, _dim, first, D);
    _size = dim2siz(D, _dim);
  }
  
  //@{
  /// Многомерные селектор и модификатор.
  /// matrix(1,2) возвратит элемент, который располагается в первой строке, втором столбце.
  /// При этом матрица должна быть двумерной.
  /// 
  /// Внимание: размерность матрицы не проверяется! 
  /// (Ее невозможно проверить в такой абстрактной реализации класса матрицы).
  /// Если матрица является D-мерной, то в эту функцию должно быть передано ровно D аргументов типа size_t.
  /// Если будет передано больше аргументов, то в расчет будут приняты только первые D аргументов.
  /// Если будет передано меньше аргументов, то вероятно возникнет ошибка обращения к памяти.
  const T& operator()(size_t first, ...) const { 
    va2array2(size_t, sub, first, D);
    return _data[sub2ind(sub)]; 
  }

  T& operator()(size_t first, ...) { 
    va2array2(size_t, sub, first, D);
    return _data[sub2ind(sub)]; 
  }
  //@}

  //@{
  /// Одномерные селектор и модификатор
  const T& operator[](size_t i) const { return _data[i]; }
  T& operator[](size_t i) { return _data[i]; }
  //@}

  //@{
  /// Преобразование к типу T *
  operator const T *() const { return _data; }
  operator T *() { return _data; }
  //@}

  /// Возвращает массив размерностей матрицы.
  const size_t *dims() const { return _dim; }
  /// Возвращает заданную размерность матрицы.
  size_t dim(int i) const { return _dim[i]; }
  /// Возвращает количество элементов матрицы.
  size_t size() const { return _size; }
  /// Возвращает указатель на начало данных матрицы.
  T *ptr() { return _data; }

private:
  /// Измерения. 
  /// Массив с размерами матрицы. Перечисляет размеры матрицы от низшего к высшему
  size_t _dim[D];

  /// Размер данных.
  size_t _size;

  /// Указатель на данные.
  /// Собственно данные матрицы. Должны иметь размер _dim[0] * _dim[1] * ...
  T *_data;

  /// Расчет линейного индекса
  size_t sub2ind(size_t sub[]) const {
    switch(D) { // для оптимизации (возможно компилятор уберет ненужные ветки)
    case 1: return sub[0];
    case 2: return _dim[1] * sub[0] + sub[1];
    case 3: return _dim[2] * (_dim[1] * sub[0] + sub[1]) + sub[2];
    default: return ::sub2ind(D, _dim, sub);
    }
  }

};

/// Функция для быстрого получения двумерной матрицы.

template<typename T>
Matrix<T,2> matrix_ptr(T *data, size_t A, size_t B) {
  return Matrix<T,2>(data, A, B);
}

/// Функция для быстрого получения трехмерной матрицы.

template<typename T>
Matrix<T,3> matrix_ptr(T *data, size_t A, size_t B, size_t C) {
  return Matrix<T,3>(data, A, B, C);
}

/// Функция для быстрого получения четырехмерной матрицы.

template<typename T>
Matrix<T,4> matrix_ptr(T *data, size_t A, size_t B, size_t C, size_t D) {
  return Matrix<T,4>(data, A, B, C, D);
}

//@{
///
/// Изменение размера матрицы - reshape.
/// 

template<typename T, int D1, int D2>
Matrix<T,D2> reshape(Matrix<T,D1>& mtx1, const size_t dims2[]) {
  size_t siz2 = dim2siz(D2, dims2);
  if(mtx1.size() != siz2) 
    throw "Matrix sizes must agree on reshape";
  return Matrix<T,D2>(mtx1.ptr(), dims2);
}

template<typename T, int D1, int D2>
Matrix<T,D2> reshape(Matrix<T,D1>& mtx1, const size_t first, ...) {
  va2array2(size_t, dims, first, D2);
  return reshape<T,D1,D2>(mtx1, dims);
}
//@}

///
/// Транспонирование двумерной матрицы.
/// Для правильного транспонирования матрицы должны быть размеров MxN и NxM, 
///  или по крайней мере в выходной матрице должно быть место.
/// Матрицы также могут быть разных типов, 
///  главное, чтобы существовал оператор "=" для преобразования между ними.
///

template<typename T1, typename T2>
bool transpose(const Matrix<T1,2>& m, Matrix<T2,2>& t) {
  // проверяем, есть ли место в выходной матрице.
  if(m.dim(0) > t.dim(1)
  || m.dim(1) > t.dim(0))
    return false;
  // собственно транспонирование
  for(size_t i = 0; i < m.dim(0); i++) {
  for(size_t j = 0; j < m.dim(1); j++) {
    t(j,i) = m(i,j);
  }}
  return true;
}

#endif//_SPL_MATRIX_
