#ifndef _SPL_COMMON_
#define _SPL_COMMON_

///
/// \file  common.h
/// \brief Общие определения библиотеки SPL.
///

#undef IN
#undef OUT
#undef INOUT

#define IN const
#define OUT
#define INOUT

#define EXPORT __declspec(dllexport) 

#define _USE_MATH_DEFINES

#include <stdio.h>

/// 
/// Пространство имен функций и классов SPL.
/// Используется для всех функций и классов, предназначенных
///  для расчета параметров речи 
///  и включенных в библиотеку SPL (Speech Parameterization Library).
///

namespace spl {

//@{
/// Функции для управления памятью в библиотеке SPL.
void EXPORT *spl_alloc_low(size_t siz);
void EXPORT spl_free(void *addr);

template<typename T>
T *spl_alloc(size_t siz) {
	return (T*) spl_alloc_low(siz * sizeof(T));
}
//@}

}

#endif//_SPL_COMMON_
