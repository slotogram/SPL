#ifndef _SPL_COMMON_
#define _SPL_COMMON_

///
/// \file  common.h
/// \brief ����� ����������� ���������� SPL.
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
/// ������������ ���� ������� � ������� SPL.
/// ������������ ��� ���� ������� � �������, ���������������
///  ��� ������� ���������� ���� 
///  � ���������� � ���������� SPL (Speech Parameterization Library).
///

namespace spl {

//@{
/// ������� ��� ���������� ������� � ���������� SPL.
void EXPORT *spl_alloc_low(size_t siz);
void EXPORT spl_free(void *addr);

template<typename T>
T *spl_alloc(size_t siz) {
	return (T*) spl_alloc_low(siz * sizeof(T));
}
//@}

}

#endif//_SPL_COMMON_
