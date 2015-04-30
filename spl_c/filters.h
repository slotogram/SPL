#ifndef _SPL_FILTERS_
#define _SPL_FILTERS_

#include "spl_types.h"

namespace spl {

///
/// ����� filters, ������� ���������� ��������� ��� ������� ��������
///

class filters {
public:

	/// ��������� ��������
	virtual bool generate() = 0;

	/// ���������� �������� � ������� ����������
	virtual bool prepare() = 0;

};


} // namespace spl

#endif//_SPL_FILTERS_