#ifndef _SPL_FILTERS_
#define _SPL_FILTERS_

#include "spl_types.h"

namespace spl {

///
/// Класс filters, который определяет интерфейс для классов фильтров
///

class filters {
public:

	/// Генерация фильтров
	virtual bool generate() = 0;

	/// Подготовка фильтров к быстрой фильтрации
	virtual bool prepare() = 0;

};


} // namespace spl

#endif//_SPL_FILTERS_