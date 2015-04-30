#pragma warning (disable: 4005) // macro (IN, OUT) redefinition

#include <windows.h>
#include "common.h"

namespace spl {

void *spl_alloc_low(size_t siz) {
	return new char[siz];
}

void spl_free(void *addr) {
	delete [] addr;
}

} // namespace spl

