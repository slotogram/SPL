#ifndef _IO_SPLIT_
#define _IO_SPLIT_

#include "io.h"

namespace io {

template<typename T>
class tee:
	public ostream<T>
{
public:
	tee(ostream<T>& _first, ostream<T>& _second):
	  first(_first), second(_second)
	{}

	virtual size_t write(const T *block, size_t count) {
		first.write(block, count);
		second.write(block, count);
		return count;
	}

	size_t pos() const { return 0; }
	bool eos() const { return false; }
	void close() {}

private:
	ostream<T>& first, &second;
};


} // namespace io

#endif//_IO_SPLIT_