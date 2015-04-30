#ifndef _IO_BUFFER_
#define _IO_BUFFER_

#include <memory>
using std::auto_ptr;

template<typename T>
struct membuf:
	public auto_ptr<T>
{
public:

	membuf():
	  auto_ptr() {}

	membuf(size_t siz):
	  auto_ptr(new T[siz]) {}

	membuf(T *header):
	  auto_ptr(header) {}

	membuf(const membuf& x):
	  auto_ptr(const_cast<membuf&>(x)) {}

	operator T *() {
		return get();
	}

};

#endif//_IO_BUFFER_