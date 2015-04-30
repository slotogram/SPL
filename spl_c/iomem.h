#ifndef _IO_MEMORY_
#define _IO_MEMORY_

///
/// \file  iomem.h
/// \brief In-memory input and/or output streams.
///

#include "io.h"

namespace io {

////////////////////////////////////////////////////////////////////////
//                      MEMORY STREAM CLASSES                         //
////////////////////////////////////////////////////////////////////////

template<typename T> class abstract_mstream;
template<typename T> class  imstream;
template<typename T> class  omstream;
template<typename T> class iomstream;

/// 
/// Abstract memory stream.
/// Contains members independent of memory stream direction.
/// Should not be used in functions, use imstream, omstream and iomstream instead.
/// Should be used only as memory streams base class.
///

template<typename T>
class abstract_mstream: 
	 virtual public abstract_stream 
{
public:

	/// Memory stream size function.
	/// Returns memory buffer size.
	size_t size() const { return _size; }

	//@{
	/// Memory stream position.
	/// It is abstract_stream::pos function implementation.
	virtual size_t pos() const { return _pos; }
	virtual size_t pos(size_t newpos) { 
		if(newpos >=0 && newpos < _size) 
			_pos = newpos; 
		return _pos; 
	}
	//@}

	/// Checks whether end of stream is reached.
	/// It is abstract_stream::eos function implementation.
	virtual bool eos() const { return _closed || _pos == _size; }

	/// Close memory stream.
	virtual void close() { _closed = true; }

protected:

	/// Abstract memory stream constructor.
	/// Sets memory buffer and its size.
	abstract_mstream(T *data, size_t siz):
	  _size(siz), _pos(0), _data(data), _closed(false) {}

	/// Memory array size (in elements).
	size_t _size;

	/// Current position of stream.
	/// This element is mutable, because input streams are usually "const".
	mutable size_t _pos;

	/// Data pointer.
	/// Pointer is constant, data may be not.
	T * const _data;

	/// Closed flag
	bool _closed;

};

///
/// Input memory stream.
/// This class can be used for object creation.
/// Object imstream can be then passed as argument.
/// It is not recommended to require it in functions, 
///  it is better to require abstract input stream instead (istream).
///

template<typename T>
class imstream: 
	public abstract_mstream<T>, 
	public istream<T>
{
public:

	/// Input memory stream constructor.
	/// Sets memory buffer from arguments and position to 0.
	imstream(const T *data, size_t siz): 
		abstract_mstream<T>(const_cast<T *>(data), siz) {}

	/// Get next element.
	virtual bool get(T& x) {
		return eos() ? false : ((x = _data[_pos++]), true);
	}

};

///
/// Output memory stream.
/// This class can be used for object creation.
/// Object omstream can be then passed as argument.
/// It is not recommended to require it in functions, 
///  it is better to require abstract output stream instead (ostream).
///

template<typename T>
class omstream: 
	public abstract_mstream<T>, 
	public ostream<T>
{
public:

	/// Input memory stream constructor.
	/// Sets memory buffer from arguments and position to 0.
	omstream(T *data, size_t siz): 
		abstract_mstream<T>(data, siz) {}

	/// Put one element.
	virtual bool put(const T& x) {
		return eos() ? false : ((_data[_pos++] = x), true);
	}
	
};

} // namespace io 

#endif//_IO_MEMORY_
