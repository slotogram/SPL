#ifndef _IO_BASE_
#define _IO_BASE_

///
/// \file  io.h
/// \brief Abstract templates of input and/or output streams.
/// 
/// These should be used as parent classes of "real" streams (in-memory, files, etc).
/// In functions, that need streams, istream and ostream should be used as function arguments.
///

#include <memory.h>
#include "common.h"

// disable warnings about iomstream multiple inheritance (inherits <xxx> via dominance)...
#pragma warning (disable: 4250) 


///
/// Alternative streams template library namespace.
/// Contains abstract classes of input and output streams
///  and their implementation in different media (memory, files, ...).
/// 

namespace io {

///
/// Stream mode.
///

enum stream_mode {
	in   = 0x01, ///< input mode - stream can input data
	out  = 0x02, ///< output mode - stream can output data
};


////////////////////////////////////////////////////////////////////////
//                     ABSTRACT STREAM CLASSES                        //
////////////////////////////////////////////////////////////////////////

class abstract_stream;
template<typename T> class  istream;
template<typename T> class  ostream;

///
/// Abstract class for all streams.
/// Contains members independent of stream element data type.
/// Should not be used in functions, use \ref istream, \ref ostream and \ref iostream instead.
/// Should be used only as base class for other streams.
///

class EXPORT abstract_stream 
{
public: 

	//@{
	/// Stream pos (position) functions.
	/// Return input stream position.
	/// Second function also sets stream position.
	/// Note: input/output streams have both input and output positions, so they may have two pos'es as well.
	virtual size_t pos() const = 0;
	virtual size_t pos(size_t newpos) { return skip(newpos - pos()); }
	//@}

	/// Stream skip function.
	/// Skips \a N elements from stream.
	/// Returns new stream position.
	/// \a N can be negative also.
	virtual size_t skip(size_t N) { return pos(pos() + N); }

	/// Stream eos (end of stream) function.
	/// Returns true if input stream reached its end.
	virtual bool eos() const = 0;

	/// Close the stream. Free all resources.
	virtual void close() = 0;

};

///
/// Abstract input stream.
/// Should be used in functions to show *any* input stream.
/// You can not construct istream object, 
/// you should construct its child classes instead, i.e. imstream.
///

template<typename T>
class istream: 
	virtual public abstract_stream
{
public:

	/// Get one element from stream.
	/// Returns true if get was successful, 
	///  otherwise returns false and \a x is not modified.
	/// Being pure virtual this function should be defined in child classes.
	virtual bool get(T& x) {
		return read(&x, 1) == 1;
	}

	//@{
	/// Main read functions.
	/// Accept buffer to read into and its size. 
	/// Return number of elements read.
	/// Usually it is \a count elements, 
	///  otherwise end of input stream (eos) has been reached.
	///
	/// Has two implementations: pure virtual and non-virtual.
	/// Pure virtual implementation should be defined in child classes.
	/// Non-virtual implementation is based on \ref get function 
	///  and need not to be defined in child classes.
	template<typename T1>
	size_t read(T1 *buf, size_t count) {
		T x;
		size_t j;
		for(j = 0; j < count; j++) {
			if(!get(x)) break;
			buf[j] = x;
		}
		return j;
	}

	virtual size_t read(T *buf, size_t count) {
		return read<T>(buf, count);
	}
	//@}

};


///
/// Abstract output stream.
/// Should be used in functions to show *any* output stream.
/// You can not construct ostream object, 
/// you should construct its child classes instead, i.e. omstream.
///

template<typename T>
class ostream: 
	virtual public abstract_stream
{
public:

	/// Put one element into stream.
	/// Returns true if put was successful, 
	///  otherwise returns false and \a x is not modified.
	virtual bool put(const T& x) {
		return write(&x, 1) == 1;
	}

	//@{
	/// Main write functions.
	/// Accept buffer \a buf to write to the stream and its size \a count.
	/// Return number of elements written.
	/// Usually it is \a count elements, 
	///  otherwise end of output stream (eos) has been reached.
	/// 
	/// Has two implementations: pure virtual and non-virtual.
	/// Pure virtual implementation should be defined in child classes.
	/// Non-virtual implementation is based on \ref put function 
	///  and need not to be defined in child classes.
	template<typename T1>
	size_t write(const T1 *buf, size_t count) {
		T x;
		size_t j;
		for(j = 0; j < count; j++) {
			x = buf[j];
			if(!put(x)) break;
		}
		return j;
	}

	virtual size_t write(const T *buf, size_t count) {
		return write<T>(buf, count);
	}
	//@}

};


///
/// Abstract input-output pipe.
/// Input-output pipe is output stream connected to input stream:
///  one side outputs, other side inputs corresponding values.
/// Should be used in functions to show *any* input/output pipes.
/// You can not construct iopipe object directly, 
/// you should construct its child classes instead.
/// 

template<typename T1, typename T2 = T1>
class iopipe 
{
public:

	//@{
	/// Get input and output streams.
	istream<T1>& input()  { return *_in;  }
	ostream<T2>& output() { return *_out; }
	operator istream<T1>&() { return *_in;  }
	operator ostream<T2>&() { return *_out; }
	//@}

protected: // is used by child classes

	void set_istream(istream<T1> *in) {
		this->_in = in;
	}

	void set_ostream(ostream<T2> *out) {
		this->_out = out;
	}

private:

	//@{
	/// Input and output stream references.
	/// Should be pointers to child class members!
	istream<T1> *_in;
	ostream<T2> *_out;
	//@}

};


} // namespace io 

#endif//_IO_BASE_
