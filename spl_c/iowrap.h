#ifndef _IO_WRAP_
#define _IO_WRAP_

///
/// \file  iowrap.h
/// \brief Abstract input and output wrappers classes.
///
/// This module defines input and output wrappers.
/// Input wrapper (\ref iwrap) is an input stream that wraps another input stream 
///  and performs extra operations (i.e. preprocessing) before passing data to user.
/// Output wrapper (\ref owrap) is an output stream that wraps another output stream
///  and performs extra operations (i.e. postprocessing) before passing data to output.
///
/// I/O wrappers can be used to implement preprocessing and postprocessing of data, 
///  e.g. bitwise I/O, type conversion, format conversion and so on.
///

#include "io.h"


namespace io {

////////////////////////////////////////////////////////////////////////
//                        BIT STREAM CLASSES                          //
////////////////////////////////////////////////////////////////////////

template<typename _Stream> class abstract_wrap;
template<typename T, typename T2> class iwrap;
template<typename T, typename T2> class owrap;


///
/// Abstract stream wrapper.
/// Contains members independent of stream wrapper direction.
/// Should not be used in functions, use \ref iwrap and \ref owrap instead.
/// Should be used only as stream wrappers base class.
///

template<typename _Stream>
class abstract_wrap: 
 virtual public abstract_stream
{
public:

 //@{
 /// Common abstract_stream members - route to underlying input stream by default.
 virtual bool   eos() const { return _understream->eos(); }
 virtual size_t pos() const { return _understream->pos(); }
 virtual size_t pos(size_t newsize) { return _understream->pos(newsize); }
 virtual void close() { _understream->close(); }
 //@}

protected:

 /// Protected default constructor.
 /// It is necessary only for virtual inheritance, as it requires default constructor (in VC++).
 abstract_wrap(): _understream(0) {}

 /// the underlying stream.
 _Stream *_understream;

};


///
/// Input stream abstract wrapper.
/// It is an input stream that wraps another input stream 
///  and performs extra operations (i.e. preprocessing) before passing data to user.
/// Template parameter T2 - is a type, which is got from input.
/// Template parameter T - is a type, which is passed to user.
///

template<typename T, typename T2 = T>
class iwrap: 
 virtual public abstract_wrap< istream<T2> >,
 public istream<T>
{
protected:
	iwrap(istream<T2>& in_str) {
		_understream = &in_str;
	}
};

/// Input wrapper, that converts (static_cast) wrapped stream input element-by-element.
template<typename T, typename T2>
class iwrapelem:
	public iwrap<T, T2>
{
public:
	typedef void (*convert_func)(const T2& x, T& y);

	iwrapelem(istream<T2>& in_str):
	  iwrap<T, T2>(in_str), convert(default_convert) 
	{}
	
	iwrapelem(istream<T2>& in_str, convert_func f):
	  iwrap<T, T2>(in_str), convert(f) 
	{}
	
	/// istream::get() implementation in \ref iwrapelem - convert stream input element-by-element.
	virtual bool get(T& x) {
		T2 x2; 
		bool r = _understream->get(x2); 
		if(r) convert(x2, x);
		return r;
	}
private:
	convert_func convert;
	
	static void default_convert(const T2& x, T& y) {
		y = static_cast<T>(x);
	}
};

/// Input wrapper, that inputs into buffer of type T[] as if it was buffer of type T2[] (reinterpret_cast).
/// It can be used to implement some input behaviour for one element type (T2), and then use it for any type necessary.
template<typename T, typename T2>
class iwrapblock:
	public iwrap<T, T2>
{
public:

	iwrapblock(istream<T2>& in_str):
	  iwrap<T, T2>(in_str) {}
	
	virtual size_t read(T *block, size_t count) {
		size_t count_translated = count * sizeof(T) / sizeof(T2);
		size_t read = _understream->read( reinterpret_cast<T2 *>(block), count_translated);
		size_t read_translated = read * sizeof(T2) / sizeof(T);
		return read_translated;
	}

};


///
/// Output stream abstract wrapper.
/// It is an output stream that wraps another output stream
///  and performs extra operations (i.e. postprocessing) before passing data to output.
/// Template parameter T - is a type, which is got from user.
/// Template parameter T2 - is a type, which is passed to output.
///

template<typename T, typename T2 = T>
class owrap: 
 virtual public abstract_wrap< ostream<T2> >,
 public ostream<T>
{
protected:
	owrap(ostream<T2>& out_str) {
		_understream = &out_str;
	}
};

/// Output wrapper, that converts (static_cast) wrapped stream output element-by-element.
template<typename T, typename T2>
class owrapelem:
	public owrap<T, T2>
{
public:
	typedef void (*convert_func)(const T& x, T2& y);

	owrapelem(ostream<T2>& out_str):
	  owrap<T, T2>(out_str), convert(default_convert) {}
	
	owrapelem(ostream<T2>& out_str, convert_func f):
	  owrap<T, T2>(out_str), convert(f) {}
	
	/// ostream::put() implementation in \ref owrapelem - convert stream output element-by-element.
	virtual bool put(const T& x) {
		T2 y;
		convert(x, y);
		return _understream->put(y);
	}
private:
	convert_func convert;

	static void default_convert(const T& x, T2& y) {
		y = static_cast<T2>(x);
	}
};

/// Output wrapper, that outputs into buffer of type T[] as if it was buffer of type T2[] (reinterpret_cast).
/// It can be used to implement some output behaviour for one element type (T2), and then use it for any type necessary.
template<typename T, typename T2>
class owrapblock:
	public owrap<T, T2>
{
public:

	owrapblock(ostream<T2>& out_str):
	  owrap<T, T2>(out_str) {}
	
	virtual size_t write(const T *block, size_t count) {
		size_t count_translated = count * sizeof(T) / sizeof(T2);
		size_t written = _understream->write( reinterpret_cast<const T2 *>(block), count_translated);
		size_t written_translated = written * sizeof(T2) / sizeof(T);
		return written_translated;
	}

};

} // namespace io

#endif//_IO_WRAP_
