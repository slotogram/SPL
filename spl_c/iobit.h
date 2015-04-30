#ifndef _IO_BITS_
#define _IO_BITS_

///
/// \file iobit.h
/// \brief Bit input and output wrappers.
///
/// This classes enable bit-by-bit input/output from/into other stream classes (memory, files, ...).
/// 

#include "iowrap.h"

namespace io {

////////////////////////////////////////////////////////////////////////
//                        BIT STREAM CLASSES                          //
////////////////////////////////////////////////////////////////////////

template<typename B, typename _Stream> class abstract_bitwrap;
template<typename B> class ibitwrap;
template<typename B> class obitwrap;

//@{
/// Standard IO-bit typedefs.
typedef ibitwrap<unsigned char>  ibitwrap8;
typedef obitwrap<unsigned char>  obitwrap8;
typedef ibitwrap<unsigned short> ibitwrap16;
typedef obitwrap<unsigned short> obitwrap16;
typedef ibitwrap<unsigned long>  ibitwrap32;
typedef obitwrap<unsigned long>  obitwrap32;
//@}


/// 
/// Abstract bit stream wrapper.
/// Contains members independent of bit stream wrapper direction.
/// Should not be used in functions, use ibitwrap, obitwrap and iobitwrap instead.
/// Should be used only as bit streams base class.
///

template<typename B, typename _Stream>
class abstract_bitwrap:
	virtual public abstract_wrap<_Stream>
{
public:
	//@{
	/// Bit stream position.
	/// It is abstract_stream::pos function implementation.
	virtual size_t pos() const { return _understream->pos() * maxbits() + _bits; }
	virtual size_t pos(size_t newpos) {
		return 0;
	}
	//@}

	/// Checks whether end of stream is reached.
	/// It is abstract_stream::eos function implementation.
	virtual bool eos() const { return _understream->eos() && _bits == maxbits(); }

protected:

	/// Abstract bit stream constructor - protected - class is to inherit only.
	abstract_bitwrap(_Stream& str): _buffer(0), _bits(0) {	}

	/// Type for bit count in internal buffer.
	typedef unsigned char bc_t;

	/// Buffer for storing of current position bits.
	/// Mutable for input streams.
	mutable B _buffer;

	/// Count of bits in internal buffer. 
	/// Mutable for input streams.
	mutable bc_t _bits;

	/// Maximum bit count in buffer.
	static bc_t maxbits() { return sizeof(B) * 8; }

	/// Bitmask.
	static B bitmask(bc_t bit) { return 1 << bit; }

}; // end of class abstract_bitwrap


///
/// Input bit stream.
/// This class can be used for object creation.
/// Object ibitwrap can be then passed as argument.
/// It is not recommended to require it in functions, 
///  it is better to require abstract input stream instead (\ref istream).
///

template<typename B = unsigned char>
class ibitwrap: 
	public abstract_bitwrap< B, istream<B> >, 
	public iwrap<bool, B>
{
public:
	/// Input bit stream constructor.
	/// Sets underlying input stream from arguments.
	ibitwrap(istream<B>& str): 
	  abstract_bitwrap< B, istream<B> >(str), iwrap<bool, B>(str) {}

	/// Get next bit.
	virtual bool get(bool &x) {
		if(_bits >= maxbits() || _bits == 0) {
			if(!_understream->get(_buffer)) return false;
			_bits = 0;
		}
		x = (_buffer & bitmask(_bits++)) != 0;
		return true;
	}
}; // end of class ibitwrap


///
/// Output bit stream.
/// This class can be used for object creation.
/// Object obitwrap can be then passed as argument.
/// It is not recommended to require it in functions, 
///  it is better to require abstract output stream instead (ostream).
///

template<typename B = unsigned char>
class obitwrap: 
	public abstract_bitwrap< B, ostream<B> >, 
	public owrap<bool, B>
{
public:
	/// Output bit stream constructor.
	/// Sets underlying output stream from arguments.
	obitwrap(ostream<B>& str): 
	  abstract_bitwrap< B, ostream<B> >(str), owrap<bool,B>(str) {}

	/// Output bit stream destructor.
	/// Flushes buffer into underlying stream.
	~obitwrap() {
		flush();
	}

	/// Put next bit.
	virtual bool put(const bool &x) {
		if(x) _buffer |= bitmask(_bits);
		if(++_bits >= maxbits()) {
			if(!_understream->put(_buffer)) return false;
			_buffer = 0;
			_bits = 0;
		}
		return true;
	}

	/// Flush bit buffer into output stream.
	/// Attention! 
	/// This function can cause null bits in the middle of stream, 
	///  if it is called before end of underlying output stream.
	bool flush() {
		if(_bits) { 
			size_t written = _understream->put(_buffer);
			_buffer = 0;
			_bits = 0;
			return written > 0;
		}
		return false;
	}

}; // end of class obitwrap


} // namespace io

#endif//_IO_BITS_
