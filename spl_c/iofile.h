#ifndef _IO_FILES_
#define _IO_FILES_

///
/// \file  iofile.h
/// \brief File input and/or output streams.
///

#include "io.h"

namespace io {

////////////////////////////////////////////////////////////////////////
//                       FILE STREAM CLASSES                          //
////////////////////////////////////////////////////////////////////////

class abstract_fstream;
template<typename T> class  ifstream;
template<typename T> class  ofstream;
template<typename T> class iofstream;

/// 
/// Abstract file stream.
/// Contains members independent of file stream element data type.
/// Should not be used in functions, use ifstream, ofstream and iofstream instead.
/// Should be used only as file streams base class.
///

class EXPORT abstract_fstream: 
	 virtual public abstract_stream 
{
public:
	/// File stream size.
	virtual size_t size() const;

	//@{
	/// Functions pos and skip - work with file stream position.
	virtual size_t pos() const;
	virtual size_t pos(size_t newpos);
	virtual size_t skip(size_t N);
	//@}

	/// Check end of file stream.
	virtual bool eos() const;

	/// Close file stream.
	virtual void close();

protected:

	/// Protected constructor - only for inheritance
	abstract_fstream(const char *filename, stream_mode mode);

	/// Destructor - automatically closes the file.
	~abstract_fstream() { close(); }

	/// Raw read.
	size_t _read(void *buf, size_t bytes);

	/// Raw write.
	size_t _write(const void *buf, size_t bytes);

	/// File handle.
	void *_data;

};

/// Input file stream.
/// Can be used for direct object instantiation.

template<typename T>
class ifstream: 
	public abstract_fstream, 
	public istream<T>
{
public:

	/// Constructor.
	ifstream(const char *filename):
	  abstract_fstream(filename, io::in) {}

	/// Get one element.
	bool get(T& x) {
		return _read(&x, sizeof(T)) == sizeof(T);
	}

	//@{
	/// Functions size, pos, skip, read - same as in \ref abstract_fstream class, 
	///  except for scaling by type size.
	virtual size_t size() const { return abstract_fstream::size() / sizeof(T); }
	virtual size_t pos () const { return abstract_fstream::pos () / sizeof(T); }
	virtual size_t pos (size_t n) { return abstract_fstream::pos (n*sizeof(T)) / sizeof(T); }
	virtual size_t skip(size_t n) { return abstract_fstream::skip(n*sizeof(T)) / sizeof(T); }

	virtual size_t read(T *buf, size_t count) {
		return _read(buf, count * sizeof(T)) / sizeof(T);
	}
	//@}

};

/// Output file stream.
/// Can be used for direct object instantiation.

template<typename T>
class ofstream:
	public abstract_fstream, 
	public ostream<T>
{
public: 

	/// Constructor.
	ofstream(const char *filename):
	  abstract_fstream(filename, io::out) {}

	/// Put one element.
	bool put(const T& x) {
		return _write(&x, sizeof(T)) == sizeof(T);
	}

	//@{
	/// Functions size, pos, skip, write - same as in \ref abstract_fstream class, 
	///  except for scaling by type size.
	virtual size_t size() const { return abstract_fstream::size() / sizeof(T); }
	virtual size_t pos () const { return abstract_fstream::pos () / sizeof(T); }
	virtual size_t pos (size_t n) { return abstract_fstream::pos (n*sizeof(T)) / sizeof(T); }
	virtual size_t skip(size_t n) { return abstract_fstream::skip(n*sizeof(T)) / sizeof(T); }

	virtual size_t write(const T *buf, size_t count) {
		return _write(buf, count * sizeof(T)) / sizeof(T);
	}
	//@}

	/// End of stream check redefinition. 
	/// Output to file never ends.
	virtual bool eos() const { return false; }

};

/// Helper template function array_to_file.
/// Writes any array to specified file.

template<typename T>
bool array_to_file(const T *x, size_t N, const char *file) {
	try {
		io::ofstream<T> out(file);
		return out.write(x, N) == N;
	} catch(...) {
		return false;
	}
}

/// Helper template function array_from_file.
/// Reads any array from specified file.

template<typename T>
bool array_from_file(T *x, size_t N, const char *file) {
	try {
		io::ifstream<T> in(file);
		return in.read(x, N) == N;
	} catch(...) {
		return false;
	}
}

} // namespace io

#endif//_IO_FILES_
