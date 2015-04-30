#ifndef _IO_MIC_
#define _IO_MIC_

#include "common.h"
#include "io.h"
#include "iowrap.h"
#include "iobuf.h"
#include <mmsystem.h>
#include <Windows.h>

namespace io {

class EXPORT mic_writer_impl {
public:
	static mic_writer_impl *create(ostream<byte>& output, int elemsize, int sample_rate);
	static void destroy(mic_writer_impl *);

	virtual void start() = 0;
	virtual void stop()  = 0;
	virtual void close() = 0;
};

template<typename T>
class mic_writer {
public:
	mic_writer(ostream<T>& output, int sample_rate): 
	  _output(output),
	  _impl(mic_writer_impl::create(_output, sizeof(T), sample_rate))
	{
	}

	~mic_writer() {
		mic_writer_impl::destroy(_impl);
	}

	virtual void start() { _impl->start(); }
	virtual void stop()  { _impl->stop();  }
	virtual void close() { _impl->close(); }

private:
	typedef unsigned char byte;
	owrapblock<byte,T> _output;
	mic_writer_impl *_impl;
};

}

#endif//_IO_MIC_