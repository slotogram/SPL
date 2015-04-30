#include "common.h"
#include "iobuf.h"
#include "membuf.h"
using namespace io;

#include <windows.h>

#include <queue>
using std::queue;

template<typename T>
class sync_queue
{
public:
	sync_queue() {
		InitializeCriticalSection(&critical_section);
	}

	~sync_queue() {
		DeleteCriticalSection(&critical_section);
	}

	void push(T& x) {
		enter();
		impl.push(x);
		leave();
	}

	T pop() {
		enter();
		T x = impl.front();
		impl.pop();
		leave();
		return x;
	}

	size_t size() const {
		return impl.size();
	}

private:
	std::queue<T> impl;
	mutable CRITICAL_SECTION critical_section;

	void enter() const {
		EnterCriticalSection(&critical_section);
	}

	void leave() const {
		LeaveCriticalSection(&critical_section);
	}
};


typedef unsigned char byte;

struct buffer_t {
	membuf<byte> data;
	size_t full_size;
	size_t fill_size;
};

typedef sync_queue<buffer_t> bufqueue_t;

class abstract_bufstream:
	virtual public abstract_stream
{
public:
	abstract_bufstream(bufqueue_t& filled_, bufqueue_t& free_, bool& closed):
	  filled_buffers(filled_), free_buffers(free_), _closed(closed), _pos(0)
	  {}

	virtual size_t pos() const {
		return _pos;
	}

	virtual size_t pos(size_t newpos) {
		throw "Not implemented";
	}

	virtual size_t skip(size_t N) {
		throw "Not implemented";
	}

	void set_buffer_filled_event(HANDLE event_) {
		buffer_filled_event = event_;
	}

protected:
	bufqueue_t& filled_buffers;
	bufqueue_t& free_buffers;
	bool& _closed;
	size_t _pos;
	HANDLE buffer_filled_event;
};

class obuf_uni:
	public ostream<byte>,
	public abstract_bufstream
{
public:
	obuf_uni(bufqueue_t& filled_, bufqueue_t& free_, bool& closed_, size_t bufsize_):
	  abstract_bufstream(filled_, free_, closed_), bufsize(bufsize_)
	{
	}

	virtual size_t write(const byte *data, size_t count);

	virtual bool eos() const {
		return _closed;
	}
	virtual size_t pos() const {
		return _pos;
	}
	virtual void close() {
		_closed = true;
	}

private:
	size_t bufsize;
};

size_t obuf_uni::write(const byte *data, size_t count) {
	if(_closed) return 0;
	size_t written = 0;
	while(written < count) {
		// get some free buffer
		buffer_t buffer;
		if(free_buffers.size() > 0) {
			buffer = free_buffers.pop();
		} else {
			buffer.data = membuf<byte>(bufsize);
			buffer.full_size = bufsize;
		}

		// fill it
		size_t bytes_to_copy = min(count - written, buffer.full_size);
		memcpy(buffer.data, data + written, bytes_to_copy);
		buffer.fill_size = bytes_to_copy;
		written += bytes_to_copy;

		// add to filled buffers
		filled_buffers.push(buffer);
		if(filled_buffers.size() == 1) {
			SetEvent(buffer_filled_event);
		}
	}
	_pos += written;
	return written;
}

class ibuf_uni:
	public istream<byte>,
	public abstract_bufstream
{
public:
	ibuf_uni(bufqueue_t& filled_, bufqueue_t& free_, bool& closed_, int max_):
	  abstract_bufstream(filled_, free_, closed_), buffer(), bufpos(0), max_buffers(max_)
	{
	}
	
	virtual size_t read(byte *data, size_t count);

	virtual bool eos() const {
		return _closed && buffer.data.get() == NULL && filled_buffers.size() == 0;
	}

	virtual void close() {
		_closed = true;
	}

private:
	buffer_t buffer;
	size_t bufpos;
	int max_buffers;
};

size_t ibuf_uni::read(byte *data, size_t count) {
	size_t read = 0;
	while(read < count) {
		// if current buffer is empty
		if(buffer.data == NULL) {
			// if there are no filled buffers
			while(filled_buffers.size() == 0) {
				if(_closed) {
					goto out;
				}
				// wait for the first one
				WaitForSingleObject(buffer_filled_event, INFINITE);
				ResetEvent(buffer_filled_event);
			}
			// get filled buffer
			buffer = filled_buffers.pop();
			bufpos = 0;
		}

		// read from it
		size_t bytes_to_copy = min(count - read, buffer.fill_size - bufpos);
		memcpy(data + read, buffer.data + bufpos, bytes_to_copy);
		bufpos += bytes_to_copy;
		read += bytes_to_copy;

		// if whole buffer is read - free it
		if(buffer.fill_size == bufpos) {
			if(free_buffers.size() < max_buffers) {
				// add it to free buffers
				free_buffers.push(buffer);
			} else {
				// just free resources
				buffer.data.reset();
			}
		}
	}
	out:
	_pos += read;
	return read;
}

// unidirectional iobuf
class iobuf_uni:
	public iobuf_impl
{
public:
	iobuf_uni():
	  max_buffers(4), bufsize(500), _closed(false),
	  _in(filled_buffers, free_buffers, _closed, max_buffers),
	  _out(filled_buffers, free_buffers, _closed, bufsize)
	{
		buffer_filled_event = CreateEvent(NULL, TRUE, FALSE, NULL);
		_in.set_buffer_filled_event(buffer_filled_event);
		_out.set_buffer_filled_event(buffer_filled_event);
	}

	~iobuf_uni() {
		CloseHandle(buffer_filled_event);
	}

	virtual istream<byte>& input() {
		return _in;
	}
	virtual ostream<byte>& output() {
		return _out;
	}

private:
	int max_buffers;
	size_t bufsize;
	HANDLE buffer_filled_event;

	bufqueue_t filled_buffers;
	bufqueue_t free_buffers;
	bool _closed;

	ibuf_uni _in;
	obuf_uni _out;
};

iobuf_impl *iobuf_impl::create() {
	return new iobuf_uni();
}

void iobuf_impl::destroy(iobuf_impl *impl) {
	delete impl;
}
