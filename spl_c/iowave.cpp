


#include "iowave.h"
#include <limits.h>

namespace spl {

iwstream::iwstream(const char *file): 
	_f(file), F(0), N(0), M(0), B(0)
{
	word w; dword d, s;
	if(!get_element(d) || d != 0x46464952) // "RIFF"
		throw "Not a RIFF file";
	get_element(d); // chunk data size
	if(!get_element(d) || d != 0x45564157) // "WAVE"
		throw "Not a WAVE file";

	// go through chunks
	do {
		// chunk id && chunk size
		if(!get_element(d) || !get_element(s)) 
			throw "Can't read from file";
		if(d == 0x20746D66) { // "fmt "
			if(!get_element(w) || w != 1) // compression code
				throw "Don't support compressed wave-files";
			if(!get_element(M) // number of channels
			|| !get_element(F) // sample rate
			|| !get_element(d) // bytes per second
			|| !get_element(w) // block size
			|| !get_element(B))// bits per sample
				throw "Can't read from file";
			B = (B+7)/8; // bits per sample -> bytes per sample
		} else if(d == 0x61746164) { // "data"
			N = s; // chunk size
			_base = _f.pos(); // запоминаем базу секции отсчетов
			break; // found data-chunk, go out
		} else {
			if(!get_element(d)) // chunk size
				throw "Can't read from file";
			_f.skip(d-8); // skip this chunk
		}
	} while(!_f.eos());

	if(!F || !N || !M || !B) 
		throw "Key fields of wave-file were not read";

}

bool iwstream::get(signal_t& x) {
	byte b; short s; long d;
	switch(B) {
	case 1: return get_slice<byte>(b)  ? (convert_sample_to_float<byte>(b, x),  true) : false;
	case 2: return get_slice<short>(s) ? (convert_sample_to_float<short>(s, x), true) : false;
	case 4: return get_slice<long>(d)  ? (convert_sample_to_float<long>(d, x),  true) : false;
	default: return false;
	}
}

size_t iwstream::pos() const {
	return (_f.pos() - _base) / M;
}

size_t iwstream::pos(size_t n) {
	return (_f.pos(_base + n*M) - _base) / M;
}

size_t iwstream::skip(size_t n) {
	_f.skip(n*M);
	return pos();
}

bool iwstream::eos() const {
	return _f.pos() == _base + N;
}

signal iwstream::getInMemory()
{
	signal *ret = new signal;
	ret->F = F;
	ret->N = N/B;
	ret->X = new signal_t[ret->N];
	int i=0;
	while (!eos())
	{		
		iwstream::get(ret->X[i]);
		i++;
	}
	if (i<ret->N) throw  new _exception();
	return *ret;

}
} // namespace spl
