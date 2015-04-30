#ifndef _SPL_MASK_HELPER_
#define _SPL_MASK_HELPER_

#include "iofile.h"
#include "mask.h"
#include "conv.h"
#include "iowrap.h"
#include "spl_types.h"
#include <algorithm>

namespace spl {

//
// Forward declarations
//

/// ������������ �������� ����������
struct mask_filters;

/// ���������� ������������� ������������� ���������� �� ����� ����������� ������ ��������.
bool mask_filters_generate(const freq_scale& s, mask_filters& f, const mask_parameters& p);

/// ���������� ������������� ������������� ����������, ���������������� ��� �������� ����������.
bool mask_filters_generate_fast(const freq_scale& s, mask_filters& f, const mask_parameters& p);

/// ��������� ������������� ����������.
size_t EXPORT mask_stream(IN mask_filters& f, io::istream<spectrum_t>& in_str, io::ostream<mask_t>& out_str);

/// ��������� ������������� ���������� - ������� �������.
size_t EXPORT mask_stream_fast(IN mask_filters& f, io::istream<spectrum_t>& in_str, io::ostream<mask_t>& out_str);



// ����� input-wrapper ��� ���������� ����� ������.
template<typename T>
class istream_block_extend;

// ����� output-wrapper ��� ������� ����� ������.
template<typename T> 
class ostream_block_reduce;

/// ������������� ��� ���������� ����� ������ � �������.
typedef istream_block_extend<spectrum_t> istream_spectrum_extend;

/// �������������� ��� ������� ����� ������ � �����.
typedef ostream_block_reduce<mask_t> ostream_mask_reduce;


/// ������������ �������� ����������

struct mask_filters {
	
	/// ���������� ������� � �������.
	int K; 

	/// ������ ���� ��������.
	int Ws; 

	/// ������ ������������� �������� ����������.
	/// ������� ����������� (K, Ws).
	double *H;

	/// ����� �����, �� ������� ������� ������������� ����� ������.
	scale_form sc_form;
};

///
/// �� ��, ��� mask_filters, 
/// ������ � ����������� ������� (��� �++)
///

struct EXPORT mask_filters_class: mask_filters 
{
	/// �����������.
	/// ������ ����� ������� mask_filters_generate(), 
	///  �������������� ������� ����� ��� �������.
	/// ��� ���������� ����� �������������.
	mask_filters_class(
		const freq_scale& s, 
		const mask_parameters& p);
	~mask_filters_class() { conv_free(H); }

	/// ��������� ������������ � ����.
	bool save(const char *file) { return io::array_to_file(H, K * Ws, file); }

	/// ��������� ���������� � ����������� �� ����� �����.
	/// ��� ��������� ����� ��������� ������� ��������� ���������� �����.
	/// ��� ��������� ���� ����� ��������� �������, ������������������ ������.
	size_t mask_stream(io::istream<spectrum_t>& in_str, io::ostream<mask_t>& out_str) {
		if(allow_fast) 
			return mask_stream_fast(*this, in_str, out_str);
		else 
			return spl::mask_stream(*this, in_str, out_str);
	}
private:
	bool allow_fast;
};



///
/// ������������� ������: ���������� ������ ������.
///
/// ��� ��������� ������������� ��������� ������������� ����������
///  �� ���� ���������� ������� ������� � ��� ������� �� �������� ������� ����.
///

template<typename T>
class istream_block_extend: public io::iwrap<T> 
{
public:
	/// �����������.
	/// ��������� ������������� �����, ������ ����� � ������ ����������.
	istream_block_extend(io::istream<T>& str, int _K, int _W):
	  io::iwrap<T>(str), K(_K), W(_W), _pos(0), _maxpos(0) { 
		// ����� _maxpos = _pos = 0, �� ���������� ��������� ������ ����
		_buf = spl_alloc<T>(W + K + W);
	}

	/// ����������.
	/// ����������� ���������� ������
	~istream_block_extend() { spl_free(_buf); }

	//@{
	/// ����������� ������� abstract_stream.
	virtual bool eos() const { return _pos >= _maxpos && _understream->eos(); }
	virtual size_t pos() const { return _understream->pos() / K * (W+K+W) + _pos; }
	virtual size_t pos(size_t newsize) { return 0; } // ���� �� ���������� �� �������������
	//@}

	/// ��������� ��������.
	virtual bool get(T& x) {
		if(_pos >= _maxpos) {
			if(!read_block()) return false;
		}
		x = _buf[_pos++];
		return true;
	}
//*
	/// ��������� ������ ���������.
	/// ���������� ������� istream::read()
	virtual size_t read(T *buf, size_t count) {
		size_t i = 0;
		// ���� ���� ���� ������
		while(i < count) {
			// ���� �������� ��������� - ������ ��������� ����
			if(_pos >= _maxpos) {
				if(!read_block()) break; // ���� �������� � ��� ��������� - �������
			}
			// �������� ����������� ��: ����������� � _buf � ����������� � buf
			size_t N = std::min(_maxpos - _pos, count - i);
			size_t r = std::copy(_buf + _pos, _buf + _pos + N, buf + i) - buf - i;
			_pos += r; i += r;
		}
		return i; // ���������� ���������� ���������� ��������
	}
//*/
private:

	bool read_block() {
		// ������ ����: K ��������� - ����������, ������� �����������
		_maxpos = _understream->read(_buf + W, K);
		// ���� ������ �� ����������� - ����� �����
		if(_maxpos == 0) return false;
		// ����������� _maxpos ���������, ������ ���������� �������:
		std::fill(_buf, _buf + W, _buf[W]);
		T *buf2 = _buf + W + _maxpos;
		std::fill(buf2, buf2 + W, buf2[-1]);
		_pos = 0;
		_maxpos = W + _maxpos + W;
		return true;
	}

	/// ������ �����.
	const int K;
	/// ������ ����������.
	const int W;
	/// ���������� ����� � ����������� ��������.
	T *_buf;
	/// ������� � ������.
	size_t _pos;
	/// ������������ ������� � ������.
	size_t _maxpos;
};

///
/// �������������� ������: ������������ ��������.
/// ������������ � ������� ��������� ������������� ����������
///  ��� �������� �� ����������� ����� ������ � �������.
///

template<typename T>
class ostream_block_reduce: public io::owrap<T>
{
public:
	/// �����������.
	/// ��������� �������� �����, ������ ����� � ������ ����������.
	ostream_block_reduce(io::ostream<T>& str, size_t _K, size_t _W):
	  io::owrap<T>(str), K(_K), W(_W), rel_pos(0) {}

	/// ������� �������.
	bool put(const T& x) {
		bool r = true;
		if(rel_pos >= W && rel_pos < W + K) {
			r = _understream->put(x);
		}
		rel_pos = (rel_pos + 1) % (W + K + W);
		return r;
	}
//*
	/// ������� ������ ���������.
	size_t write(const T *buf, size_t count) {
		size_t i = 0, j = 0;
		while(i < count) {
			// ���� ������� �� ���, ��� �����
			size_t N;
			if(rel_pos < W) {
				N = std::min(W - rel_pos, count - i);
			} else if(rel_pos >= K+W) {
				N = std::min(W+K+W+W - rel_pos, count - i);
			} else {
				N = 0;
			}
			i += N;
			rel_pos = (rel_pos + N) % (W + K + W);
			if(i >= count) break;
			N = std::min(W + K - rel_pos, count - i);
			size_t M = _understream->write(buf + i, N);
			i += M; j += M;
			// ���� �������� ������, ��� ����� - �������
			if(M < N) break;
			rel_pos += M;
		}
		return j; // ������� ����� ���������� �� ����� ���������.
	}
//*/
private:
	size_t rel_pos;
	/// ������ �����.
	const size_t K;
	/// ������ ����������.
	const size_t W;
};


} // namespace spl

#endif//_SPL_MASK_HELPER_