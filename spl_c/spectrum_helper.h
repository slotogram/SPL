#ifndef _SPL_SPECTRUM_HELPER_
#define _SPL_SPECTRUM_HELPER_

#include "spectrum.h"
#include "iowrap.h"
#include <algorithm>

namespace spl {

///
/// ���������, ������������ ��� ��������� �������� (\ref spectrum_filters).
///

struct filter_parameters {
	
	/// ������� ������������� �������.
	freq_t F;

	/// ��������, ������������ �������� ���������� (������ ���� ��������).
	double ksi;

};

/// 
/// ������������ �������������� �������, ������������� ��� ��������� �������.
/// 

struct spectrum_filters {

	/// ���������� ������� ����������.
	int K;
	
	/// ������ ���� ����������.
	/// ������ ���� ���������� ������������ ��� ��������� ������� � ������� �� ������������ �������� ksi � ������� �������������.
	/// ��� ������ ���� ����������, ��� ��������� �������� ������, �� ��� ������ ��� ����������.
	size_t Ws;

	/// ������ ���������� � �������� ������������� ����������.
	/// ������� ����������� (K, CONV_WIN_SIZ, 2)
	/// �������� � ���������� ������������ ���������� ��������������� ������������ ��� ����������.
	/// ��������� �������� � ���������� ������������ ������� �� ���� �� pi/2, 
	/// �� ����� ������ ��� ����� � ���� ����������, 
	/// �������� �� �� ��� ��������� ������������ ����� ���������� � ����.
	double *H;

};

/// 
/// ��������� ������������� ���������� �� ����� ����������� ������ ��������.
/// � ��������� �������� \a f ��� ������ ���� �������� ����� ��� ��������� ������������� ����������.
/// 

bool spectrum_filters_generate(
	IN freq_scale& s,        ///< ����� ����������� ������ ��������
	OUT spectrum_filters& f, ///< ������������ ����������
	IN filter_parameters& p  ///< ��������� ��������� ��������
);


///
/// �� ��, ��� spectrum_filters, 
/// ������ � ����������� ������� (��� �++)
///

struct EXPORT spectrum_filters_class: spectrum_filters 
{
	/// �����������.
	/// ������ ����� ������� filters_generate(), 
	///  �������������� ������� ����� ��� �������.
	/// ��� ���������� ����� �������������.
	spectrum_filters_class(
		const freq_scale& s, 
		const filter_parameters& p);
	/// ��������� ������������ �������� � ����.
	bool save(const char *file);
	~spectrum_filters_class();
};

///
/// �������, ������������ ����� ������� �� \a N ������� ��������
///

template<typename T>
class iwstream_extend: public io::iwrap<T> {
public:

	iwstream_extend(io::istream<T>& str, size_t N):
	  io::iwrap<T>(str), _N(N), _pos(0) {}

	//@{
	/// ����������� ������� abstract_stream.
	virtual bool eos() const { return _pos >= _N; }
	virtual size_t pos() const { return _understream->pos() + _pos; }
	virtual size_t pos(size_t newpos) { return 0; } // ���� �� ���������� �� �������������
	//@}

	/// ��������� �������� 
	virtual bool get(T& x) {
		if(_understream->eos()) {
			if(!eos()) {
				x = 0;
				_pos++;
			} else {
				return false;
			}
		} else {
			_understream->get(x);
		}
		return true;
	}

	/// ��������� ������ ���������
	virtual size_t read(T *buf, size_t count) {
		size_t r1 = _understream->read(buf, count);
		size_t r2 = 0;
		if(r1 < count) {
			r2 = std::min(count-r1, _N-_pos);
			std::fill(buf + r1, buf + r1 + r2, 0);
			_pos += r2;
		}
		return r1 + r2;
	}

private:
	size_t _pos, _N;

};

} // namespace spl

#endif//_SPL_SPECTRUM_HELPER_