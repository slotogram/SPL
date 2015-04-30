#ifndef _IO_WAVE_
#define _IO_WAVE_

#include "iofile.h"
#include "spl_types.h"
#include <limits>

namespace spl {

/// 
/// ����� wave-�����.
/// ��������� wave-���� � ������������ � ���� ��������� io::istream<signal_t>.
/// ������������� ����� ������� freq() ��� ��������� ������� �������������.
///

class EXPORT iwstream:
	public io::istream<signal_t>
{
public:
	/// �����������. ��������� ���� � wav-�����. 
	iwstream(const char *file);
	/// �������� ��������� ������.
	/// � �������������� �������� (��������, ������) ������ ������� �������� ����� �� ������� �� ������.
	virtual bool get(signal_t& x);
	/// ������� ������� � ������. ���������� �������� �����.
	virtual size_t pos() const;
	/// ���������� ������� � ������.
	virtual size_t pos(size_t N);
	/// ���������� \a N ��������.
	virtual size_t skip(size_t N);
	/// ��������� ����� ������.
	virtual bool eos() const;

	/// �������� ������� �������������.
	freq_t freq() const { return (freq_t)F; }
	//��������� ������ � ������.
	signal getInMemory();

	/// ������� �����.
	virtual void close() { _f.close(); }

private:
	typedef unsigned char byte;
	typedef unsigned short word;
	typedef unsigned long dword;
	io::ifstream<byte> _f;

	template<typename T>
	bool get_element(T& x) {
		return _f.read((byte*)&x, sizeof(T)) > 0;
	}

	// �������, ���������� ���� ��������� ����� ������ �� ���� �������.
	// �������� ��� ������ - ����� �� ���� �������.
	template<typename T>
	bool get_slice(T& x);
	
	// WAVE FILE FORMAT
	size_t _base; ///< ���� ������ ��������
	dword N; ///< ������ ������ ��������
	word M; ///< ���������� ����� �� ��������� ������� (������ - ����, ������)
	word B; ///< ���������� ���� � �����
	dword F; ///< ������� �������������
};

template<typename T>
bool iwstream::get_slice(T& x) {
	x = 0; 
	T y;
	for(int i = 0; i < M; i++) {
		if(!get_element(y)) return false;
		x += y;
	}
	return true;
}

template<typename T>
void convert_sample_to_float(const T& x, signal_t& y) {
	T m = std::numeric_limits<T>::min();
	T M = std::numeric_limits<T>::max();
	y = ( 2*signal_t(x) - (signal_t(M) + 1 + m) ) / (signal_t(M) + 1 - m);
}

template<typename T>
void convert_float_to_sample(const signal_t& x, T& y) {
	// TODO implement this
}

} // namespace spl

#endif//_IO_WAVE_
