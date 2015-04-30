#include "iofile.h"

#ifdef _DEBUG
# include <stdio.h>
#endif

//*
#define _IO_FILE_STDIO_
/*/
#define _IO_FILE_WINDOWS_
//*/

#ifdef _IO_FILE_WINDOWS_
# include <windows.h>
# define READ_MODE GENERIC_READ, FILE_SHARE_WRITE, 0, OPEN_EXISTING, 0, 0
# define WRITE_MODE GENERIC_WRITE, FILE_SHARE_READ, 0, CREATE_ALWAYS, 0, 0
# define OPEN_FILE CreateFile
# define CLOSE_FILE CloseFile
# define INV_FILE INVALID_HANDLE_VALUE
# define FH (HANDLE)_data
#endif

#ifdef _IO_FILE_STDIO_
# include <stdio.h>
# include <io.h>
# define READ_MODE "rb"
# define WRITE_MODE "wb"
# define OPEN_FILE fopen
# define CLOSE_FILE fclose
# define INV_FILE 0
# define FH (FILE*)_data
#endif



namespace io {

abstract_fstream::abstract_fstream(
	const char *filename, 
	stream_mode mode)
{
	switch(mode) {
	case io::in:
		_data = OPEN_FILE(filename, READ_MODE);
		break;
	case io::out:
		_data = OPEN_FILE(filename, WRITE_MODE);
		break;
	default:
		throw "Invalid mode";
	}
	if(_data == INV_FILE) {
#ifdef _DEBUG
		printf("abstract_fstream: Can not open file '%s' for %s\n", filename, mode == io::in ? "input" : "output");
#endif
		throw "Can not open file";
	}
}

void abstract_fstream::close() {
	CLOSE_FILE(FH);
}

#ifdef _IO_FILE_WINDOWS_

size_t abstract_fstream::size() const {
	LARGE_INTEGER r;
	if(GetFileSizeEx(FH, &r)) 
		return (size_t)r.QuadPart;
	return 0;
}

size_t abstract_fstream::pos() const {
	LARGE_INTEGER _null, _pos;
	_null.QuadPart = 0;
	if(SetFilePointerEx(FH, _null, &_pos, FILE_CURRENT))
		return (size_t)_pos.QuadPart;
	return 0;
}

size_t abstract_fstream::pos(size_t N) {
	LARGE_INTEGER _new, _pos;
	_new.QuadPart = N;
	if(SetFilePointerEx(FH, _new, &_pos, FILE_BEGIN))
		return (size_t)_pos.QuadPart;
	return pos();
}

size_t abstract_fstream::skip(size_t N) {
	LARGE_INTEGER _new, _pos;
	_new.QuadPart = N;
	if(SetFilePointerEx(FH, _new, &_pos, FILE_CURRENT))
		return (size_t)_pos.QuadPart;
	return pos();
}

bool abstract_fstream::eos() const {
	return pos() >= size();
}

size_t abstract_fstream::_read(void *buf, size_t bytes) {
	DWORD read;
	ReadFile(FH, buf, (DWORD)bytes, &read, 0);
	return read;
}

size_t abstract_fstream::_write(const void *buf, size_t bytes) {
	DWORD written;
	WriteFile(FH, buf, (DWORD)bytes, &written, 0);
	return written;
}

#endif // _IO_FILE_WINDOWS_

#ifdef _IO_FILE_STDIO_

size_t abstract_fstream::size() const {
	return (size_t)_filelengthi64(_fileno(FH));
}

size_t abstract_fstream::pos() const {
	fseek(FH, 0, SEEK_CUR);
	return ftell(FH);
}

size_t abstract_fstream::pos(size_t newpos) {
	fseek(FH, (long)newpos, SEEK_SET);
	return ftell(FH);
}

size_t abstract_fstream::skip(size_t N) {
	fseek(FH, (long)N, SEEK_CUR);
	return ftell(FH);
}

bool abstract_fstream::eos() const {
	return feof(FH) != 0;
}

size_t abstract_fstream::_read(void *buf, size_t bytes) {
	return fread(buf, 1, bytes, FH);
}

size_t abstract_fstream::_write(const void *buf, size_t bytes) {
	return fwrite(buf, 1, bytes, FH);
}

#endif // _IO_FILE_STDIO_


} // namespace io
