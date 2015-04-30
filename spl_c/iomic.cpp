#include <windows.h>
#include <mmsystem.h>
#include "iomic.h"
#include "iofile.h"

using namespace io;

class mic_writer_winmm_impl:
	public mic_writer_impl
{
public:

	mic_writer_winmm_impl(ostream<byte>& output, int elemsize, int sample_rate):
	  _output(output), _elemsize(elemsize), _freq(sample_rate), _bufsize(BUFSIZE), _closed(false)
	{
		thread = CreateThread(NULL, 10000, threadProc, this, 0, &threadId);
		initialized_event = CreateEvent(NULL, TRUE, FALSE, NULL);
		deinitialized_event = CreateEvent(NULL, TRUE, FALSE, NULL);
	}

	~mic_writer_winmm_impl() {
		close();
		WaitForSingleObject(deinitialized_event, INFINITE);
		CloseHandle(thread);
		CloseHandle(initialized_event);
		CloseHandle(deinitialized_event);
	}

	// methods, that are called in creator thread
	virtual void start() { send(MIC_START); }
	virtual void stop () { send(MIC_STOP ); }
	virtual void close() { send(MIC_CLOSE); }

private:

	// methods, that are called in worker thread
	void _open();
	void _init(WAVEHDR&);
	void _start();
	void _process(WAVEHDR *);
	void _stop();
	void _deinit(WAVEHDR&);
	void _close();

	// thread message constants
	enum ThreadMsg {
		MIC_START = WM_USER + 1,
		MIC_DATA,
		MIC_STOP,
		MIC_CLOSE
	};

	// send a message to worker thread
	void send(UINT uMsg) {
		if(_closed) return;
		WaitForSingleObject(initialized_event, INFINITE);
		PostThreadMessage(threadId, uMsg, 0, 0);
	}

	ostream<byte>& _output;
	int _elemsize;
	int _freq;
	int _bufsize;
	bool _closed;

	// thread details
	HANDLE thread;
	HANDLE initialized_event;
	HANDLE deinitialized_event;
	DWORD threadId;

	HWAVEIN hWaveIn;
	// buffers
	static const int NUM_BUFFERS = 5;
	static const int BUFSIZE = 1024;
	WAVEHDR headers[NUM_BUFFERS];

	static DWORD WINAPI threadProc(LPVOID lpParam);

#ifdef _DEBUG
	HMMIO _debug;
#endif
};

DWORD WINAPI mic_writer_winmm_impl::threadProc(LPVOID lpParam) {
	mic_writer_winmm_impl *that = static_cast<mic_writer_winmm_impl *>(lpParam);
	that->_open();

	MSG msg;

	// prepare message queue
	PeekMessage(&msg, NULL, WM_USER, WM_USER, PM_NOREMOVE);

	// set initialized event
	SetEvent(that->initialized_event);

	// message queue:
	while(GetMessage(&msg, NULL, 0, 0)) {
		switch(msg.message) {
		case MM_WIM_DATA:
			that->_process(reinterpret_cast<WAVEHDR *>(msg.lParam));
			break;
		case MIC_START:
			that->_start();
			break;
		case MIC_STOP:
			that->_stop();
			break;
		case MIC_CLOSE:
			that->_close();
			return 1;
		}
	}
	return 0;
}

void mic_writer_winmm_impl::_open() {
	 MMRESULT result;

	// open device
	WAVEFORMATEX fmt;
	fmt.wFormatTag = WAVE_FORMAT_PCM;     // simple, uncompressed format
	fmt.nChannels = 1;                    //  1=mono, 2=stereo
	fmt.nSamplesPerSec = _freq;           // 44100
	fmt.nAvgBytesPerSec = _freq * 2;      // = nSamplesPerSec * n.Channels * wBitsPerSample/8
	fmt.nBlockAlign = _elemsize;          // = n.Channels * wBitsPerSample/8
	fmt.wBitsPerSample = _elemsize * 8;   //  16 for high quality, 8 for telephone-grade
	fmt.cbSize = 0;
	result = waveInOpen(&hWaveIn, WAVE_MAPPER, &fmt, threadId, 0L, WAVE_FORMAT_DIRECT | CALLBACK_THREAD);
	if(result != MMSYSERR_NOERROR) {
		throw "Could not open wave input";
	}

	// init and add buffers
	for(int i = 0; i < NUM_BUFFERS; i++) {
		_init(headers[i]);
		waveInAddBuffer(hWaveIn, &headers[i], sizeof(WAVEHDR));
	}

#ifdef _DEBUG
	MMIOINFO info1;
	MMCKINFO info2, info3;
	ZeroMemory(&info1, sizeof(MMIOINFO));
	ZeroMemory(&info2, sizeof(MMCKINFO));
	ZeroMemory(&info3, sizeof(MMCKINFO));
	_debug = mmioOpen("debug.wav", &info1, MMIO_CREATE | MMIO_WRITE);
	info2.fccType = mmioFOURCC('W', 'A', 'V', 'E');
	mmioCreateChunk(_debug, &info2, MMIO_CREATERIFF);
#endif
}

void mic_writer_winmm_impl::_init(WAVEHDR& hdr) {
	hdr.lpData = (LPSTR) new byte[_bufsize];
	hdr.dwBufferLength = _bufsize;
	hdr.dwBytesRecorded = 0;
	hdr.dwUser = 0;
	hdr.dwFlags = 0;
	hdr.dwLoops = 0;
	waveInPrepareHeader(hWaveIn, &hdr, sizeof(WAVEHDR));
}

void mic_writer_winmm_impl::_start() {
	waveInStart(hWaveIn);
}

void mic_writer_winmm_impl::_process(WAVEHDR *hdr) {
	// write to output
	_output.write((byte *)hdr->lpData, hdr->dwBytesRecorded);
	// add it back to wave input again
	waveInAddBuffer(hWaveIn, hdr, sizeof(WAVEHDR));
#ifdef _DEBUG
	mmioWrite(_debug, hdr->lpData, hdr->dwBytesRecorded);
#endif
}

void mic_writer_winmm_impl::_stop() {
	waveInStop(hWaveIn);
}

void mic_writer_winmm_impl::_deinit(WAVEHDR& hdr) {
	byte *buf = reinterpret_cast<byte *>(hdr.lpData);
	waveInUnprepareHeader(hWaveIn, &hdr, sizeof(WAVEHDR));
	delete [] buf;
}

void mic_writer_winmm_impl::_close() {
	// set closed flag
	_closed = true;

	// force stop
	_stop();

	// free buffers
	for(int i = 0; i < NUM_BUFFERS; i++) {
		_deinit(headers[i]);
	}

	// close device
	waveInClose(hWaveIn);

#ifdef _DEBUG
	mmioClose(_debug, NULL);
#endif
}

mic_writer_impl *mic_writer_impl::create(ostream<byte>& output, int elemsize, int sample_rate) {
	return new mic_writer_winmm_impl(output, elemsize, sample_rate);
}

void mic_writer_impl::destroy(mic_writer_impl *impl) {
	delete impl;
}

