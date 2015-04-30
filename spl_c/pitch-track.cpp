#define NOMINMAX
#include <windows.h>

#include "scale.h"
#include "spectrum.h"
#include "mask.h"
#include "vocal.h"

#include "membuf.h"
#include "iobuf.h"
#include "iomic.h"
#include "iofile.h"
#include "iowave.h"
#include "iosplit.h"
#include "resource.h"

using namespace spl;
using namespace io;

struct spec_params {
	freq_scale *sc;                   ///< Шкала частот
	io::istream<signal_t> *in_str;    ///< Входной поток (сигнал)
	io::ostream<spectrum_t> *out_str; ///< Выходной поток (спектр)
	freq_t F;                         ///< Частота дискретизации
	double ksi;                       ///< Величина, определяющая точность вычислений (размер окна)
};

DWORD WINAPI specThreadProc(LPVOID lpParam) {
	spec_params *params = (spec_params *) lpParam;
	return filter_stream(*params->sc, *params->in_str, *params->out_str, params->F, params->ksi);
}

struct mask_params {
	freq_scale *sc;                  ///< Шкала резонансных частот
	io::istream<spectrum_t> *in_str; ///< Входной поток  (спектр)
	io::ostream<mask_t> *out_str;    ///< Выходной поток (маска)
	mask_parameters p;               ///< Параметры маскировки
};

DWORD WINAPI maskThreadProc(LPVOID lpParam) {
	mask_params *params = (mask_params *) lpParam;
	return mask_stream(*params->sc, *params->in_str, *params->out_str, params->p);
}

struct pitch_params {
	freq_scale *sc;              ///< шкала частот
	io::istream<mask_t> *in_str; ///< входной поток (маска)
	io::ostream<short> *out_str; ///< выходной поток (номера каналов ЧОТ)
	mask_parameters *pm;         ///< параметры генерации маскирующей функции
	pitch_parameters p;          ///< параметры генерации шаблона маски
};

DWORD WINAPI pitchThreadProc(LPVOID lpParam) {
	pitch_params *params = (pitch_params *) lpParam;
	return pitch_track(*params->sc, *params->in_str, *params->out_str, *params->pm, params->p);
}

struct vocal_params {
	io::istream<short> *in_str;
	io::ostream<short> *out_str;
	freq_t F;
	vocal_parameters p;
};

DWORD WINAPI vocalThreadProc(LPVOID lpParam) {
	vocal_params *params = (vocal_params *) lpParam;
	return vocal_segment(*params->in_str, *params->out_str, params->F, params->p);
}

struct track_params {
	freq_scale *sc;             // шкала частот
	mic_writer<short> *imic;    // поток микрофона для управления
	io::istream<short> *in_str; // входной поток (номера каналов ЧОТ)
	HWND hWnd;
	bool started;
};

#define MSG_TRACK (WM_USER + 1)

DWORD WINAPI trackThreadProc(LPVOID lpParam) {
	track_params *params = (track_params *) lpParam;
	freq_scale& scale = *params->sc;
	io::istream<short>& input = *params->in_str;
	short prev_chan = -1, chan;
	while(input.get(chan)) {
		if(prev_chan != chan) {
			prev_chan = chan;
			HWND hWnd = params->hWnd;
			if(hWnd) {
				PostMessage(hWnd, MSG_TRACK, 0, chan);
			}
		}
	}
	return true;
}

#define GWL_USERDATA (-21)

INT_PTR CALLBACK trackDialogProc(HWND hWnd, UINT msg, WPARAM wParam, LPARAM lParam) {
	switch(msg) {
	case WM_INITDIALOG:
		ShowWindow(hWnd, SW_SHOW);
		SetWindowLong(hWnd, GWL_USERDATA, lParam);
		((track_params *)lParam)->hWnd = hWnd;
		return TRUE;
	case WM_COMMAND:
		switch(LOWORD(wParam)) {
		case IDC_BTN_START_STOP:
			if(HIWORD(wParam) == BN_CLICKED) {
				track_params *track_p = (track_params *) GetWindowLong(hWnd, GWL_USERDATA);
				if(track_p->started) {
					track_p->imic->stop();
					track_p->started = false;
					SetDlgItemText(hWnd, IDC_BTN_START_STOP, TEXT("Start"));
				} else {
					track_p->imic->start();
					track_p->started = true;
					SetDlgItemText(hWnd, IDC_BTN_START_STOP, TEXT("Stop"));
				}
			}
			return TRUE;
		case IDCANCEL:
			DestroyWindow(hWnd);
			return TRUE;
		}
		return FALSE;
	case MSG_TRACK:
		{
			if(lParam == -1) {
				SetDlgItemText(hWnd, IDC_TXT_PITCH, TEXT("-"));
			} else {
				TCHAR pitch_str[20];
				track_params *track_p = (track_params *) GetWindowLong(hWnd, GWL_USERDATA);
				swprintf(pitch_str, 20, TEXT("%.0lf"), track_p->sc->Fr[lParam]);
				SetWindowText(GetDlgItem(hWnd, IDC_TXT_PITCH), pitch_str);
			}
		}
		return TRUE;
	}
	return FALSE;
}

HANDLE createSimpleThread(LPTHREAD_START_ROUTINE startAddress, LPVOID lpParam) {
	return CreateThread(NULL, 10000, startAddress, lpParam, NULL, NULL);
}

int WINAPI WinMain(HINSTANCE hInstance, HINSTANCE hPrevInstance, LPSTR lpCmdLine, int nShowCmd)
{
	// constants
	const int num_channels = 256;
	const int sample_rate = 16000;
	const double spec_accuracy = 0.001;
	const double mask_accuracy = 0.001;

	// generate scale
	membuf<freq_t> scale_mem(num_channels);
	freq_scale scale;
	scale.K = num_channels;
	scale.Fr = scale_mem;
	scale_parameters scale_p;
	scale_p.a = 0;
	scale_p.b = num_channels - 1;
	scale_p.i = 0;
	scale_p.j = num_channels - 1;
	scale_p.Fi = 50;
	scale_p.Fj = 2000;
	scale_generate_model(scale, scale_p);

	// initiate pipe:
	iobuf<signal_t> sig_buf = iobuf<signal_t>();
	iobuf<spectrum_t> spec_buf = iobuf<spectrum_t>();
	iobuf<mask_t> mask_buf = iobuf<mask_t>();
	iobuf<short> pitch_buf = iobuf<short>();
	iobuf<short> vocal_buf = iobuf<short>();

	// files for debugging purposes:
	ofstream<short> mic_debug("microphone-debug.bin");
	ofstream<signal_t> sig_debug("signal-debug.bin");
	ofstream<spectrum_t> spec_debug("spectrum-debug.bin");
	ofstream<mask_t> mask_debug("mask-debug.bin");
	ofstream<short> pitch_debug("pitch-debug.bin");
	ofstream<short> vocal_debug("vocal-debug.bin");

	// tee
	tee<signal_t> sig_out(sig_buf.output(), sig_debug);
	tee<signal_t> spec_out(spec_buf.output(), spec_debug);
	tee<mask_t> mask_out(mask_buf.output(), mask_debug);
	tee<short> pitch_out(pitch_buf.output(), pitch_debug);
	tee<short> vocal_out(vocal_buf.output(), vocal_debug);

	// microphone input conversion:
	owrapelem<short, signal_t> mic_conv = owrapelem<short, signal_t>(sig_out, convert_sample_to_float);

	tee<short> mic_out(mic_conv, mic_debug);

	// input from microphone:
	mic_writer<short> imic = mic_writer<short>(mic_out, sample_rate);

	// spectrum thread:
	spec_params spec_p;
	spec_p.sc = &scale;
	spec_p.in_str = &sig_buf.input();
	spec_p.out_str = &spec_out;
	spec_p.F = sample_rate;
	spec_p.ksi = spec_accuracy;
	HANDLE specThread = createSimpleThread(specThreadProc, &spec_p);

	// mask thread:
	mask_params mask_p;
	mask_p.sc = &scale;
	mask_p.in_str = &spec_buf.input();
	mask_p.out_str = &mask_out;
	mask_p.p = DEFAULT_MASK_PARAMETERS;
	HANDLE maskThread = createSimpleThread(maskThreadProc, &mask_p);

	// pitch thread:
	pitch_params pitch_p;
	pitch_p.sc = &scale;
	pitch_p.in_str = &mask_buf.input();
	pitch_p.out_str = &pitch_out;
	pitch_p.pm = &mask_p.p;
	pitch_p.p.F1 = 70;
	pitch_p.p.F2 = 400;
	pitch_p.p.Nh = 2;
	HANDLE pitchThread = createSimpleThread(pitchThreadProc, &pitch_p);

	// vocal thread:
	vocal_params vocal_p;
	vocal_p.in_str = &pitch_buf.input();
	vocal_p.out_str = &vocal_out;
	vocal_p.F = sample_rate;
	vocal_p.p = DEFAULT_VOCAL_PARAMETERS;
//	HANDLE vocalThread = createSimpleThread(vocalThreadProc, &vocal_p);

	// track thread:
	track_params track_p;
	track_p.sc = &scale;
	track_p.imic = &imic;
	track_p.in_str = &pitch_buf.input();
	track_p.started = false;
	track_p.hWnd = NULL;
	HANDLE trackThread = createSimpleThread(trackThreadProc, &track_p);

	// track window:
	DialogBoxParam(hInstance, (LPCWSTR) IDD_TRACK_DIALOG, NULL, trackDialogProc, (LPARAM) &track_p);

	// destroy everything:
	CloseHandle(trackThread);
//	CloseHandle(vocalThread);
	CloseHandle(pitchThread);
	CloseHandle(maskThread);
	CloseHandle(specThread);

	return 0;
}
