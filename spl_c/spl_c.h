#ifndef _SPL_C_
#define _SPL_C_

#ifndef _SPL_COMMON_
#define EXPORT
#endif

#include "wrap_struct.h"

#ifndef _SPL_TYPES_
typedef double freq_t;
typedef double signal_t;
typedef double spectrum_t;
typedef unsigned char mask_t;
#else
using spl::freq_t;
using spl::signal_t;
using spl::spectrum_t;
using spl::mask_t;
#endif

#ifdef __cplusplus
extern "C" {
#endif

#define SCALE_ARGS (int K, freq_t F1, freq_t F2, freq_t *sc)
declare_wrap_func(void, spl_scale_generate_linear, 4, SCALE_ARGS);
declare_wrap_func(void, spl_scale_generate_log, 4, SCALE_ARGS);
declare_wrap_func(void, spl_scale_generate_mel, 4, SCALE_ARGS);
declare_wrap_func(void, spl_scale_generate_bark, 4, SCALE_ARGS);
declare_wrap_func(void, spl_scale_generate_model, 4, SCALE_ARGS);
#undef SCALE_ARGS

declare_wrap_func(unsigned, spl_filter_memory_linear, 8, (int K, int N, freq_t F, freq_t F1, freq_t F2, const signal_t *s, spectrum_t *sp, double ksi));
declare_wrap_func(unsigned, spl_filter_memory, 7, (int K, int N, freq_t F, const freq_t *sc, const signal_t *s, spectrum_t *sp, double ksi));
declare_wrap_func(unsigned, spl_filter_binary_file, 6, (int K, freq_t F, const freq_t *sc, const char *in_file, const char *out_file, double ksi));
declare_wrap_func(unsigned, spl_filter_wave_file, 5, (int K, const freq_t *sc, const char *in_file, const char *out_file, double ksi));

declare_wrap_func(unsigned, spl_load_signal, 4, (char *file, int N, signal_t *buf, freq_t *F));

declare_wrap_func(unsigned, spl_mask_function_generate, 4, (int K, int Ws, freq_t *sc, double *H));

declare_wrap_func(unsigned, spl_mask_memory, 6, (int K, int N, const freq_t *sc, const spectrum_t *sp, mask_t *m, double ksi));
declare_wrap_func(unsigned, spl_mask_binary_file, 6, (int K, const freq_t *sc, const char *in_file, const char *out_file, bool out_bits, double ksi));
declare_wrap_func(unsigned, spl_mask_memory_ext, 10, (int K, int N, const freq_t *sc, const spectrum_t *sp, mask_t *m, double ksi, double delta, double rho, bool border_effect, bool allow_fast));
declare_wrap_func(unsigned, spl_mask_binary_file_ext, 10, (int K, const freq_t *sc, const char *in_file, const char *out_file, bool out_bits, double ksi, double delta, double rho, bool border_effect, bool allow_fast));

declare_wrap_func(unsigned, spl_pitch_memory, 6, (int K, int N, const freq_t *sc, const mask_t *m, short *p, double ksi));
declare_wrap_func(unsigned, spl_pitch_binary_file, 7, (int K, int N, const freq_t *sc, const char *in_file, const char *out_file, bool in_bits, double ksi));
declare_wrap_func(unsigned, spl_pitch_memory_ext, 12, (int K, int N, const freq_t *sc, const mask_t *m, short *p, double ksi, double delta, double rho, bool border_effect, double F1, double F2, int Nh));
declare_wrap_func(unsigned, spl_pitch_binary_file_ext, 13, (int K, int N, const freq_t *sc, const char *in_file, const char *out_file, bool in_bits, double ksi, double delta, double rho, bool border_effect, double F1, double F2, int Nh));

#ifdef __cplusplus
}
#endif

#endif//_SPL_C_
