




























































































































typedef double freq_t;
typedef double signal_t;
typedef double spectrum_t;
typedef unsigned char mask_t;












struct _spl_scale_generate_linear_args { int K ; freq_t F1 ; freq_t F2 ; freq_t *sc; };  void spl_scale_generate_linear (int K, freq_t F1, freq_t F2, freq_t *sc);  void spl_scale_generate_linear_str ( struct _spl_scale_generate_linear_args * );;
struct _spl_scale_generate_log_args { int K ; freq_t F1 ; freq_t F2 ; freq_t *sc; };  void spl_scale_generate_log (int K, freq_t F1, freq_t F2, freq_t *sc);  void spl_scale_generate_log_str ( struct _spl_scale_generate_log_args * );;
struct _spl_scale_generate_mel_args { int K ; freq_t F1 ; freq_t F2 ; freq_t *sc; };  void spl_scale_generate_mel (int K, freq_t F1, freq_t F2, freq_t *sc);  void spl_scale_generate_mel_str ( struct _spl_scale_generate_mel_args * );;
struct _spl_scale_generate_bark_args { int K ; freq_t F1 ; freq_t F2 ; freq_t *sc; };  void spl_scale_generate_bark (int K, freq_t F1, freq_t F2, freq_t *sc);  void spl_scale_generate_bark_str ( struct _spl_scale_generate_bark_args * );;
struct _spl_scale_generate_model_args { int K ; freq_t F1 ; freq_t F2 ; freq_t *sc; };  void spl_scale_generate_model (int K, freq_t F1, freq_t F2, freq_t *sc);  void spl_scale_generate_model_str ( struct _spl_scale_generate_model_args * );;


struct _spl_filter_memory_linear_args { int K ; int N ; freq_t F ; freq_t F1 ; freq_t F2 ; const signal_t *s ; spectrum_t *sp ; double ksi; };  unsigned spl_filter_memory_linear (int K, int N, freq_t F, freq_t F1, freq_t F2, const signal_t *s, spectrum_t *sp, double ksi);  unsigned spl_filter_memory_linear_str ( struct _spl_filter_memory_linear_args * );;
struct _spl_filter_memory_args { int K ; int N ; freq_t F ; const freq_t *sc ; const signal_t *s ; spectrum_t *sp ; double ksi; };  unsigned spl_filter_memory (int K, int N, freq_t F, const freq_t *sc, const signal_t *s, spectrum_t *sp, double ksi);  unsigned spl_filter_memory_str ( struct _spl_filter_memory_args * );;
struct _spl_filter_binary_file_args { int K ; freq_t F ; const freq_t *sc ; const char *in_file ; const char *out_file ; double ksi; };  unsigned spl_filter_binary_file (int K, freq_t F, const freq_t *sc, const char *in_file, const char *out_file, double ksi);  unsigned spl_filter_binary_file_str ( struct _spl_filter_binary_file_args * );;
struct _spl_filter_wave_file_args { int K ; const freq_t *sc ; const char *in_file ; const char *out_file ; double ksi; };  unsigned spl_filter_wave_file (int K, const freq_t *sc, const char *in_file, const char *out_file, double ksi);  unsigned spl_filter_wave_file_str ( struct _spl_filter_wave_file_args * );;

struct _spl_load_signal_args { char *file ; int N ; signal_t *buf ; freq_t *F; };  unsigned spl_load_signal (char *file, int N, signal_t *buf, freq_t *F);  unsigned spl_load_signal_str ( struct _spl_load_signal_args * );;

struct _spl_mask_function_generate_args { int K ; int Ws ; freq_t *sc ; double *H; };  unsigned spl_mask_function_generate (int K, int Ws, freq_t *sc, double *H);  unsigned spl_mask_function_generate_str ( struct _spl_mask_function_generate_args * );;

struct _spl_mask_memory_args { int K ; int N ; const freq_t *sc ; const spectrum_t *sp ; mask_t *m ; double ksi; };  unsigned spl_mask_memory (int K, int N, const freq_t *sc, const spectrum_t *sp, mask_t *m, double ksi);  unsigned spl_mask_memory_str ( struct _spl_mask_memory_args * );;
struct _spl_mask_binary_file_args { int K ; const freq_t *sc ; const char *in_file ; const char *out_file ; bool out_bits ; double ksi; };  unsigned spl_mask_binary_file (int K, const freq_t *sc, const char *in_file, const char *out_file, bool out_bits, double ksi);  unsigned spl_mask_binary_file_str ( struct _spl_mask_binary_file_args * );;
struct _spl_mask_memory_ext_args { int K ; int N ; const freq_t *sc ; const spectrum_t *sp ; mask_t *m ; double ksi ; double delta ; double rho ; bool border_effect ; bool allow_fast; };  unsigned spl_mask_memory_ext (int K, int N, const freq_t *sc, const spectrum_t *sp, mask_t *m, double ksi, double delta, double rho, bool border_effect, bool allow_fast);  unsigned spl_mask_memory_ext_str ( struct _spl_mask_memory_ext_args * );;
struct _spl_mask_binary_file_ext_args { int K ; const freq_t *sc ; const char *in_file ; const char *out_file ; bool out_bits ; double ksi ; double delta ; double rho ; bool border_effect ; bool allow_fast; };  unsigned spl_mask_binary_file_ext (int K, const freq_t *sc, const char *in_file, const char *out_file, bool out_bits, double ksi, double delta, double rho, bool border_effect, bool allow_fast);  unsigned spl_mask_binary_file_ext_str ( struct _spl_mask_binary_file_ext_args * );;

struct _spl_pitch_memory_args { int K ; int N ; const freq_t *sc ; const mask_t *m ; short *p ; double ksi; };  unsigned spl_pitch_memory (int K, int N, const freq_t *sc, const mask_t *m, short *p, double ksi);  unsigned spl_pitch_memory_str ( struct _spl_pitch_memory_args * );;
struct _spl_pitch_binary_file_args { int K ; int N ; const freq_t *sc ; const char *in_file ; const char *out_file ; bool in_bits ; double ksi; };  unsigned spl_pitch_binary_file (int K, int N, const freq_t *sc, const char *in_file, const char *out_file, bool in_bits, double ksi);  unsigned spl_pitch_binary_file_str ( struct _spl_pitch_binary_file_args * );;
struct _spl_pitch_memory_ext_args { int K ; int N ; const freq_t *sc ; const mask_t *m ; short *p ; double ksi ; double delta ; double rho ; bool border_effect ; double F1 ; double F2 ; int Nh; };  unsigned spl_pitch_memory_ext (int K, int N, const freq_t *sc, const mask_t *m, short *p, double ksi, double delta, double rho, bool border_effect, double F1, double F2, int Nh);  unsigned spl_pitch_memory_ext_str ( struct _spl_pitch_memory_ext_args * );;
struct _spl_pitch_binary_file_ext_args { int K ; int N ; const freq_t *sc ; const char *in_file ; const char *out_file ; bool in_bits ; double ksi ; double delta ; double rho ; bool border_effect ; double F1 ; double F2 ; int Nh; };  unsigned spl_pitch_binary_file_ext (int K, int N, const freq_t *sc, const char *in_file, const char *out_file, bool in_bits, double ksi, double delta, double rho, bool border_effect, double F1, double F2, int Nh);  unsigned spl_pitch_binary_file_ext_str ( struct _spl_pitch_binary_file_ext_args * );;






