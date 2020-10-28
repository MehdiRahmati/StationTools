#include <stdio.h>
#include <time.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
//#include <complex.h>
#include <fftw3.h>
//#include <complex.h>
#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif
#include "mysac.h"

#define PI 3.14159265358979311600
#define NPTS_SAC_FILE_MAX 17280000
#define N_MAX_ZEROES 25
#define N_MAX_POLES 25
#define ZERO_MACHINE_VALUE 0.00000001

typedef struct {
        int N;
        double *re;
        double *im;
} SPECTRUM;


typedef struct {
	double real;
	double imag;
} COMPLEX_RP;


typedef struct {
	int n_zeroes;
	int n_poles;
	double scale;
	COMPLEX_RP poles[N_MAX_POLES];
	COMPLEX_RP zeroes[N_MAX_ZEROES];
} POLE_ZERO;

void linear_interpolate_psd(float *freqs, int nfreqs, float *input, float *out, float minPer, float maxPer, float incPer);
int compute_nodes(float max, float min, float inc);
COMPLEX_RP generate_response(POLE_ZERO *pz, double freq);
void decon_response_function(fftw_complex *data_spect, fftw_complex *response_spect, fftw_complex *spect_out, int nspec);
COMPLEX_RP average_spectral_ratio(fftw_complex *signal_spectrum,fftw_complex *untapered_spectrum, int nspec);
//void fftw_forward(double *input, SPECTRUM *spectrum, int npts);
//void fftw_backward(SPECTRUM *spectrum, double *array_out);
void cos_taper(double *arr_in, double *arr_out, int npts);
void remove_mean(double *arr_in, double *arr_out, int npts);
void remove_trend(double *arr_in, double *arr_out, int npts);
COMPLEX_RP conjugate(COMPLEX_RP *comp);
COMPLEX_RP complex_multiply(COMPLEX_RP *comp1, COMPLEX_RP *comp2);
COMPLEX_RP complex_divide(COMPLEX_RP *comp1, COMPLEX_RP *comp2);
COMPLEX_RP complex_add(COMPLEX_RP *comp1, COMPLEX_RP *comp2);
COMPLEX_RP complex_subtract(COMPLEX_RP *comp1, COMPLEX_RP *comp2);
COMPLEX_RP complex_exp(COMPLEX_RP *comp);
double complex_amplitude(COMPLEX_RP *comp);
void initialize_spectrum(SPECTRUM *spect, int npts);
void write_sac (char *fname, float *sig, SAC_HD *SHD);
int read_sac (char *fname, float *sig, SAC_HD *SHD, int nmax);
int read_sac_pole_zero_file(char *filename, POLE_ZERO *pz_out);
void print_sac_pole_zero_file(POLE_ZERO *pz);
void zero_pole_zero(POLE_ZERO *pz);
void get_amplitude_spectrum(SPECTRUM *spect_in, double *arr_out);
void time_derivative(double *in, double *out, double delta, int npts);
void bw6plp(double *x, double *arr_out, int np, double dt, double lpfc);
void decimate_signal(double *input_signal, double *arr_out, int npts, double dt_in, double dt_out);
int gcd(int a, int b);
double find_common_sample_rate(double dt_1, double dt_2);
double find_time_shift(int sec_1, int msec_1, int sec_2, int msec_2);
void phase_shift(double t_shift, double dt,int np, fftw_complex *spect_in);
void signal_to_ground_acceleration(double *in, double *out, int npts, double delta, POLE_ZERO *pz, int acc_flag,
								   //the following are for memory usage. By passing the arrays in, we avoid re-allocating and freeing memory
								   double *signal, fftw_complex *signal_spectrum, fftw_complex *response_spectrum,
								   double *fftw_dbl, fftw_complex *fftw_cmplx,
								   fftw_plan plan_forward, fftw_plan plan_backward);
double recursive_mean(double new_value, double old_mean, int n);
double recursive_standard_deviation(double new_value, double current_mean, double prev_value, int n);
void running_mean_psd(double *arr_in, double *arr_out, int npts);
void calculate_psd(fftw_complex *spect, double *psd_out, double dt, int npts);
void smooth_1d_array(double *arr_in, double *arr_out, int npts, int n_smooth);
void hann_taper(double *arr_in, double *arr_out, int npts);
void hamming_taper(double *arr_in, double *arr_out, int npts);
void get_phase_spectrum(SPECTRUM *spect_in, double *arr_out);
int smooth_2d_array(double **dat, double **out, int n_x, int n_y, int n_smooth);
int compute_log10_nodes(double max, double min);
void log10_linear_interpolate(double *freqs, int nfreqs, double *input, double *freqsOut, double *out, int *nFreqsOut, double minPer, double maxPer);
void calculate_noise(double *out, fftw_complex *input, fftw_complex *cross, float dt, int npts);
void finish_coherence(fftw_complex *numerMean, fftw_complex *numerSD, fftw_complex *denom, int npts, double *meanOut, double *sdOut);
void compute_coherence(fftw_complex *spect1, fftw_complex *spect2, int npts, int j, fftw_complex *numerMean, fftw_complex *numerSD, fftw_complex *denom);
void finish_incoherence(fftw_complex *numerMean, fftw_complex *numerSD, fftw_complex *denom, int npts, float dt, double *meanOut1, double *sdOut1, double *meanOut2, double *sdOut2);
void finish_two_psds(fftw_complex *amps, int npts, float dt, double *out1, double *out2);
void compute_tri_noise(fftw_complex *spect1, fftw_complex *spect2, fftw_complex *spect3, int npts, int j,
					   fftw_complex *P11, fftw_complex *P12, fftw_complex *P13,
					   fftw_complex *P21, fftw_complex *P22, fftw_complex *P23,
					   fftw_complex *P31, fftw_complex *P32, fftw_complex *P33);

void finish_tri_noise(fftw_complex *P11, fftw_complex *P12, fftw_complex *P13,
					  fftw_complex *P21, fftw_complex *P22, fftw_complex *P23,
					  fftw_complex *P31, fftw_complex *P32, fftw_complex *P33,
					  int npts, float dt, fftw_complex *N11, fftw_complex *N22, fftw_complex *N33);



