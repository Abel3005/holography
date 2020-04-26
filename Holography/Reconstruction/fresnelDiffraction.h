
#pragma once
#define M_PI 3.14
#pragma once
#define lambda 660e-9
#define Nx  10000
#define Ny  10000
// Pixel Pitch
#define dx  2e-6
#define dy  2e-6
// Segment
#define NumSeg_x  2
#define NumSeg_y  2
// Off - axis angle
#define off_axis_x  1
#define off_axis_y  1
// Reconstruction Distance from SLM
#define z_obj  0.1 // 1. to Object Plane
#define z_obs  0.2 // 2. to Observation Palne


typedef double fftw_complex[2];
void fresnel_fft(fftw_complex *u1, int N, double lamda, double z, double p);
void response(fftw_complex *h, int N, double lamda, double z, double p);
fftw_complex* fresnel_direct(fftw_complex *u, int N, double lamda, double z, double p);
void intensity(fftw_complex *u, int N);
void fft(fftw_complex *u1, fftw_complex *u2, int N);
void ifft(fftw_complex *u1, fftw_complex *u2, int N);
void fft_shift(fftw_complex *u, int N);
void mul_complex(fftw_complex *a, fftw_complex *b, fftw_complex *c, int N);
void mul_dbl(fftw_complex *u, double a, int N);
void fresnel_fft(fftw_complex *u1, int N, double lamda, double z, double p);
void propAssimple(fftw_complex* u1, int N, double lamda, double z, double p);
