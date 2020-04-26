#include <iostream>
#include <fftw3.h>

#include "fresnelDiffraction.h"


void intensity(fftw_complex* u, int N) {
	for (int i = 0; i < N * N; i++) {
		double re = u[i][0];
		double im = u[i][1];
		u[i][0] = re * re + im * im;
	}
}
void response(fftw_complex* h, int N, double lamda, double z, double p) {
	fftw_complex tmp;
	tmp[0] = 0.0;
	tmp[1] = 0.0;
	for (double n = 0; n < N; n++) {
		for (double m = 0; m < N; m++) {
			int idx = m + n * N;
			double dfx = (m - N / 2) * p;
			double dfy = (n - N / 2) * p;
			double phase = (dfx * dfx + dfy * dfy) * M_PI / (lamda * z);
			h[idx][0] = cos(phase);
			h[idx][1] = sin(phase);

		}
	}
}
void fft(fftw_complex* u1, fftw_complex* u2, int N) {
	fftw_plan plan = fftw_plan_dft_2d(
		N, N, u1, u2, FFTW_FORWARD, FFTW_ESTIMATE
	);
	fftw_execute(plan);
	fftw_destroy_plan(plan);
}
void ifft(fftw_complex* u1, fftw_complex* u2, int N) {
	fftw_plan plan = (fftw_plan)fftw_plan_dft_2d(
		N, N, u1, u2, FFTW_BACKWARD, FFTW_ESTIMATE
	);
	fftw_execute(plan);
	fftw_destroy_plan(plan);
}
void fft_shift(fftw_complex* u, int N) {
	int Nh = N / 2;
	for (int i = 0; i < Nh; i++) {
		for (int j = 0; j < Nh; j++) {

			fftw_complex tmp1, tmp2;
			int adr1, adr2;

			adr1 = j + i * N;
			adr2 = (j + Nh) + (i + Nh) * N;
			tmp1[0] = u[adr1][0]; tmp1[1] = u[adr1][1];
			tmp2[0] = u[adr2][0]; tmp2[1] = u[adr2][1];
			u[adr1][0] = tmp2[0]; u[adr1][1] = tmp2[1];
			u[adr2][0] = tmp1[0]; u[adr2][1] = tmp1[1];

			adr1 = (j + Nh) + i * N;
			adr2 = (j + Nh) + (i + Nh) * N;
			tmp1[0] = u[adr1][0]; tmp1[1] = u[adr1][1];
			tmp2[0] = u[adr2][0]; tmp2[1] = u[adr2][1];

			u[adr1][0] = tmp2[0]; u[adr1][1] = tmp2[1];
			u[adr2][0] = tmp1[0]; u[adr2][1] = tmp1[1];
		}
	}
}
void ifft_shift(fftw_complex* u, int N) {
	int Nh = (N + 1) / 2;
	for (int i = 0; i < Nh; i++) {
		for (int j = 0; j < Nh; j++) {

			fftw_complex tmp1, tmp2;
			int adr1, adr2;

			adr1 = j + i * N;
			adr2 = (j + Nh) + (i + Nh) * N;
			tmp1[0] = u[adr1][0]; tmp1[1] = u[adr1][1];
			tmp2[0] = u[adr2][0]; tmp2[1] = u[adr2][1];
			u[adr1][0] = tmp2[0]; u[adr1][1] = tmp2[1];
			u[adr2][0] = tmp1[0]; u[adr2][1] = tmp1[1];

			adr1 = (j + Nh) + i * N;
			adr2 = (j + Nh) + (i + Nh) * N;
			tmp1[0] = u[adr1][0]; tmp1[1] = u[adr1][1];
			tmp2[0] = u[adr2][0]; tmp2[1] = u[adr2][1];

			u[adr1][0] = tmp2[0]; u[adr1][1] = tmp2[1];
			u[adr2][0] = tmp1[0]; u[adr2][1] = tmp1[1];
		}
	}
}
void mul_complex(fftw_complex* a, fftw_complex* b, fftw_complex* c, int N) {
	fftw_complex tmp;
	for (int i = 0; i < N * N; i++) {
		tmp[0] = a[i][0] * b[i][0] - a[i][1] * b[i][1];
		tmp[1] = a[i][0] * b[i][0] + a[i][1] * b[i][0];
		c[i][0] = tmp[0];
		c[i][1] = tmp[1];

	}
}
void mul_dbl(fftw_complex* u, double a, int N) {
	for (int i = 0; i < N * N; i++) {
		u[i][0] *= a;
		u[i][1] *= a;
	}
}
void fresnel_fft(fftw_complex* u1, int N, double lamda, double z, double p)
{
	fft_shift(u1, N);
	fft(u1, u1, N);

	fftw_complex* u2 = (fftw_complex*)malloc(sizeof(fftw_complex) * N * N);
	response(u2, N, lamda, z, p);
	fft_shift(u2, N);
	fft(u2, u2, N);

	mul_complex(u1, u2, u1, N);
	ifft(u1, u1, N);
	fft_shift(u1, N);

	mul_dbl(u1, 1.0 / (N * N), N);
	std::free(u2);
}
void propAssimple(fftw_complex* u1, int N, double lamda, double z, double p)
{
	fftw_complex* h = (fftw_complex*)malloc(sizeof(fftw_complex) * N * N);
	fft_shift(u1, N);
	fft(u1, u1, N);
	fft_shift(u1, N);
	response(h, N, lamda, z, p);
	mul_complex(u1, h, u1, N);
	ifft_shift(u1, N);
}