#include <iostream>
#include <fftw3.h>
#include <opencv2/opencv.hpp>
#include "PointCloud_init.h"
#include "FresnelDiffraction.h"
using namespace cv;
using namespace std;
#define H_ROW 1024
#define H_COL 1024
int main()
{
	double p = 1 / (H_ROW * dx);
	/*Mat img = imread("lena.bmp", IMREAD_GRAYSCALE);
	Mat image1(H_COL, H_ROW, CV_8UC1);
	cv::imshow("img",img);

	int len = img.rows * img.cols;
	std::cout << len << endl;
	//cout << len << endl;
	fftw_complex* imgsrc1 = (fftw_complex*)malloc((sizeof(fftw_complex)) * 512 *512);
	fftw_complex* imgsrc = (fftw_complex*)malloc((sizeof(fftw_complex))*H_ROW*H_COL);
	fftw_complex* u_lens = (fftw_complex*)malloc((sizeof(fftw_complex)) *H_ROW * H_COL);
	uchar* temp = img.data;
	int tempInt = 0;


	for (int i = tempInt; i < H_ROW*H_COL; i++) {
		imgsrc[i][0] = 0;
		imgsrc[i][1] = 0;
	}
	for (int i = 0; i < img.cols; i++) {
		for (int j = 0; j < img.rows; j++) {
			int matLength = (i * img.rows) + j;

			imgsrc1[matLength][0] = img.at<uchar>(i, j);
		}
	}
	fresnel_fft(imgsrc1, 512, 633e-9, 10e-2, 10e-6);
	response(u_lens, 512, 633e-9, -5e-2, 10e-6);
	mul_complex(imgsrc1, u_lens, imgsrc1, 512);
	fresnel_fft(imgsrc1, 512, 633e-9, 10e-2, 10e-6);
	intensity(imgsrc1, 512);

	for (int i = 0; i < img.cols; i++) {
		for (int j = 0; j < img.rows; j++) {
			int matLength = (i * img.rows) + j;

			img.at<uchar>(i, j) = imgsrc1[matLength][0];
		}
	}
	imshow("b", img);
	printf("%d", img.rows);
	for (int i = 0; i < img.cols; i++) {
		for (int j = 0; j < img.rows; j++) {
			int abc = len;
			int matLength = (i+256) * H_ROW + (j+256);
			imgsrc[matLength][0] = img.at<uchar>(i, j);
			imgsrc[matLength][1] = 0;
			tempInt++;
		}
	}
	fftw_complex* median2 = (fftw_complex*)malloc(sizeof(fftw_complex) * H_COL * H_ROW);
	memcpy(median2, imgsrc, sizeof(fftw_complex) * H_COL * H_ROW);
	propAssimple(median2, H_ROW, lambda, z_obj, p);
	intensity(median2, H_ROW);

	for (int i = 0; i < image1.cols; i++) {
		for (int j = 0; j < image1.rows; j++) {
			int matLength = (i * image1.rows) + j;

			image1.at<uchar>(i, j) = median2[matLength][0];
		}
	}
	imshow("a", image1);*/
	FILE* infp;
	int npoint = 0;
	fopen_s(&infp, "Input/tetrahedron_01K.txt", "rb");
	fscanf_s(infp, "%d", &npoint);
	float* x = (float*)malloc(sizeof(float) * npoint);
	float* y = (float*)malloc(sizeof(float) * npoint);
	float* z = (float*)malloc(sizeof(float) * npoint);
	clock_t start, end;

	for (int i = 0; i < npoint; i++) {
		fscanf_s(infp, "%f", &(x[i]));
		fscanf_s(infp, "%f", &(y[i]));
		fscanf_s(infp, "%f", &(z[i]));
	}
	fclose(infp);
	int NxSeg = Nx / NumSeg_x;
	int NySeg = Ny / NumSeg_y;
	double k0 = 2 * M_PI / lambda;
	double cutOff_deg_X = 2 * asin(lambda / (2 * dx));
	double cutOff_deg_Y = 2 * asin(lambda / (2 * dy));
	int ii = 1;

	fftw_complex* objectWave_pc = (fftw_complex*)malloc(sizeof(fftw_complex) * Nx * Ny);
	memset(objectWave_pc, 0, sizeof(fftw_complex) * Nx * Ny);
	start = clock();
	for (int s = 0; s < npoint; s++) {



		if (z[s] > 0)
			ii = -1;
		for (int i = 1; i <= Ny; i++) {
			for (int j = 1; j <= Nx; j++) {
				double hx = ((j - Nx) / 2) * dx;
				double hy = ((i - Ny) / 2) * dx;
				double hz = 0;
				double Rhj = sqrt((hx - x[s]) * (hx - x[s]) + (hy - y[s]) * (hy - y[s]) + (hz - z[s]) * (hz - z[s]));
				bool cutOff_X = abs((hx - x[s]) / (hz - z[s])) < tan(cutOff_deg_X / 2 * M_PI / 180);
				bool cutOff_Y = abs((hy - y[s]) / (hz - z[s])) < tan(cutOff_deg_Y / 2 * M_PI / 180);
				bool cutOff = cutOff_X && cutOff_Y;
				if (!cutOff) {
					objectWave_pc[(i - 1) * Nx + (j - 1)][0] += (hz - z[s]) * cos(ii * (k0 * Rhj)) / (Rhj * Rhj);
					objectWave_pc[(i - 1) * Ny + (j - 1)][1] += (hz - z[s]) * sin(ii * (k0 * Rhj)) / (Rhj * Rhj);
				}

			}
		}
	}
	end = clock();
	std::cout << "time : " << end - start << endl;
	propAssimple(objectWave_pc, Nx, lambda, z_obj, p);
	intensity(objectWave_pc, Nx);
	Mat image = Mat::zeros((Ny / 20), (Nx / 20), CV_8UC1);
	for (int i = 0; i < (Ny / 20); i++) {
		for (int j = 0; j < (Nx / 20); j++) {
			int matLength = i * (Nx * 20) + (j * 20);

			image.at<uchar>(j, i) = objectWave_pc[matLength][0];
		}
	}
	std::free(x);
	std::free(y);
	std::free(z);
	cv::imshow("re", image);
	getchar();

	//fftw_complex* imgsrc = (fftw_complex*)malloc((sizeof(fftw_complex)) * 512 * 512);
	//fftw_complex* u_lens = (fftw_complex*)malloc((sizeof(fftw_complex)) * 512 * 512);
	//uchar* temp = img.data;
	//int tempInt = 0;



	/*for (int i = 0; i < img.cols; i++) {
		for (int j = 0; j < img.rows; j++) {
			int abc = len;
			int matLength = i * img.rows + j;
			imgsrc[matLength][0] = img.at<uchar>(j, i);
			imgsrc[matLength][1] = 0;
			tempInt++;
		}
	}*/

	/*fresnel_fft(imgsrc,H_ROW,633e-9 , 10e-2, 10e-6);
	fftw_complex* median = (fftw_complex*)malloc(sizeof(fftw_complex) * H_COL*H_ROW);
	memcpy(median, imgsrc, sizeof(fftw_complex) * H_COL * H_ROW);
	intensity(median, H_ROW);

	for (int i = 0; i < image1.cols; i++) {
		for (int j = 0; j < image1.rows; j++) {
			int matLength = i * image1.rows + j;

			image1.at<uchar>(i, j) = median[matLength][0];
		}
	}
	for (int i = 0; i < img.cols; i++) {
		for (int j = 0; j < img.rows; j++) {
			int matLength = i * img.rows + j;

			img.at<uchar>(i, j) = imgsrc[matLength][0];
		}
	}

	cv::imshow("median", image1);
	response(u_lens, H_ROW, 633e-9, -5e-2, 10e-6);
	mul_complex(imgsrc, u_lens, imgsrc, H_ROW);

	fresnel_fft(imgsrc, H_ROW, 633e-9, 10e-2, 10e-6);
	intensity(imgsrc, H_ROW);
	for (int i = 0; i < image1.cols; i++) {
		for (int j = 0; j < image1.rows; j++) {
			int matLength = i * image1.rows + j;

			image1.at<uchar>(i, j) = imgsrc[matLength][0];
		}
	}
	cv::imshow("result",image1);
	cv::waitKey();*/
	return 0;
}
/*fftw_complex* fresnel_direct(fftw_complex *u, int N, double lamda, double z, double p) {
	fftw_complex *u2 = (fftw_complex*)malloc(sizeof(fftw_complex)*N*N);
	for (int n2 = 0; n2 < N; n2++) {
		for (int m2 = 0; m2 < N; m2++) {
			fftw_complex tmp;
			tmp[0] = 0.0;
			tmp[1] = 0.0;
			for (int n1 = 0; n1 < N; n1++) {
				for (int m1 = 0; m1 < N; m1++) {
					int idx1 = m1 + n1 * N;
					double dx = ((m2 - N / 2) - (m1 - N / 2))*p;
					double dy = ((n2 - N / 2) - (n1 - N / 2))*p;
					double phase = (dx*dx + dy * dy)*M_PI / (lamda*z);
					fftw_complex e, t;
					e[0] = cos(phase);
					e[1] = sin(phase);
					t[0] = u[idx1][0];
					t[1] = u[idx1][1];
					tmp[0] += t[0] * e[0] - t[1] * e[1];
					tmp[1] += t[0] * e[1] + t[1] * e[0];
				}
			}
			int idx2 = m2 + n2 * N;
			u2[idx2][0] = tmp[0];
			u2[idx2][1] = tmp[1];
		}

	}
	free(u);
	return u2;
}*/

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