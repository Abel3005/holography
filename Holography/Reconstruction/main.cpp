#include <iostream>
#include <opencv2/opencv.hpp>
#include "fresnelDiffraction.h"
using namespace cv;
int main()
{
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
#pragma omp parallel for num_threads(2)
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
	std::cout << "time : " << end - start << std::endl;
	
	//propAssimple(objectWave_pc, Nx, lambda, z_obj, p);
	intensity(objectWave_pc, Nx);
	std::cout << objectWave_pc << std::endl;
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
}