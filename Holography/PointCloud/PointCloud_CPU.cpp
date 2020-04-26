#include <iostream>
#include "Pointcloud_init.h"
#include <cmath>
#include <time.h>


#define M_PI 3.14159265358979
typedef double fftw_complex[2];
using namespace std;
void mul_complex(fftw_complex* a, fftw_complex* b, fftw_complex* c, int N);
int main() {
    //%%%%%%%%%%%%%% Read Point Cloud Data%%%%%%%%%%%%%%%%%%%%%%%%%
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

        /*% exp(1i * k0 * r) for forward prop.i.e.z < 0
        % exp(-1i * k0 * r) for back    prop.i.e.z > 0
        %
        %obj.wave is back propagated when z > 0
        %*/

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
    cout << "time : " << end - start << endl;
    free(x);
    free(y);
    free(z);
    getchar();

    return 0;
}