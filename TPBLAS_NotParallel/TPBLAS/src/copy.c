#include "mnblas.h"
#include "complexe2.h"
#include <omp.h>

void mncblas_scopy(const int N, const float *X, const int incX, float *Y, const int incY) {
    for (int i = 0; i < N; i += incX) {
        int j = (i / incX) * incY; // Calculer j en fonction de i
        if (j < N) {
            Y[j] = X[i];
        }
    }
}

void mncblas_dcopy(const int N, const double *X, const int incX, double *Y, const int incY) {
    for (int i = 0; i < N; i += incX) {
        int j = (i / incX) * incY; // Calculer j en fonction de i
        if (j < N) {
            Y[j] = X[i];
        }
    }
}

void mncblas_ccopy(const int N, const void *X, const int incX, void *Y, const int incY) {
    const complexe_float_t *complexX = (const complexe_float_t *)X;
    complexe_float_t *complexY = (complexe_float_t *)Y;

    for (int ix = 0; ix < N; ix += incX) {
        int iy = (ix / incX) * incY; // Calculer iy en fonction de ix
        if (iy < N) {
            complexY[iy].real = complexX[ix].real;
            complexY[iy].imaginary = complexX[ix].imaginary;
        }
    }
}

void mncblas_zcopy(const int N, const void *X, const int incX, void *Y, const int incY) {
    const complexe_double_t *complexX = (const complexe_double_t *)X;
    complexe_double_t *complexY = (complexe_double_t *)Y;

    for (int ix = 0; ix < N; ix += incX) {
        int iy = (ix / incX) * incY; // Calculer iy en fonction de ix
        if (iy < N) {
            complexY[iy].real = complexX[ix].real;
            complexY[iy].imaginary = complexX[ix].imaginary;
        }
    }
}
