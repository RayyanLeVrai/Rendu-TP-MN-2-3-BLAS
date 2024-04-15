#include "mnblas.h"
#include "complexe2.h"
#include <omp.h>
#include <stdio.h>
#include <math.h>

float mnblas_sasum(const int N, const float *X, const int incX) {
    float sum = 0.0;
    #pragma omp parallel for reduction(+:sum)
    for (int i = 0; i < N; i += incX) {
        sum += fabs(X[i]);
    }
    return sum;
}

double mnblas_dasum(const int N, const double *X, const int incX) {
    double sum = 0.0;
    #pragma omp parallel for reduction(+:sum)
    for (int i = 0; i < N; i += incX) {
        sum += fabs(X[i]);
    }
    return sum;
}

float mnblas_scasum(const int N, const void *X, const int incX) {
    float sum = 0.0;
    const complexe_float_t *cX = (const complexe_float_t *) X;
    #pragma omp parallel for reduction(+:sum)
    for (int i = 0; i < N; i += incX) {
        sum += sqrt(cX[i].real * cX[i].real + cX[i].imaginary * cX[i].imaginary);
    }
    return sum;
}

double mnblas_dzasum(const int N, const void *X, const int incX) {
    double sum = 0.0;
    const complexe_double_t *zX = (const complexe_double_t *) X;
    #pragma omp parallel for reduction(+:sum)
    for (int i = 0; i < N; i += incX) {
        sum += sqrt(zX[i].real * zX[i].real + zX[i].imaginary * zX[i].imaginary);
    }
    return sum;
}
