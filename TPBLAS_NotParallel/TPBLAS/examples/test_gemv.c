#include "mnblas.h"
#include "complexe2.h"
#include "flop.h"
#include <stdio.h>
#include <omp.h>

#define M 128
#define N 128

typedef float vfloat[N];
typedef double vdouble[N];
typedef complexe_float_t vcomplex[N];
typedef complexe_double_t vzcomplex[N];

float A_f[M * N], X_f[N], Y_f[M];
double A_d[M * N], X_d[N], Y_d[M];
complexe_float_t A_c[M * N], X_c[N], Y_c[M];
complexe_double_t A_z[M * N], X_z[N], Y_z[M];

void init_vector_float(float *V, int size, float value) {
    for (int i = 0; i < size; i++) {
        V[i] = value;
    }
}

void init_vector_double(double *V, int size, double value) {
    for (int i = 0; i < size; i++) {
        V[i] = value;
    }
}

void init_vector_complex_float(complexe_float_t *V, int size, float real, float imag) {
    for (int i = 0; i < size; i++) {
        V[i].real = real;
        V[i].imaginary = imag;
    }
}

void init_vector_complex_double(complexe_double_t *V, int size, double real, double imag) {
    for (int i = 0; i < size; i++) {
        V[i].real = real;
        V[i].imaginary = imag;
    }
}

int main() {
    struct timespec start, end;
    init_nano();

    
    init_vector_float(A_f, M * N, 1.0f);
    init_vector_float(X_f, N, 2.0f);
    init_vector_float(Y_f, M, 0.0f);

    init_vector_double(A_d, M * N, 1.0);
    init_vector_double(X_d, N, 2.0);
    init_vector_double(Y_d, M, 0.0);

    init_vector_complex_float(A_c, M * N, 1.0, 0.0);
    init_vector_complex_float(X_c, N, 2.0, 0.0);
    init_vector_complex_float(Y_c, M, 0.0, 0.0);

    init_vector_complex_double(A_z, M * N, 1.0, 0.0);
    init_vector_complex_double(X_z, N, 2.0, 0.0);
    init_vector_complex_double(Y_z, M, 0.0, 0.0);

    // Tests mncblas_sgemv
    TOP_NANO(start);
    mncblas_sgemv(1, 1, M, N, 1.0f, A_f, N, X_f, 1, 0.0f, Y_f, 1);
    TOP_NANO(end);
    printf("sgemv float time %e seconds\n", diff_nano(&start, &end));

    // Tests mncblas_dgemv
    TOP_NANO(start);
    mncblas_dgemv(1, 1, M, N, 1.0, A_d, N, X_d, 1, 0.0, Y_d, 1);
    TOP_NANO(end);
    printf("sgemv double time %e seconds\n", diff_nano(&start, &end));

    // Tests mncblas_cgemv
    TOP_NANO(start);
    mncblas_cgemv(1, 1, M, N, &A_c[0], A_c, N, X_c, 1, &Y_c[0], Y_c, 1);
    TOP_NANO(end);
    printf("sgemv complex float time %e seconds\n", diff_nano(&start, &end));

    // Tests mncblas_zgemv
    TOP_NANO(start);
    mncblas_zgemv(1, 1, M, N, &A_z[0], A_z, N, X_z, 1, &Y_z[0], Y_z, 1);
    TOP_NANO(end);
    printf("sgemv complex double time %e seconds\n", diff_nano(&start, &end));

    return 0;
}
