#include "mnblas.h"
#include "complexe2.h"
#include "flop.h"
#include <stdio.h>
#include <omp.h>

#define VECSIZE 65536
#define NB_FOIS 10
#define M 128
#define N 128
#define K 128

typedef float vfloat[VECSIZE];
typedef double vdouble[VECSIZE];
typedef complexe_float_t vcomplex[VECSIZE];
typedef complexe_double_t zcomplex[VECSIZE];



float A_f[M*K], B_f[K*N], C_f[M*N];
double A_d[M*K], B_d[K*N], C_d[M*N];
complexe_float_t A_c[M*K], B_c[K*N], C_c[M*N];
complexe_double_t A_z[M*K], B_z[K*N], C_z[M*N];

void matrix_init_f(float* matrix, int rows, int cols, float value) {
    for (int i = 0; i < rows*cols; i++) {
        matrix[i] = value;
    }
}

void matrix_init_d(double* matrix, int rows, int cols, double value) {
    for (int i = 0; i < rows*cols; i++) {
        matrix[i] = value;
    }
}

void matrix_init_c(complexe_float_t* matrix, int rows, int cols, float real, float imag) {
    for (int i = 0; i < rows*cols; i++) {
        matrix[i].real = real;
        matrix[i].imaginary = imag;
    }
}

void matrix_init_z(complexe_double_t* matrix, int rows, int cols, double real, double imag) {
    for (int i = 0; i < rows*cols; i++) {
        matrix[i].real = real;
        matrix[i].imaginary = imag;
    }
}

int main() {
    struct timespec start, end;
    init_nano();

    matrix_init_f(A_f, M, K, 1.0);
    matrix_init_f(B_f, K, N, 2.0);
    matrix_init_f(C_f, M, N, 0.0);

    matrix_init_d(A_d, M, K, 1.0);
    matrix_init_d(B_d, K, N, 2.0);
    matrix_init_d(C_d, M, N, 0.0);

    matrix_init_c(A_c, M, K, 1.0, 0.0);
    matrix_init_c(B_c, K, N, 2.0, 0.0);
    matrix_init_c(C_c, M, N, 0.0, 0.0);

    matrix_init_z(A_z, M, K, 1.0, 0.0);
    matrix_init_z(B_z, K, N, 2.0, 0.0);
    matrix_init_z(C_z, M, N, 0.0, 0.0);

    // Tests mncblas_sgemm
    for (int i = 0; i < NB_FOIS; i++) {
        TOP_NANO(start);
        mncblas_sgemm(1,1,1, M, N, K, 1.0, A_f, K, B_f, N, 0.0, C_f, N);
        TOP_NANO(end);
        printf("sgemm nano time %e seconds\n", diff_nano(&start, &end));
    }

    // Tests mncblas_dgemm
    for (int i = 0; i < NB_FOIS; i++) {
        TOP_NANO(start);
        mncblas_dgemm(1,1,1, M, N, K, 1.0, A_d, K, B_d, N, 0.0, C_d, N);
        TOP_NANO(end);
        printf("dgemm nano time %e seconds\n", diff_nano(&start, &end));
    }

    // Tests mncblas_cgemm
    for (int i = 0; i < NB_FOIS; i++) {
        TOP_NANO(start);
        mncblas_cgemm(1,1,1, M, N, K, &A_c[0], A_c, K, B_c, N, &C_c[0], C_c, N);
        TOP_NANO(end);
        printf("cgemm nano time %e seconds\n", diff_nano(&start, &end));
    }

    // Tests mncblas_zgemm
    for (int i = 0; i < NB_FOIS; i++) {
        TOP_NANO(start);
        mncblas_zgemm(1,1,1, M, N, K, &A_z[0], A_z, K, B_z, N, &C_z[0], C_z, N);
        TOP_NANO(end);
        printf("zgemm nano time %e seconds\n", diff_nano(&start, &end));
    }

    return 0;
}
