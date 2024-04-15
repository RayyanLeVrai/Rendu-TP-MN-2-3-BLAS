#include "mnblas.h"
#include "complexe2.h"
#include "flop.h"
#include <stdio.h>

#define VECSIZE 65536
#define NB_FOIS 10

typedef float vfloat[VECSIZE];
typedef double vdouble[VECSIZE];
typedef complexe_float_t vcomplex[VECSIZE];
typedef complexe_double_t zcomplex[VECSIZE];

vfloat Xf, Yf;
vdouble Xd, Yd;
vcomplex Xc, Yc;
zcomplex Xz, Yz;

float alpha_f = 2.0;
double alpha_d = 2.0;
complexe_float_t alpha_c = {2.0, 3.0};
complexe_double_t alpha_z = {2.0, 3.0};

void vector_init_f(vfloat V, float val) {
    for (int i = 0; i < VECSIZE; i++) {
        V[i] = val;
    }
}

void vector_init_d(vdouble V, double val) {
    for (int i = 0; i < VECSIZE; i++) {
        V[i] = val;
    }
}

void vector_init_c(vcomplex V, float real, float imag) {
    for (int i = 0; i < VECSIZE; i++) {
        V[i].real = real;
        V[i].imaginary = imag;
    }
}

void vector_init_z(zcomplex V, double real, double imag) {
    for (int i = 0; i < VECSIZE; i++) {
        V[i].real = real;
        V[i].imaginary = imag;
    }
}

int main() {
    struct timespec start, end;
    init_nano();

    // Tests mnblas_saxpy
    vector_init_f(Xf, 1.0);
    vector_init_f(Yf, 0.0);
    for (int i = 0; i < NB_FOIS; i++) {
        TOP_NANO(start);
        mnblas_saxpy(VECSIZE, alpha_f, Xf, 1, Yf, 1);
        TOP_NANO(end);
        printf("saxpy nano time %e seconds\n", diff_nano(&start, &end));
    }

    // Tests mnblas_daxpy
    vector_init_d(Xd, 1.0);
    vector_init_d(Yd, 0.0);
    for (int i = 0; i < NB_FOIS; i++) {
        TOP_NANO(start);
        mnblas_daxpy(VECSIZE, alpha_d, Xd, 1, Yd, 1);
        TOP_NANO(end);
        printf("daxpy nano time %e seconds\n", diff_nano(&start, &end));
    }

    // Tests mnblas_caxpy
    vector_init_c(Xc, 1.0, 0.0);
    vector_init_c(Yc, 0.0, 0.0);
    for (int i = 0; i < NB_FOIS; i++) {
        TOP_NANO(start);
        mnblas_caxpy(VECSIZE, &alpha_c, Xc, 1, Yc, 1);
        TOP_NANO(end);
        printf("caxpy nano time %e seconds\n", diff_nano(&start, &end));
    }

    // Tests mnblas_zaxpy
    vector_init_z(Xz, 1.0, 0.0);
    vector_init_z(Yz, 0.0, 0.0);
    for (int i = 0; i < NB_FOIS; i++) {
        TOP_NANO(start);
        mnblas_zaxpy(VECSIZE, &alpha_z, Xz, 1, Yz, 1);
        TOP_NANO(end);
        printf("zaxpy nano time %e seconds\n", diff_nano(&start, &end));
    }

    return 0;
}
