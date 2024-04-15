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

vfloat Xf;
vdouble Xd;
vcomplex Xc;
zcomplex Xz;

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
    float sasum;
    double dasum;
    float scasum;
    double dzasum;

    init_nano();

    // Tests mnblas_sasum
    vector_init_f(Xf, 1.0);
    for (int i = 0; i < NB_FOIS; i++) {
        TOP_NANO(start);
        sasum = mnblas_sasum(VECSIZE, Xf, 1);
        TOP_NANO(end);
        printf("sasum nano time %e seconds\n", diff_nano(&start, &end));
    }

    // Tests mnblas_dasum
    vector_init_d(Xd, 1.0);
    for (int i = 0; i < NB_FOIS; i++) {
        TOP_NANO(start);
        dasum = mnblas_dasum(VECSIZE, Xd, 1);
        TOP_NANO(end);
        printf("dasum nano time %e seconds\n", diff_nano(&start, &end));
    }

    // Tests mnblas_scasum
    vector_init_c(Xc, 1.0, 0.0);
    for (int i = 0; i < NB_FOIS; i++) {
        TOP_NANO(start);
        scasum = mnblas_scasum(VECSIZE, Xc, 1);
        TOP_NANO(end);
        printf("scasum nano time %e seconds\n", diff_nano(&start, &end));
    }

    // Tests mnblas_dzasum
    vector_init_z(Xz, 1.0, 0.0);
    for (int i = 0; i < NB_FOIS; i++) {
        TOP_NANO(start);
        dzasum = mnblas_dzasum(VECSIZE, Xz, 1);
        TOP_NANO(end);
        printf("dzasum nano time %e seconds\n", diff_nano(&start, &end));
    }

    return 0;
}
