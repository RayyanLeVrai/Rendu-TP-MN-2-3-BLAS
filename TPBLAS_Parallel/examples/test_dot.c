#include "mnblas.h"
#include "complexe.h"
#include "flop.h"
#include <stdio.h>
#define VECSIZE    65536
#define NB_FOIS    10

typedef float vfloat [VECSIZE];
typedef double vdouble [VECSIZE];
typedef complexe_float_t vcomplex [VECSIZE];
typedef complexe_double_t zcomplex [VECSIZE];

vfloat vec1, vec2;
vdouble vec1d, vec2d;
vcomplex vec1c, vec2c;
zcomplex vec1z, vec2z;

complexe_float_t cdotu, cdotc;
complexe_double_t zdotu, zdotc;

void vector_init_float (vfloat V, float x) {
    for (unsigned int i = 0; i < VECSIZE; i++) {
        V [i] = x;
    }
}

void vector_init_double (vdouble V, double x) {
    for (unsigned int i = 0; i < VECSIZE; i++) {
        V [i] = x;
    }
}

void vector_init_complex_float (vcomplex V, float real, float imag) {
    for (unsigned int i = 0; i < VECSIZE; i++) {
        V[i].real = real;
        V[i].imaginary = imag;
    }
}

void vector_init_complex_double (zcomplex V, double real, double imag) {
    for (unsigned int i = 0; i < VECSIZE; i++) {
        V[i].real = real;
        V[i].imaginary = imag;
    }
}

int main (int argc, char **argv) {
    struct timespec start, end;
    float res;
    double resd;
    init_nano();

    vector_init_float(vec1, 1.0);
    vector_init_float(vec2, 2.0);

    vector_init_double(vec1d, 1.0);
    vector_init_double(vec2d, 2.0);

    vector_init_complex_float(vec1c, 1.0, 0.0);
    vector_init_complex_float(vec2c, 2.0, 0.0);

    vector_init_complex_double(vec1z, 1.0, 0.0);
    vector_init_complex_double(vec2z, 2.0, 0.0);

    for (int i = 0; i < NB_FOIS; i++) {
        TOP_NANO(start);
        res = mncblas_sdot(VECSIZE, vec1, 1, vec2, 1);
        TOP_NANO(end);
        printf("sdot nano time %e seconds\n", diff_nano(&start, &end));

        TOP_NANO(start);
        resd = mncblas_ddot(VECSIZE, vec1d, 1, vec2d, 1);
        TOP_NANO(end);
        printf("ddot nano time %e seconds\n", diff_nano(&start, &end));

        TOP_NANO(start);
        mncblas_cdotu_sub(VECSIZE, vec1c, 1, vec2c, 1, &cdotu);
        TOP_NANO(end);
        printf("cdotu nano time %e seconds\n", diff_nano(&start, &end));

        TOP_NANO(start);
        mncblas_cdotc_sub(VECSIZE, vec1c, 1, vec2c, 1, &cdotc);
        TOP_NANO(end);
        printf("cdotc nano time %e seconds\n", diff_nano(&start, &end));

        TOP_NANO(start);
        mncblas_zdotu_sub(VECSIZE, vec1z, 1, vec2z, 1, &zdotu);
        TOP_NANO(end);
        printf("zdotu nano time %e seconds\n", diff_nano(&start, &end));

        TOP_NANO(start);
        mncblas_zdotc_sub(VECSIZE, vec1z, 1, vec2z, 1, &zdotc);
        TOP_NANO(end);
        printf("zdotc nano time %e seconds\n", diff_nano(&start, &end));
    }

    printf("==========================================================\n");

    printf("sdot result = %f\n", res);
    printf("ddot result = %f\n", resd);
    printf("cdotu result = (%f, %f)\n", cdotu.real, cdotu.imaginary);
    printf("cdotc result = (%f, %f)\n", cdotc.real, cdotc.imaginary);
    printf("zdotu result = (%f, %f)\n", zdotu.real, zdotu.imaginary);
    printf("zdotc result = (%f, %f)\n", zdotc.real, zdotc.imaginary);

    return 0;
}
