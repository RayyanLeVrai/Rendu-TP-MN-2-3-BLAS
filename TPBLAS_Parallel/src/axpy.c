#include "mnblas.h"
#include "complexe2.h"
#include <omp.h>
#include <stdio.h>
#include <stddef.h>

void mnblas_saxpy(const int N, const float alpha, const float *X,
                 const int incX, float *Y, const int incY)
{
    if( N<0 || alpha == 0.0)
    {
        return;
    }

    if (X==NULL || Y==NULL)
    {
        return;
    }

    int i;

    if (incX == 1 && incY ==1)
    {
        #pragma omp parallel for
        for (i=0 ; i<N ; i++)
        {
            Y[i] = alpha*X[i] + Y[i];
        }
    }  

    else {
        int ix = 0;
        int iy = 0;

        if (incX < 0) {
            ix = (-N + 1) * incX;
        }
        if (incY < 0) {
            iy = (-N + 1) * incY;
        }
        #pragma omp parallel for
        for (i = 0; i < N; i++) {
            Y[iy] = alpha * X[ix] + Y[iy];
            ix += incX;
            iy += incY;
        }
    }
}


void mnblas_daxpy(const int N, const double alpha, const double *X,
                 const int incX, double *Y, const int incY)
{
    if( N<0 || alpha == 0.0)
    {
        return;
    }

    if (X==NULL || Y==NULL)
    {
        return;
    }

    int i;

    if (incX == 1 && incY ==1)
    {
        #pragma omp parallel for
        for (i=0 ; i<N ; i++)
        {
            Y[i] = alpha*X[i] + Y[i];
        }
    }  

    else {
        int ix = 0;
        int iy = 0;

        if (incX < 0) {
            ix = (-N + 1) * incX;
        }
        if (incY < 0) {
            iy = (-N + 1) * incY;
        }
        #pragma omp parallel for
        for (i = 0; i < N; i++) {
            Y[iy] = alpha * X[ix] + Y[iy];
            ix += incX;
            iy += incY;
        }
    }
}                 

void mnblas_caxpy(const int N, const void *alpha, const void *X, const int incX, void *Y, const int incY) {
	

	// on convertit les pointeurs
	 	   
    const complexe_float_t *alpha_c = (const complexe_float_t *)alpha;
    const complexe_float_t *X_c = (const complexe_float_t *)X;
    complexe_float_t *Y_c = (complexe_float_t *)Y;

    #pragma omp parallel for
    for (int i = 0; i < N; i++) {
        int ix = i * incX;
        int iy = i * incY;
        
        
        complexe_float_t temp;
        temp.real = alpha_c->real * X_c[ix].real - alpha_c->imaginary * X_c[ix].imaginary;
        temp.imaginary = alpha_c->real * X_c[ix].imaginary + alpha_c->imaginary * X_c[ix].real;
        
        Y_c[iy].real += temp.real;
        Y_c[iy].imaginary += temp.imaginary;
    }
}

void mnblas_zaxpy(const int N, const void *alpha, const void *X, const int incX, void *Y, const int incY) {
    
    const complexe_double_t *alpha_c = (const complexe_double_t *)alpha;
    const complexe_double_t *X_c = (const complexe_double_t *)X;
    complexe_double_t *Y_c = (complexe_double_t *)Y;

    #pragma omp parallel for
    for (int i = 0; i < N; i++) {
        int ix = i * incX;
        int iy = i * incY;
        
        complexe_double_t temp;
        temp.real = alpha_c->real * X_c[ix].real - alpha_c->imaginary * X_c[ix].imaginary;
        temp.imaginary = alpha_c->real * X_c[ix].imaginary + alpha_c->imaginary * X_c[ix].real;
        
        Y_c[iy].real += temp.real;
        Y_c[iy].imaginary += temp.imaginary;
    }
}                 
