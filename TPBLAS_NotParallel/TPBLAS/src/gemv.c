#include "mnblas.h"
#include "complexe2.h" 
#include <omp.h>


void mncblas_sgemv(const MNCBLAS_LAYOUT layout,
                   const MNCBLAS_TRANSPOSE TransA, const int M, const int N,
                   const float alpha, const float *A, const int lda,
                   const float *X, const int incX, const float beta,
                   float *Y, const int incY) {

        for (int i = 0; i < M; i++) {
            Y[i * incY] *= beta;
        }

        if (TransA == MNCblasNoTrans) {
            for (int i = 0; i < M; i++) {
                float temp = 0.0;
                for (int j = 0; j < N; j++) {
                    temp += A[i * lda + j] * X[j * incX];
                }
                Y[i * incY] += alpha * temp;
            }
        } else {
            for (int i = 0; i < N; i++) {
                float temp = 0.0;
                for (int j = 0; j < M; j++) {
                    temp += A[j * lda + i] * X[j * incX];
                }
                Y[i * incY] += alpha * temp;
            }
        }
    
}



void mncblas_dgemv(const MNCBLAS_LAYOUT layout,
                   const MNCBLAS_TRANSPOSE TransA, const int M, const int N,
                   const double alpha, const double *A, const int lda,
                   const double *X, const int incX, const double beta,
                   double *Y, const int incY) {
    for (int i = 0; i < M; i++) {
        Y[i * incY] *= beta;
    }

    if (TransA == MNCblasNoTrans) {
        for (int i = 0; i < M; i++) {
            double temp = 0.0;
            for (int j = 0; j < N; j++) {
                temp += A[i * lda + j] * X[j * incX];
            }
            Y[i * incY] += alpha * temp;
        }
    } else {
        for (int i = 0; i < N; i++) {
            double temp = 0.0;
            for (int j = 0; j < M; j++) {
                temp += A[j * lda + i] * X[j * incX];
            }
            Y[i * incY] += alpha * temp;
        }
    }
}



void mncblas_cgemv(MNCBLAS_LAYOUT layout,
                   MNCBLAS_TRANSPOSE TransA, const int M, const int N,
                   const void *alpha, const void *A, const int lda,
                   const void *X, const int incX, const void *beta,
                   void *Y, const int incY) {

    const complexe_float_t *A_c = (const complexe_float_t *)A;
    const complexe_float_t *X_c = (const complexe_float_t *)X;
    complexe_float_t *Y_c = (complexe_float_t *)Y;
    complexe_float_t alpha_c = *(const complexe_float_t *)alpha;
    complexe_float_t beta_c = *(const complexe_float_t *)beta;

    for (int i = 0; i < M; i++) {
        complexe_float_t y_val = Y_c[i * incY];
        Y_c[i * incY].real = y_val.real * beta_c.real - y_val.imaginary * beta_c.imaginary;
        Y_c[i * incY].imaginary = y_val.real * beta_c.imaginary + y_val.imaginary * beta_c.real;
    }
    for (int i = 0; i < M; i++) {
        complexe_float_t temp = {0.0, 0.0};
        for (int j = 0; j < N; j++) {
            complexe_float_t a_val = A_c[i * lda + j];
            complexe_float_t x_val = X_c[j * incX];
            
            float real_part = a_val.real * x_val.real - a_val.imaginary * x_val.imaginary;
            float imaginary_part = a_val.real * x_val.imaginary + a_val.imaginary * x_val.real;

            
            temp.real += real_part;
            temp.imaginary += imaginary_part;
        }

    
        float temp_real = temp.real * alpha_c.real - temp.imaginary * alpha_c.imaginary;
        float temp_imaginary = temp.real * alpha_c.imaginary + temp.imaginary * alpha_c.real;

        
        Y_c[i * incY].real += temp_real;
        Y_c[i * incY].imaginary += temp_imaginary;
    }
}




void mncblas_zgemv(MNCBLAS_LAYOUT layout,
                   MNCBLAS_TRANSPOSE TransA, const int M, const int N,
                   const void *alpha, const void *A, const int lda,
                   const void *X, const int incX, const void *beta,
                   void *Y, const int incY) {

    const complexe_double_t *A_c = (const complexe_double_t *)A;
    const complexe_double_t *X_c = (const complexe_double_t *)X;
    complexe_double_t *Y_c = (complexe_double_t *)Y;
    complexe_double_t alpha_c = *(const complexe_double_t *)alpha;
    complexe_double_t beta_c = *(const complexe_double_t *)beta;

    for (int i = 0; i < M; i++) {
        complexe_double_t y_val = Y_c[i * incY];
        Y_c[i * incY].real = y_val.real * beta_c.real - y_val.imaginary * beta_c.imaginary;
        Y_c[i * incY].imaginary = y_val.real * beta_c.imaginary + y_val.imaginary * beta_c.real;
    }
    for (int i = 0; i < M; i++) {
        complexe_double_t temp = {0.0, 0.0};
        for (int j = 0; j < N; j++) {
            complexe_double_t a_val = A_c[i * lda + j];
            complexe_double_t x_val = X_c[j * incX];
            
            double real_part = a_val.real * x_val.real - a_val.imaginary * x_val.imaginary;
            double imaginary_part = a_val.real * x_val.imaginary + a_val.imaginary * x_val.real;

            
            temp.real += real_part;
            temp.imaginary += imaginary_part;
        }

        
        float temp_real = temp.real * alpha_c.real - temp.imaginary * alpha_c.imaginary;
        float temp_imaginary = temp.real * alpha_c.imaginary + temp.imaginary * alpha_c.real;

        
        Y_c[i * incY].real += temp_real;
        Y_c[i * incY].imaginary += temp_imaginary;
    }
}