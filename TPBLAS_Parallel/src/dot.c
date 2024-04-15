#include "mnblas.h"
#include "complexe2.h"
#include <omp.h>

float mncblas_sdot(const int N, const float *X, const int incX, 
                 const float *Y, const int incY)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;
  float dot = 0.0 ;

  #pragma omp parallel for
  for (i = 0 ; i < N ; i += incX)
    {
      dot += X [i] * Y [j] ;
      j+=incY ;
    }

  return dot ;
}

double mncblas_ddot(const int N, const double *X, const int incX, 
                 const double *Y, const int incY)
{
  register unsigned int i = 0;
  register unsigned int j = 0;

  double dot = 0.0;
  #pragma omp parallel for
  for (i=0 ; i<N ; i+=incX)
  {
    dot += X [i] * Y[j];
    j+=incY;
  }
  
  return dot;
}

void mncblas_cdotu_sub(const int N, const void *X, const int incX, const void *Y, const int incY, void *dotu) {
    
    const complexe_float_t *x = (const complexe_float_t *)X;
    const complexe_float_t *y = (const complexe_float_t *)Y;
    complexe_float_t *dot_product = (complexe_float_t *)dotu;

    
    dot_product->real = 0.0f;
    dot_product->imaginary = 0.0f;

    #pragma omp parallel for
    for (int i = 0; i < N; i++) {
        dot_product->real += x[i * incX].real * y[i * incY].real - x[i * incX].imaginary * y[i * incY].imaginary;
        dot_product->imaginary += x[i * incX].real * y[i * incY].imaginary + x[i * incX].imaginary * y[i * incY].real;
    }
}


void mncblas_cdotc_sub(const int N, const void *X, const int incX, const void *Y, const int incY, void *dotc) {
    
    const complexe_float_t *x = (const complexe_float_t *)X;
    const complexe_float_t *y = (const complexe_float_t *)Y;
    complexe_float_t *dot_product = (complexe_float_t *)dotc;

    
    dot_product->real = 0.0f;
    dot_product->imaginary = 0.0f;

    #pragma omp parallel for
    for (int i = 0; i < N; i++) {
        dot_product->real += x[i * incX].real * y[i * incY].real + x[i * incX].imaginary * y[i * incY].imaginary;
        dot_product->imaginary += x[i * incX].real * y[i * incY].imaginary - x[i * incX].imaginary * y[i * incY].real;
    }
}

void   mncblas_zdotu_sub(const int N, const void *X, const int incX,
                       const void *Y, const int incY, void *dotu)
{
   complexe_double_t dot = {0.0, 0.0};
    #pragma omp parallel for
    for (int i = 0; i < N; i++) {
        complexe_double_t x = ((complexe_double_t *)X)[i * incX];
        complexe_double_t y = ((complexe_double_t *)Y)[i * incY];

        
        complexe_double_t temp;
        temp.real = x.real * y.real - x.imaginary * y.imaginary;
        temp.imaginary = x.real * y.imaginary + x.imaginary * y.real;
        dot.real += temp.real;
        dot.imaginary += temp.imaginary;
    }

    *((complexe_double_t *)dotu) = dot;
}
  
void   mncblas_zdotc_sub(const int N, const void *X, const int incX,
                       const void *Y, const int incY, void *dotc)
{
 complexe_double_t dot = {0.0, 0.0};
    #pragma omp parallel for
    for (int i = 0; i < N; i++) {
        complexe_double_t x = ((complexe_double_t *)X)[i * incX];
        complexe_double_t y = ((complexe_double_t *)Y)[i * incY];

  
        complexe_double_t temp;
        temp.real = x.real * y.real - x.imaginary * (-y.imaginary); // Notice the conjugation of y
        temp.imaginary = x.real * (-y.imaginary) + x.imaginary * y.real;
        dot.real += temp.real;
        dot.imaginary += temp.imaginary;
    }

    *((complexe_double_t *)dotc) = dot;
}
