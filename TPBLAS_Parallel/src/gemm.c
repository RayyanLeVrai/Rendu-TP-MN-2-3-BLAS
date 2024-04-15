#include "mnblas.h"
#include "complexe2.h"
#include <omp.h>


void mncblas_sgemm(MNCBLAS_LAYOUT layout, MNCBLAS_TRANSPOSE TransA,
                 MNCBLAS_TRANSPOSE TransB, const int M, const int N,
                 const int K, const float alpha, const float *A,
                 const int lda, const float *B, const int ldb,
                 const float beta, float *C, const int ldc)
{   
    #pragma omp parallel for
    for (int i = 0; i < M; ++i) {
        for (int j = 0; j < N; ++j) {
            C[i*ldc + j] *= beta;
        }
    }
    #pragma omp parallel for
    for (int i = 0; i < M; ++i) {
        for (int k = 0; k < K; ++k) {
            float temp = alpha * A[i*lda + k];
            for (int j = 0; j < N; ++j) {
                C[i*ldc + j] += temp * B[k*ldb + j];
            }
        }
    }
}                 




void mncblas_dgemm(MNCBLAS_LAYOUT layout, MNCBLAS_TRANSPOSE TransA,
                 MNCBLAS_TRANSPOSE TransB, const int M, const int N,
                 const int K, const double alpha, const double *A,
                 const int lda, const double *B, const int ldb,
                 const double beta, double *C, const int ldc)
{
    #pragma omp parallel for
    for (int i = 0; i < M; ++i) {
        for (int j = 0; j < N; ++j) {
            C[i*ldc + j] *= beta;
        }
    }
    #pragma omp parallel for
    for (int i = 0; i < M; ++i) {
        for (int k = 0; k < K; ++k) {
            double temp = alpha * A[i*lda + k];
            for (int j = 0; j < N; ++j) {
                C[i*ldc + j] += temp * B[k*ldb + j];
            }
        }
    }
}







void mncblas_cgemm(MNCBLAS_LAYOUT layout, MNCBLAS_TRANSPOSE TransA,
                 MNCBLAS_TRANSPOSE TransB, const int M, const int N,
                 const int K, const void *alpha, const void *A,
                 const int lda, const void *B, const int ldb,
                 const void *beta, void *C, const int ldc)
{
    // Conversion des pointeurs void en pointeurs vers complexe_float_t
    const complexe_float_t *A_c = (const complexe_float_t *)A;
    const complexe_float_t *B_c = (const complexe_float_t *)B;
    complexe_float_t *C_c = (complexe_float_t *)C;
    const complexe_float_t alpha_c = *(const complexe_float_t *)alpha;
    const complexe_float_t beta_c = *(const complexe_float_t *)beta;

    // Parcours des lignes de la matrice résultante C
    #pragma omp parallel for
    for (int i = 0; i < M; i++) {
        // Parcours des colonnes de la matrice résultante C
        for (int j = 0; j < N; j++) {
            // Initialisation du nombre complexe temporaire pour le calcul de C[i,j]
            complexe_float_t temp = {0.0, 0.0};

            // Calcul de la contribution de chaque élément
            for (int k = 0; k < K; k++) {
                // Récupération des éléments A[i,k] et B[k,j]
                complexe_float_t a = A_c[i * lda + k];
                complexe_float_t b = B_c[k * ldb + j];

                // Multiplication des nombres complexes a et b
                complexe_float_t prod = {a.real * b.real - a.imaginary * b.imaginary,
                                         a.real * b.imaginary + a.imaginary * b.real};
                // Accumulation du produit dans temp
                temp.real += prod.real;
                temp.imaginary += prod.imaginary;
            }

            // Récupération de l'élément actuel de C pour le mettre à jour
            complexe_float_t c_old = C_c[i * ldc + j];
            // Mise à jour de C[i,j] avec le résultat temporel, multiplié par alpha, plus beta*C[i,j]
            C_c[i * ldc + j].real = alpha_c.real * temp.real - alpha_c.imaginary * temp.imaginary +
                                    beta_c.real * c_old.real - beta_c.imaginary * c_old.imaginary;
            C_c[i * ldc + j].imaginary = alpha_c.real * temp.imaginary + alpha_c.imaginary * temp.real +
                                         beta_c.real * c_old.imaginary + beta_c.imaginary * c_old.real;
        }
    }
}











void mncblas_zgemm(MNCBLAS_LAYOUT layout, MNCBLAS_TRANSPOSE TransA,
                 MNCBLAS_TRANSPOSE TransB, const int M, const int N,
                 const int K, const void *alpha, const void *A,
                 const int lda, const void *B, const int ldb,
                 const void *beta, void *C, const int ldc)
{

    const complexe_double_t *A_c = (const complexe_double_t *)A;
    const complexe_double_t *B_c = (const complexe_double_t *)B;
    complexe_double_t *C_c = (complexe_double_t *)C;
    const complexe_double_t alpha_c = *(const complexe_double_t *)alpha;
    const complexe_double_t beta_c = *(const complexe_double_t *)beta;

    
    for (int i = 0; i < M; i++) {
        
        for (int j = 0; j < N; j++) {
            
            complexe_double_t temp = {0.0, 0.0};

            
            for (int k = 0; k < K; k++) {
                
                complexe_double_t a = A_c[i * lda + k];
                complexe_double_t b = B_c[k * ldb + j];

                
                complexe_double_t prod = {a.real * b.real - a.imaginary * b.imaginary,
                                         a.real * b.imaginary + a.imaginary * b.real};
            
                temp.real += prod.real;
                temp.imaginary += prod.imaginary;
            }

            
            complexe_double_t c_old = C_c[i * ldc + j];

            C_c[i * ldc + j].real = alpha_c.real * temp.real - alpha_c.imaginary * temp.imaginary +
                                    beta_c.real * c_old.real - beta_c.imaginary * c_old.imaginary;
            C_c[i * ldc + j].imaginary = alpha_c.real * temp.imaginary + alpha_c.imaginary * temp.real +
                                         beta_c.real * c_old.imaginary + beta_c.imaginary * c_old.real;
        }
    }
}                