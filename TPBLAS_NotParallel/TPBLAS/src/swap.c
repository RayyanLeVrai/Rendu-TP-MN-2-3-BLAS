#include "mnblas.h"
#include "complexe2.h"

void mncblas_sswap(const int N, float *X, const int incX, float *Y, const int incY) {
  for (int i = 0; i < N; i += incX) {
    int j = (i / incX) * incY;
    if (j < N) {
      float save = Y[j];
      Y[j] = X[i];
      X[i] = save;
    }
  }
}

void mncblas_dswap(const int N, double *X, const int incX, double *Y, const int incY) {
  for (int i = 0; i < N; i += incX) {
    int j = (i / incX) * incY;
    if (j < N) {
      double save = Y[j];
      Y[j] = X[i];
      X[i] = save;
    }
  }
}

void mncblas_cswap(const int N, void *X, const int incX, void *Y, const int incY) {
  for (int i = 0; i < N; i += incX * 2) {
    int j = (i / (incX * 2)) * (incY * 2);
    if (j < N) {
      float *cX = (float *)X;
      float *cY = (float *)Y;
      float save[2];
      save[0] = cY[j];
      save[1] = cY[j + 1];
      cY[j] = cX[i];
      cY[j + 1] = cX[i + 1];
      cX[i] = save[0];
      cX[i + 1] = save[1];
    }
  }
}

void mncblas_zswap(const int N, void *X, const int incX, void *Y, const int incY) {
  for (int i = 0; i < N; i += incX * 2) {
    int j = (i / (incX * 2)) * (incY * 2);
    if (j < N) {
      double *zX = (double *)X;
      double *zY = (double *)Y;
      double save[2];
      save[0] = zY[j];
      save[1] = zY[j + 1];
      zY[j] = zX[i];
      zY[j + 1] = zX[i + 1];
      zX[i] = save[0];
      zX[i + 1] = save[1];
    }
  }
}
