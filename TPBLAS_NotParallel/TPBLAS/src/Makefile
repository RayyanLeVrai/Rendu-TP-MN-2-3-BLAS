# une bibliotheque statique (libmnblas.a)
# une bibliotheque dynamique (libmnblasdyn.so)

all: libmnblas.a libmnblasdyn.so

LIB_DIR=../lib
INC_DIR=../include
OPTIM=-O2 -fopenmp

OPTIONS=-O2 -Wall -fPIC $(OPTIM) -I$(INC_DIR)
# Bibliothèque statique

libmnblas.a: copy.o swap.o dot.o axpy.o complexe.o gemv.o gemm.o asum.o
	rm -f libmnblas.a $(LIB_DIR)/libmnblas.a
	ar -r libmnblas.a copy.o swap.o dot.o axpy.o complexe.o gemv.o gemm.o asum.o
	cp libmnblas.a $(LIB_DIR)

# Bibliothèque dynamique

libmnblasdyn.so: copy.o swap.o dot.o axpy.o complexe.o gemv.o gemm.o asum.o
	rm -f libmnblasdyn.so $(LIB_DIR)/libmnblasdyn.so
	gcc -shared -o libmnblasdyn.so copy.o swap.o dot.o axpy.o complexe.o gemv.o gemm.o asum.o
	cp libmnblasdyn.so $(LIB_DIR)

complexe.o: $(INC_DIR)/complexe.h complexe.c
	gcc $(OPTIONS) -c complexe.c

# BLAS 1

copy.o : copy.c $(INC_DIR)/mnblas.h
	gcc $(OPTIONS) -c copy.c

swap.o : swap.c $(INC_DIR)/mnblas.h
	gcc  $(OPTIONS) -c swap.c

dot.o : dot.c $(INC_DIR)/mnblas.h
	gcc  $(OPTIONS) -c dot.c

axpy.o : axpy.c $(INC_DIR)/mnblas.h
	gcc $(OPTIONS) -c axpy.c

asum.o : asum.c $(INC_DIR)/mnblas.h
	gcc $(OPTIONS) -c asum.c

# BLAS 2

gemv.o : gemv.c $(INC_DIR)/mnblas.h
	gcc $(OPTIONS) -c gemv.c
		
# BLAS 3

gemm.o : gemm.c $(INC_DIR)/mnblas.h
	gcc $(OPTIONS) -c gemm.c

clean:
	rm -f *.o *.a *.so *~  $(INC_DIR)/*~ $(LIB_DIR)/libmnblasdyn.so  $(LIB_DIR)/libmnblas.a ../lib/libmnblas.a ../lib/libmnblasdyn.so
