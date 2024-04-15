all: test_complexe4 test_complexe3 test_complexe2 test_complexe test_dot test_dot_dyn test_axpy test_asum test_copy test_swap test_gemm test_gemv

INC_DIR=../include
LIB_DIR=../lib

# Include the math library `-lm` for both static and dynamic linking scenarios
LIBST=-lmnblas -lm
LIBDYN=-lmnblasdyn -lm

OPTIONS_COMPIL=-Wall -O2 -fPIC -I$(INC_DIR)
OPTIONS_LINK_STATIC=-fopenmp -L$(LIB_DIR) $(LIBST)
OPTIONS_LINK_DYN=-fopenmp -L$(LIB_DIR) $(LIBDYN)

test_dot: test_dot.o flop.o $(LIB_DIR)/libmnblas.a
	gcc -o test_dot test_dot.o flop.o $(OPTIONS_LINK_STATIC)

flop.o: flop.c flop.h
	gcc $(OPTIONS_COMPIL) -c flop.c

test_dot_dyn: test_dot.o flop.o
	gcc -o test_dot_dyn test_dot.o flop.o $(OPTIONS_LINK_DYN)

test_dot.o: test_dot.c $(INC_DIR)/mnblas.h
	gcc $(OPTIONS_COMPIL) -c test_dot.c

test_axpy.o: test_axpy.c $(INC_DIR)/mnblas.h
	gcc $(OPTIONS_COMPIL) -c test_axpy.c

test_axpy: test_axpy.o flop.o
	gcc -o test_axpy test_axpy.o flop.o $(OPTIONS_LINK_STATIC)

test_asum.o: test_asum.c $(INC_DIR)/mnblas.h
	gcc $(OPTIONS_COMPIL) -c test_asum.c

test_asum: test_asum.o flop.o
	gcc -o test_asum test_asum.o flop.o $(OPTIONS_LINK_STATIC)

test_copy.o: test_copy.c $(INC_DIR)/mnblas.h
	gcc $(OPTIONS_COMPIL) -c test_copy.c

test_copy: test_copy.o flop.o
	gcc -o test_copy test_copy.o flop.o $(OPTIONS_LINK_STATIC)

test_swap.o: test_swap.c $(INC_DIR)/mnblas.h
	gcc $(OPTIONS_COMPIL) -c test_swap.c

test_swap: test_swap.o flop.o
	gcc -o test_swap test_swap.o flop.o $(OPTIONS_LINK_STATIC)

test_complexe.o: test_complexe.c $(INC_DIR)/complexe.h
	gcc $(OPTIONS_COMPIL) -c test_complexe.c

test_complexe: test_complexe.o flop.o
	gcc -o test_complexe test_complexe.o flop.o $(OPTIONS_LINK_STATIC)

test_complexe2.o: test_complexe2.c $(INC_DIR)/complexe2.h
	gcc $(OPTIONS_COMPIL) -c test_complexe2.c

test_complexe2: test_complexe2.o flop.o
	gcc -o test_complexe2 test_complexe2.o flop.o $(OPTIONS_LINK_STATIC)

test_complexe3.o: test_complexe3.c $(INC_DIR)/complexe2.h
	gcc $(OPTIONS_COMPIL) -c test_complexe3.c

test_complexe3: test_complexe3.o flop.o
	gcc -o test_complexe3 test_complexe3.o flop.o $(OPTIONS_LINK_STATIC)

test_complexe4.o: test_complexe4.c $(INC_DIR)/complexe2.h
	gcc $(OPTIONS_COMPIL) -c test_complexe4.c

test_complexe4: test_complexe4.o flop.o
	gcc -o test_complexe4 test_complexe4.o flop.o $(OPTIONS_LINK_STATIC)

test_gemm.o: test_gemm.c $(INC_DIR)/mnblas.h
	gcc $(OPTIONS_COMPIL) -c test_gemm.c

test_gemm: test_gemm.o flop.o
	gcc -o test_gemm test_gemm.o flop.o $(OPTIONS_LINK_STATIC)

test_gemv.o: test_gemv.c $(INC_DIR)/mnblas.h
	gcc $(OPTIONS_COMPIL) -c test_gemv.c

test_gemv: test_gemv.o flop.o
	gcc -o test_gemv test_gemv.o flop.o $(OPTIONS_LINK_STATIC)

clean:
	rm -f *.o test_dot test_dot_dyn test_axpy test_asum test_copy test_swap test_complexe test_complexe2 test_complexe3 test_complexe4 test_gemm test_gemv *~
