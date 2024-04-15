#include <stdio.h>
#include <stdlib.h>

#include "mnblas.h"
#include "complexe2.h"

#include "flop.h"

#define NB_FOIS 512

int main (int argc, char **argv)
{
    complexe_float_t c1 = {1.0, 2.0};
    complexe_float_t c2 = {3.0, 6.0};

    complexe_double_t cd1 = {10.0, 7.0};
    complexe_double_t cd2 = {25.0, 32.0};

    complexe_double_t c1_d; 
    complexe_float_t c2_f;  
    complexe_double_t c3_d; 
    complexe_float_t c4_f;  
    complexe_double_t c5_d; 
    complexe_float_t c6_f;  

    struct timespec start, end;
    int i;

    printf("\n");

    printf("\n");
    printf("Addition double :\n");
    init_nano();
    TOP_NANO(start);

    for (i = 0; i < NB_FOIS; i++) {
        c1_d = add_complexe_double(cd1, cd2);
    }

    TOP_NANO(end);
    printf ("calcul complexe nano %d\ntemps %e seconde\n ", NB_FOIS*2, diff_nano (&start, &end)) ;


    printf("\n");
    printf("Addition float :\n");
    init_nano();
    TOP_NANO(start);

    for (i = 0; i < NB_FOIS; i++) {
        c2_f = add_complexe_float(c1, c2);
    }

    TOP_NANO(end);
    printf ("calcul complexe nano %d\ntemps %e seconde\n ", NB_FOIS*2, diff_nano (&start, &end)) ;


    printf("\n");
    printf("Multiplication double :\n");
    init_nano();
    TOP_NANO(start);

    for (i = 0; i < NB_FOIS; i++) {
        c3_d = mult_complexe_double(cd1, cd2);
    }

    TOP_NANO(end);
    printf ("calcul complexe nano %d\ntemps %e seconde\n ", NB_FOIS*2, diff_nano (&start, &end)) ;


    printf("\n");
    printf("Multiplication float :\n");
    init_nano();
    TOP_NANO(start);

    for (i = 0; i < NB_FOIS; i++) {
        c4_f = mult_complexe_float(c1, c2);
    }

    TOP_NANO(end);
    printf ("calcul complexe nano %d\ntemps %e seconde\n ", NB_FOIS*2, diff_nano (&start, &end)) ;


    printf("\n");
    printf("Division double:\n");
    init_nano();
    TOP_NANO(start);

    for (i = 0; i < NB_FOIS; i++) {
        c5_d = div_complexe_double(cd1, cd2);
    }

    TOP_NANO(end);
    printf ("calcul complexe nano %d\ntemps %e seconde\n ", NB_FOIS*2, diff_nano (&start, &end)) ;


    printf("\n");
    printf("Division float:\n");
    init_nano();
    TOP_NANO(start);

    for (i = 0; i < NB_FOIS; i++) {
        c6_f = div_complexe_float(c1, c2);
    }

    TOP_NANO(end);
    printf ("calcul complexe nano %d\ntemps %e seconde\n ", NB_FOIS*2, diff_nano (&start, &end)) ;


    printf("\n");

    exit(0);
}
