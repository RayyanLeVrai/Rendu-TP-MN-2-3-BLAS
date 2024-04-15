#include <stdio.h>
#include <x86intrin.h>

#include "flop.h"

// duree du cycle processeur en nano seconde (10e9)

static double residu_nano ;

double diff_nano (struct timespec *start,
	     struct timespec *end)
 {
   return (end->tv_sec - start->tv_sec) +
     1e-9*(end->tv_nsec - start->tv_nsec) ;
 }

void init_nano ()
{
  struct timespec start, end ;

  TOP_NANO(start) ;

  TOP_NANO(end) ;

  residu_nano = diff_nano (&start, &end) ; ;
  printf ("residu_nano = %e seconde \n", residu_nano) ;
}



