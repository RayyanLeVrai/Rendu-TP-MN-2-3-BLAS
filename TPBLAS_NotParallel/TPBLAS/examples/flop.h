#include <time.h>
#include <sys/time.h>
#include <x86intrin.h>

#define TOP_NANO(x)  clock_gettime (CLOCK_MONOTONIC, &x) ;

double diff_nano (struct timespec *start, struct timespec *end) ;

void init_nano () ;

