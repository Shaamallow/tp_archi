// utils.h
#ifndef OPERATION_H
#define OPERATION_H

#include <immintrin.h>
#include <math.h>

// Sequential distance computation
double dist(float *U, float *V, int n);

// Vectorized (AVX) version assuming U and V are aligned and n is multiple of 8
double vect_dist(float *U, float *V, int n);

// Generalized vectorized version with relaxed alignment
double vect_dist_gen(float *U, float *V, int n);

#endif
