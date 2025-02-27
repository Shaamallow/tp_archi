#include "operation.h"
#include "utils.h"
#include <immintrin.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>

#define N (1024 * 1024) // Large value for benchmarking
#define NUM_THREADS 8   // Number of threads for parallel execution

typedef struct {
  float *U;
  float *V;
  int start;
  int end;
  double result;
} ThreadData;

pthread_mutex_t mutex =
    PTHREAD_MUTEX_INITIALIZER; // Mutex for result accumulation

// Threaded scalar distance computation with mutex
void *thread_func(void *arg) { pthread_exit(NULL); }

// Threaded vectorized distance computation
void *thread_func_vec(void *arg) { pthread_exit(NULL); }

// Parallel computation using pthreads
void distPar(double *U, double *V, int n, int nb_threads, int mode,
             double *result) {}
int main() {
  double *U = (double *)aligned_alloc(32, N * sizeof(double));
  double *V = (double *)aligned_alloc(32, N * sizeof(double));

  initialize_random(U, N);
  initialize_random(V, N);

  double seq_result, par_result_scalar, seq_result_vectorial,
      par_result_vectorial;
  double start, end;

  // Sequential scalar execution
  start = now();
  seq_result = dist(U, V, N);
  end = now();
  double time_seq = end - start;

  // Parallel scalar execution
  start = now();
  distPar(U, V, N, NUM_THREADS, 0, &par_result_scalar);
  end = now();
  double time_par_scalar = end - start;

  // Sequential vectorized execution
  start = now();
  seq_result_vectorial = vect_dist(U, V, N);
  end = now();
  double time_seq_vectorial = end - start;

  // Parallel vectorized execution
  start = now();
  distPar(U, V, N, NUM_THREADS, 1, &par_result_vectorial);
  end = now();
  double time_par_vectorial = end - start;

  printf("Sequential Scalar Result: %f\n", seq_result);
  printf("Parallel Scalar Result: %f\n", par_result_scalar);
  printf("Speedup Scalar: %fx\n", time_seq / time_par_scalar);

  printf("Sequential Vectorial Result: %f\n", seq_result_vectorial);
  printf("Parallel Vectorial Result: %f\n", par_result_vectorial);
  printf("Speedup Vectorial: %fx\n", time_seq_vectorial / time_par_vectorial);

  free(U);
  free(V);
  return 0;
}
