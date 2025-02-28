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

// Threaded scalar distance computation
void *thread_func(void *arg) {
  ThreadData *data = (ThreadData *)arg;
  data->result = 0.0;

  // Add offset to the pointers to reuse the dist() functions
  float *U = data->U + data->start;
  float *V = data->V + data->start;

  // Compute the distance
  data->result = dist(U, V, data->end - data->start);

  pthread_exit(NULL);
}

// Threaded vectorized distance computation
void *thread_func_vec(void *arg) {
  ThreadData *data = (ThreadData *)arg;
  data->result = 0.0;

  // Add offset to the pointers to reuse vect_dist() functions
  float *U = data->U + data->start;
  float *V = data->V + data->start;

  // Compute the distance
  data->result = vect_dist(U, V, data->end - data->start);

  pthread_exit(NULL);
}

// Parallel computation using pthreads
void distPar(float *U, float *V, int n, int nb_threads, int mode,
             double *result) {
  pthread_t threads[nb_threads];
  ThreadData data[nb_threads];
  int chunk = n / nb_threads;

  if (mode == 0) {
    for (int i = 0; i < nb_threads; i++) {
      data[i].U = U;
      data[i].V = V;
      data[i].start = i * chunk;
      data[i].end = (i == nb_threads - 1) ? n : (i + 1) * chunk;
      pthread_create(&threads[i], NULL, thread_func, &data[i]);
    }
  } else {
    for (int i = 0; i < nb_threads; i++) {
      data[i].U = U;
      data[i].V = V;
      data[i].start = i * chunk;
      data[i].end = (i == nb_threads - 1) ? n : (i + 1) * chunk;
      pthread_create(&threads[i], NULL, thread_func_vec, &data[i]);
    }
  }

  for (int i = 0; i < nb_threads; i++) {
    pthread_join(threads[i], NULL);
  }

  *result = 0.0;
  for (int i = 0; i < nb_threads; i++) {
    *result += data[i].result;
  }
}

int main() {
  float *U = (float *)aligned_alloc(32, N * sizeof(float));
  float *V = (float *)aligned_alloc(32, N * sizeof(float));

  double *U_double = (double *)aligned_alloc(32, N * sizeof(double));
  double *V_double = (double *)aligned_alloc(32, N * sizeof(double));

  initialize_random(U, N);
  initialize_random(V, N);

  copy_cast_array(U, U_double, N);
  copy_cast_array(V, V_double, N);

  double seq_result, par_result_scalar, seq_result_vectorial,
      par_result_vectorial, seq_result_vectorial_double;
  double start, end;

  // Sequential scalar execution
  start = now();
  seq_result = dist(U, V, N);
  end = now();
  double time_seq_scalar = end - start;

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

  // Fixed rounding error in Sequential Vectorized execution
  start = now();
  seq_result_vectorial_double = vect_dist_double(U_double, V_double, N);
  end = now();
  double time_seq_vectorial_double = end - start;

  // Parallel vectorized execution
  start = now();
  distPar(U, V, N, NUM_THREADS, 1, &par_result_vectorial);
  end = now();
  double time_par_vectorial = end - start;

  printf("Sequential Scalar Result: %fs\n", time_seq_scalar);
  printf("Sequential Vectorial Result: %fs\n", time_seq_vectorial);
  printf("Sequential Vectorial Double Result: %fs\n",
         time_seq_vectorial_double);
  printf("Parallel Scalar Result: %fs\n", time_par_scalar);
  printf("Parallel Vectorial Result: %fs\n", time_par_vectorial);

  // Get results values
  printf("\n---- Checking Acceleration ----\n");
  printf("Vectorial over Sequential: %fx\n",
         time_seq_scalar / time_seq_vectorial);
  printf("Vectorial over Sequential (with double): %fx\n",
         time_seq_scalar / time_seq_vectorial_double);
  printf("Vectorial over Threading: %fx\n",
         time_par_scalar / time_par_vectorial);
  printf("Threading over Scalar: %fx\n", time_seq_scalar / time_par_scalar);
  printf("Threading over Vectorial: %fx\n",
         time_seq_vectorial / time_par_vectorial);

  printf("\n---- Checking Results ----\n");
  printf("Sequential Scalar Result: %f\n", seq_result);
  printf("Parallel Scalar Result: %f\n", par_result_scalar);
  printf("Sequential Vectorial Result: %f\n", seq_result_vectorial);
  printf("Parallel Vectorial Result: %f\n", par_result_vectorial);
  printf("Sequential Vectorial Double Result: %f\n",
         seq_result_vectorial_double);

  free(U);
  free(V);
  return 0;
}
