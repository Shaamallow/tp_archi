#include <immintrin.h>
#include <math.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#define N (1024 * 1024) // Large value for benchmarking
#define NUM_THREADS 8   // Number of threads for parallel execution

typedef struct {
  float *U;
  float *V;
  int start;
  int end;
  double result;
} ThreadData;

// Function to get the current time in seconds
double now() {
  struct timeval t;
  double ft;
  gettimeofday(&t, NULL);
  ft = t.tv_sec + t.tv_usec / 1000000.0;
  return ft;
}

// Sequential distance computation
double dist(float *U, float *V, int n) {
  double sum = 0.0;
  for (int i = 0; i < n; i++) {
    sum +=
        sqrt((U[i] * U[i] + V[i] * V[i]) / (1 + (U[i] * V[i]) * (U[i] * V[i])));
  }
  return sum;
}

// Vectorized (AVX) version assuming U and V are aligned and n is multiple of 8
double vect_dist(float *U, float *V, int n) {
  __m256 sum_vec = _mm256_setzero_ps();

  for (int i = 0; i < n; i += 8) {
    __m256 u = _mm256_load_ps(&U[i]);
    __m256 v = _mm256_load_ps(&V[i]);

    __m256 u2 = _mm256_mul_ps(u, u);
    __m256 v2 = _mm256_mul_ps(v, v);
    __m256 uv = _mm256_mul_ps(u, v);
    __m256 uv2 = _mm256_mul_ps(uv, uv);
    __m256 num = _mm256_add_ps(u2, v2);
    __m256 denom = _mm256_add_ps(_mm256_set1_ps(1.0f), uv2);
    __m256 frac = _mm256_div_ps(num, denom);
    __m256 res = _mm256_sqrt_ps(frac);

    sum_vec = _mm256_add_ps(sum_vec, res);
  }

  float temp[8];
  _mm256_store_ps(temp, sum_vec);
  double sum = 0.0;
  for (int i = 0; i < 8; i++) {
    sum += temp[i];
  }
  return sum;
}

// Generalized vectorized version with relaxed alignment
double vect_dist_gen(float *U, float *V, int n) {
  double sum = 0.0;
  int i;
  for (i = 0; i <= n - 8; i += 8) {
    __m256 u = _mm256_loadu_ps(&U[i]);
    __m256 v = _mm256_loadu_ps(&V[i]);

    __m256 u2 = _mm256_mul_ps(u, u);
    __m256 v2 = _mm256_mul_ps(v, v);
    __m256 uv = _mm256_mul_ps(u, v);
    __m256 uv2 = _mm256_mul_ps(uv, uv);
    __m256 num = _mm256_add_ps(u2, v2);
    __m256 denom = _mm256_add_ps(_mm256_set1_ps(1.0f), uv2);
    __m256 frac = _mm256_div_ps(num, denom);
    __m256 res = _mm256_sqrt_ps(frac);

    float temp[8];
    _mm256_storeu_ps(temp, res);
    for (int j = 0; j < 8; j++) {
      sum += temp[j];
    }
  }
  return sum;
}

// Threaded distance computation
void *thread_func(void *arg) {
  ThreadData *data = (ThreadData *)arg;
  double sum = 0.0;

  // TODO: REWRITE THIS CODE TO USE THE PREVIOUS SCALAR IMPL
  for (int i = data->start; i < data->end; i++) {
    sum += sqrt((data->U[i] * data->U[i] + data->V[i] * data->V[i]) /
                (1 + (data->U[i] * data->V[i]) * (data->U[i] * data->V[i])));
  }

  // change this to use mutex and lock instead
  data->result = sum;
  pthread_exit(NULL);
}

void *thread_func_vec(void *arg) {
  ThreadData *data = (ThreadData *)arg;
  double sum = 0.0;

  // TODO: REWRITE THIS FUNCTION TO USE VECTORIAL CALCULATION

  pthread_exit(NULL);
}

// Parallel computation using pthreads
void distPar(float *U, float *V, int n, int nb_threads, int mode,
             double *result) {
  pthread_t threads[nb_threads];
  ThreadData thread_data[nb_threads];
  int chunk_size = n / nb_threads;

  if (mode == 0) {
    // TODO: USE scalar calculation
    for (int i = 0; i < nb_threads; i++) {
      thread_data[i].U = U;
      thread_data[i].V = V;
      thread_data[i].start = i * chunk_size;
      thread_data[i].end = (i == nb_threads - 1) ? n : (i + 1) * chunk_size;
      pthread_create(&threads[i], NULL, thread_func, &thread_data[i]);
    }
  }

  if (mode == 1) {
    // TODO: USE vectorial calculation
    for (int i = 0; i < nb_threads; i++) {
      thread_data[i].U = U;
      thread_data[i].V = V;
      thread_data[i].start = i * chunk_size;
      thread_data[i].end = (i == nb_threads - 1) ? n : (i + 1) * chunk_size;
      pthread_create(&threads[i], NULL, thread_func_vec, &thread_data[i]);
    }
  }

  double total = 0.0;
  for (int i = 0; i < nb_threads; i++) {
    pthread_join(threads[i], NULL);
    total += thread_data[i].result;
  }
  *result = total;
}

int main() {
  float *U = (float *)aligned_alloc(32, N * sizeof(float));
  float *V = (float *)aligned_alloc(32, N * sizeof(float));

  srand(time(NULL));
  for (int i = 0; i < N; i++) {
    U[i] = (float)rand() / RAND_MAX;
    V[i] = (float)rand() / RAND_MAX;
  }

  // TODO: Rewrite, cleanup, remove dups
  double seq_result, par_result;

  // NOTE: Sequential, scalar
  double start_time_seq = now();
  seq_result = dist(U, V, N);
  double end_time_seq = now();
  seq_result = end_time_seq - start_time_seq;

  // NOTE:  Parallel, scalar
  double start_time_vect = now();
  distPar(U, V, N, NUM_THREADS, 0, &par_result);
  double end_time_vect = now();
  par_result = end_time_vect - start_time_vect;

  // NOTE:  Sequential, vectorial
  // NOTE:  Parallel, vectorial

  printf("Sequential Result: %f\n", seq_result);
  printf("Parallel Result: %f\n", par_result);
  printf("Speedup: %fx\n", seq_result / par_result);

  free(U);
  free(V);

  return 0;
}
