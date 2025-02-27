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

pthread_mutex_t mutex =
    PTHREAD_MUTEX_INITIALIZER; // Mutex for result accumulation

// Function to get the current time in seconds
double now() {
  struct timeval t;
  gettimeofday(&t, NULL);
  return t.tv_sec + t.tv_usec / 1000000.0;
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

  float temp[8] __attribute__((aligned(32)));
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

// Threaded scalar distance computation with mutex
void *thread_func(void *arg) {
  ThreadData *data = (ThreadData *)arg;
  double local_sum = 0.0;

  for (int i = data->start; i < data->end; i++) {
    local_sum +=
        sqrt((data->U[i] * data->U[i] + data->V[i] * data->V[i]) /
             (1 + (data->U[i] * data->V[i]) * (data->U[i] * data->V[i])));
  }

  // Use mutex to safely update the global sum
  pthread_mutex_lock(&mutex);
  data->result += local_sum;
  pthread_mutex_unlock(&mutex);

  pthread_exit(NULL);
}

// Threaded vectorized distance computation
void *thread_func_vec(void *arg) {
  ThreadData *data = (ThreadData *)arg;
  __m256 sum_vec = _mm256_setzero_ps();

  for (int i = data->start; i < data->end; i += 8) {
    __m256 u = _mm256_load_ps(&data->U[i]);
    __m256 v = _mm256_load_ps(&data->V[i]);

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

  float temp[8] __attribute__((aligned(32)));
  _mm256_store_ps(temp, sum_vec);
  double local_sum = 0.0;
  for (int i = 0; i < 8; i++) {
    local_sum += temp[i];
  }

  // Use mutex to safely update the global sum
  pthread_mutex_lock(&mutex);
  data->result += local_sum;
  pthread_mutex_unlock(&mutex);

  pthread_exit(NULL);
}

// Parallel computation using pthreads
void distPar(float *U, float *V, int n, int nb_threads, int mode,
             double *result) {
  pthread_t threads[nb_threads];
  ThreadData thread_data[nb_threads];
  int chunk_size = n / nb_threads;
  *result = 0.0;

  for (int i = 0; i < nb_threads; i++) {
    thread_data[i].U = U;
    thread_data[i].V = V;
    thread_data[i].start = i * chunk_size;
    thread_data[i].end = (i == nb_threads - 1) ? n : (i + 1) * chunk_size;
    thread_data[i].result = 0.0;

    if (mode == 0) {
      pthread_create(&threads[i], NULL, thread_func, &thread_data[i]);
    } else {
      pthread_create(&threads[i], NULL, thread_func_vec, &thread_data[i]);
    }
  }

  for (int i = 0; i < nb_threads; i++) {
    pthread_join(threads[i], NULL);
    *result += thread_data[i].result;
  }
}

int main() {
  float *U = (float *)aligned_alloc(32, N * sizeof(float));
  float *V = (float *)aligned_alloc(32, N * sizeof(float));

  srand(time(NULL));
  for (int i = 0; i < N; i++) {
    U[i] = (float)rand() / RAND_MAX;
    V[i] = (float)rand() / RAND_MAX;
  }

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
