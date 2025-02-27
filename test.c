#include <immintrin.h>
#include <math.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#define SIZE 16 // Must be a multiple of 8 for AVX

// Function to get the current time in seconds
double now() {
  struct timeval t;
  gettimeofday(&t, NULL);
  return t.tv_sec + t.tv_usec / 1000000.0;
}

double vectorial_compute(float *U, float *V) {
  __m256 sum_vec = _mm256_setzero_ps();

  for (int i = 0; i < SIZE; i += 8) {
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

  // Get the sum of the elements in the sum_vec

  // Get manual memory allocation aligned
  float temp[8] __attribute__((aligned(32)));
  _mm256_store_ps(temp, sum_vec);

  double sum = 0.0;
  for (int i = 0; i < 8; i++) {
    sum += temp[i];
  }

  return sum;
}

double sequential_compute(float *U, float *V) {
  double sum = 0.0;
  for (int i = 0; i < SIZE; i++) {
    sum +=
        sqrt((U[i] * U[i] + V[i] * V[i]) / (1 + (U[i] * V[i]) * (U[i] * V[i])));
  }
  return sum;
}

int main() {
  // Allocate aligned memory for AVX operations
  float *U = (float *)aligned_alloc(32, SIZE * sizeof(float));
  float *V = (float *)aligned_alloc(32, SIZE * sizeof(float));

  if (!U || !V) {
    printf("Memory allocation failed!\n");
    return 1;
  }

  // Initialize U and V with random float values between 0 and 1
  srand(time(NULL));
  for (int i = 0; i < SIZE; i++) {
    U[i] = (float)rand() / RAND_MAX;
    V[i] = (float)rand() / RAND_MAX;
  }

  double seq_sum = sequential_compute(U, V);
  double vec_sum = vectorial_compute(U, V);

  printf("Seq Sum: %f\n", seq_sum);
  printf("Vec Sum: %f\n", vec_sum);

  // Free allocated memory
  free(U);
  free(V);

  return 0;
}
