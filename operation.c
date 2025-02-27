#include <immintrin.h>
#include <math.h>

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
