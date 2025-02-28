// utils.c
#include "utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>

// Function to initialize an array with random float values between 0 and 1
void initialize_random(float *arr, int size) {
  srand(time(NULL));
  for (int i = 0; i < size; i++) {
    arr[i] = (float)rand() / RAND_MAX;
  }
}

void copy_cast_array(float *src_arr, double *target_arr, int size) {
  for (int i = 0; i < size; i++) {
    target_arr[i] = (double)src_arr[i];
  }
}

// Function to print an array
void print_array(float *arr, int size) {
  for (int i = 0; i < size; i++) {
    printf("%f ", arr[i]);
    if ((i + 1) % 8 == 0)
      printf("\n");
  }
}

// Function to get the current time in seconds
double now() {
  struct timeval t;
  gettimeofday(&t, NULL);
  return t.tv_sec + t.tv_usec / 1000000.0;
}
