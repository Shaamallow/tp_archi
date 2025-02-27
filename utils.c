// utils.c
#include "utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

// Function to initialize an array with random float values between 0 and 1
void initialize_random(double *arr, int size) {
  for (int i = 0; i < size; i++) {
    arr[i] = (float)rand() / RAND_MAX;
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
