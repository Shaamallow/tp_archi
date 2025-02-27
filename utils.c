// utils.c
#include "utils.h"
#include <stdio.h>
#include <stdlib.h>

// Function to initialize an array with random float values between 0 and 1
void initialize_random(float *arr, int size) {
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
