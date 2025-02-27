// utils.h
#ifndef UTILS_H
#define UTILS_H

#include <immintrin.h>

// Function to allocate aligned memory and initialize with random values
void initialize_random(double *arr, int size);

// Function to print the first N elements of an array
void print_array(float *arr, int size);

// Function to get the current time in seconds
double now();

#endif // UTILS_H
