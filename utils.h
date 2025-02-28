// utils.h
#ifndef UTILS_H
#define UTILS_H

#include <immintrin.h>

// Function to allocate aligned memory and initialize with random values
void initialize_random(float *arr, int size);

//  Function to copy and cast an array of floats to an array of doubles
void copy_cast_array(float *src_arr, double *target_arr, int size);

// Function to print the first N elements of an array
void print_array(float *arr, int size);

// Function to get the current time in seconds
double now();

#endif // UTILS_H
