// my_functions.h
#include <string>
#ifndef FUNCTIONS_H
#define FUNCTIONS_H
#define SAMPLE_MAX 48000
#define DIMENSION 1
extern double series[DIMENSION][SAMPLE_MAX];

// Deklarasi fungsi
std::string native_lib_main(double data[DIMENSION][SAMPLE_MAX]);
void get_data(double data[DIMENSION][SAMPLE_MAX]);

#endif // FUNCTIONS_H
