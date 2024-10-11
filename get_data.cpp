#include <iostream>
#include <algorithm>
#include "functions.h"

double series[DIMENSION][SAMPLE_MAX];

void get_data(double data[DIMENSION][SAMPLE_MAX])
{
    for (int i = 0; i < DIMENSION; ++i)
    {
        for (int j = 0; j < SAMPLE_MAX; ++j)
        {
            series[i][j] = data[i][j];
        }
    }
}