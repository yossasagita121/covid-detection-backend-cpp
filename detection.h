//
// Created by Armein Z R Langi on 11/24/2023.
//

#ifndef BATUK_DETECTION_H
#define BATUK_DETECTION_H

#include <string>
#include <string>
// #include <cstdio>
#include <cstdlib>
#include <climits>
#include <ctime>
#include <cmath>

#define SAMPLE_MAX 48000
/* output is written every WHEN seconds */
#define WHEN 120
/* Size of the field for box assisted neighbour searching
   (has to be a power of 2)*/
#define NMAX 256
/* Size of the box for the scramble routine */
#define SCBOX 4096

#define DIMENSION 1
#define EMBEDDING 10
#define DIM 1
#define EMBED 10
#define HOWOFTEN 100
#define DELAY 1
#define SIZE_OF_LONG 4

extern double series[DIMENSION][SAMPLE_MAX];
typedef struct
{
    float dim = 0;
    float size = 0;
    float dispersi = 0;

} singularity;

int imax = NMAX - 1, howoften1, imin, validscale = 1;
unsigned long length = SAMPLE_MAX;
long box[NMAX][NMAX], boxc1[NMAX];
long list[SAMPLE_MAX];
long listc1[SAMPLE_MAX];
long scr[SAMPLE_MAX - (EMBED - 1) * DELAY];
long oscr[SAMPLE_MAX - (EMBED - 1) * DELAY];
// long scr[SAMPLE_MAX - (EMBED - 1) * DELAY];
double found[DIM * EMBED][HOWOFTEN];
double norm[HOWOFTEN];
double epsm[HOWOFTEN];
double EPSMAX = 1.0, EPSMIN = 1.e-3;
double epsfactor, epsinv, lneps, lnfac;
long nmax;
unsigned long MINDIST = 0, MAXFOUND = 1000;
double localdimension[DIMENSION * EMBEDDING], localsize[DIMENSION * EMBEDDING];
unsigned long get_input(char *filename);
void scramble(void);
// int proses(void);
void make_c2_1(int n);
void make_c2_dim(int n);
void rescale_data(double *x, unsigned long l, double *min, double *interval);
void get_data(void);
int proses(void);
float measure(singularity &s);
// void show_options(char *progname);
#endif // DETECTION_H

#ifndef _TISEAN_CEC_H
#define _TISEAN_CEC_H

/* These are the codes for the routines subtree */
#define RESCALE_DATA_ZERO_INTERVAL 11
#define CHECK_ALLOC_NOT_ENOUGH_MEMORY 12
#define CHECK_OPTION_NOT_UNSIGNED 13
#define CHECK_OPTION_NOT_INTEGER 14
#define CHECK_OPTION_NOT_FLOAT 15
#define CHECK_OPTION_NOT_TWO 16
#define CHECK_OPTION_C_NO_VALUE 17
#define TEST_OUTFILE_NO_WRITE_ACCESS 18
#define SOLVELE_SINGULAR_MATRIX 19
#define GET_SERIES_NO_LINES 20
#define GET_MULTI_SERIES_WRONG_TYPE_OF_C 21
#define GET_MULTI_SERIES_NO_LINES 22
#define VARIANCE_VAR_EQ_ZERO 23
#define EIG2_TOO_MANY_ITERATIONS 24
#define CHECK_OPTION_NOT_THREE 25

/* These are the codes for the main routines */
#define LYAP_SPEC_NOT_ENOUGH_NEIGHBORS 50
#define LYAP_SPEC_DATA_TOO_SHORT 51
#define AR_MODEL_TOO_MANY_POLES 52
#define EXTREMA_STRANGE_COMPONENT 53
#define FALSE_NEAREST_NOT_ENOUGH_POINTS 54
#define FSLE__TOO_LARGE_MINEPS 55
#define GHKSS__TOO_MANY_NEIGHBORS 56
#define NSTAT_Z__INVALID_STRING_FOR_OPTION 57
#define NSTAT_Z__NOT_UNSIGNED_FOR_OPTION 58
#define NSTAT_Z__TOO_LARGE_FOR_OPTION 59
#define NSTAT_Z__OPTION_NOT_SET 60
#define NSTAT_Z__TOO_MANY_PIECES 61
#define NSTEP__ESCAPE_REGION 62
#define POINCARE__WRONG_COMPONENT 63
#define POINCARE__OUTSIDE_REGION 64
#define POLYBACK__WRONG_PARAMETER_FILE 65
#define POLYNOMP__WRONG_PARAMETER_FILE 66
#define RESCALE__WRONG_INTERVAL 67
#define SAV_GOL__UNDERDETERMINED 68
#define SAV_GOL__TOO_LARGE_DERIVATIVE 69
#define MAKENOISE__FLAGS_REQUIRED 70
#define ZEROTH__STEP_TOO_LARGE 71
#define LYAP_K__MAXITER_TOO_LARGE 72
#define DELAY_WRONG_FORMAT_F 73
#define DELAY_DIM_NOT_EQUAL_F_M 74
#define DELAY_DIM_NOT_EQUAL_F_m 75
#define DELAY_WRONG_FORMAT_D 76
#define DELAY_WRONG_NUM_D 77
#define DELAY_INCONS_d_D 78
#define DELAY_SMALL_ZERO 79
#define DELAY_INCONS_m_M 80
#define ONESTEP_TOO_FEW_POINTS 81
#define MEM_SPEC_TOO_MANY_POLES 82

/* Global stuff */
#define VECTOR_TOO_LARGE_FOR_LENGTH 100

#endif // BATUK_DETECTION_H
