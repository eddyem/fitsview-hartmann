// CUtools.h - common definitions and functions for CUDA.c and NOCUDA.c
#ifndef _CUTOOLS_H_
#define _CUTOOLS_H_
#ifdef _CUDA_CU_  // file was included from CUDA.cu
	#define EXTERN extern "C"
#else
	#define EXTERN
#endif
#ifndef _GNU_SOURCE
	#define _GNU_SOURCE
#endif
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

typedef struct{
	unsigned char FilterType;	// filter type
	int w;						// filter width
	int h;						// height
	double sx;					// x half-width
	double sy;					// y half-width (sx, sy - for Gaussian-type filters)
} Filter;
// FilterType
enum{
	 LAPGAUSS			// laplasian of gaussian
	,GAUSS				// gaussian
	,SOBELH				// Sobel horizontal
	,SOBELV				// -//- vertical
	,SIMPLEGRAD			// simple gradient (by Sobel)
	,PREWITTH			// Prewitt (horizontal) - simple derivative
	,PREWITTV			// -//- (vertical)
	,MEDIAN				// median
	,STEP				// "posterisation"
};
// f->h for STEP filter
enum{
	 UNIFORM
	,LOG
	,EXP
	,SQRT
	,POW
};

// export section
EXTERN int fill_hough_lines(float *ima, float min, float max, int imW, int imH, int Rmax, int angles, float *hough);
EXTERN int DiffFilter(float *ima, float **result, Filter *f, int sizex, int sizey);
EXTERN int MedFilter(float *ima, float **result, Filter *f, int sizex, int sizey);
EXTERN int fillIsoScale(Filter *f, float **scale, float min, float wd);
EXTERN int StepFilter(float *ima, float **result, Filter *f, int sizex, int sizey, float min, float max, float **scale);
EXTERN int GradFilterSimple(float *ima, float **result, Filter *f, int sizex, int sizey);
#endif // _CUTOOLS_H_
