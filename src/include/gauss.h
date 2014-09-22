#ifndef _GAUSS_H_
#define _GAUSS_H_
#include "fitsview.h"
#include "tracking.h"
#ifndef GSL_FOUND
	#define GSLERR g_err(_("Install GSL library and remake sources for this functions"))
#endif
struct data {
	size_t n;
	double *x;
	double *y;
	double *dy;
};

double gaussian_pt(double C, double A, double sigma, double x0, double x);
void gaussian_v(double C, double A, double sigma, double x0, Points *pts);
void gauss_fit(Points *pts, double *C_, double *A_, double *sigma_, double *x0_);
/*void circle_fit(struct data *data, double *Xcenter,
				double *Ycenter, double *Radius);*/
#endif // _GAUSS_H_
