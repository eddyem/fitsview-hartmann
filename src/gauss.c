//      gauss.c - funtions for gaussian approximation
//
//      Copyright 2011 Edward V. Emelianoff <eddy@sao.ru>
//
//      This program is free software; you can redistribute it and/or modify
//      it under the terms of the GNU General Public License as published by
//      the Free Software Foundation; either version 2 of the License, or
//      (at your option) any later version.
//
//      This program is distributed in the hope that it will be useful,
//      but WITHOUT ANY WARRANTY; without even the implied warranty of
//      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//      GNU General Public License for more details.
//
//      You should have received a copy of the GNU General Public License
//      along with this program; if not, write to the Free Software
//      Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
//      MA 02110-1301, USA.
#include "gtk.h"
#include "gauss.h"
#ifdef GSL_FOUND // cmake found GSL library
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>


#undef DBG
#undef EBUG
#define DBG(...)


/*
 * Returns value of gaussian value in point x
 * x0, C, A, sigma - are parameters (middle, vertical shift, amplitude and std)
 */
double gaussian_pt(double C, double A, double sigma, double x0, double x){
	double s2 = sigma*sigma;
	double ss = A/sqrt(2*M_PI*s2);
	double E = x-x0;
	return C + ss*exp(-E*E/2/s2);
}

/*
 * Fills the ordinates of structure Points by Gaussian
 * (Abscissa does not move, so that this structure can be used
 * for multiple rendering of the same area but with different
 * parameters)
 */
void gaussian_v(double C, double A, double sigma, double x0, Points *pts){
	int i, n = pts->n;
	Point *pt = pts->data;
	double s2 = sigma*sigma * 2;
	double ss = A/sqrt(2*M_PI*s2);
	double E;
	for(i = 0; i < n; i++, pt++){
		E = pt->x - x0;
		pt->y = C + ss*exp(-E*E/s2);
		//DBG("pt %d (%g, %g)", i, pt->x, pt->y);
	}
}

int gauss_f(const gsl_vector *parms, void *data, gsl_vector *f){
	size_t n = ((struct data *)data)->n;
	double *x = ((struct data *)data)->x;
	double *y = ((struct data *)data)->y;
	double *dy = ((struct data *)data)->dy;
	double x0 = gsl_vector_get(parms, 3);
	double sigma = gsl_vector_get(parms, 2);
	double A = gsl_vector_get(parms, 1) * 0.65; // without 0.65 A would be wrong (WTF?)
	double C = gsl_vector_get(parms, 0);
	size_t i;
	double mul = A/fabs(sigma)/sqrt(2.*M_PI);
	for (i = 0; i < n; i++){
		/* Model Yi = C + A/sqrt(2*pi*sigma^2) * exp(-(x-b)^2/2/sigma^2) */
		double E = (x[i]-x0)/sigma;
		double Yi = C + mul * exp(-E * E/2.);
		gsl_vector_set(f, i, (Yi - y[i])/dy[i]);
	}
	return GSL_SUCCESS;
}

int gauss_df(const gsl_vector *parms, void *data, gsl_matrix *J){
	size_t n = ((struct data *)data)->n;
	double *x = ((struct data *)data)->x;
	double *dy = ((struct data *) data)->dy;
	double A = gsl_vector_get(parms, 1);
	double sigma = gsl_vector_get(parms, 2);
	double x0 = gsl_vector_get(parms, 3);
	double S2PI = 1. / sqrt(2.*M_PI);
	size_t i;
	for(i = 0; i < n; i++){
		/* Jacobian matrix J(i,j) = dfi / dxj, */
		/* and the xj are the parameters (C, A,sigma,x0) */
		double s = dy[i];
		double xbs = (x[i]-x0)/sigma;
		double y = exp(-xbs*xbs/2.) *S2PI / fabs(sigma) / s;
		double Ays = A*y/sigma;
		// df/dx0 = f * (x-x0)/sigma^2
		gsl_matrix_set(J, i, 3, Ays*xbs);
		 // df/dsigma = -f(1-(x-x0)^2/sigma^2)/sigma
		gsl_matrix_set(J, i, 2, Ays*(xbs*xbs-1.));
		gsl_matrix_set(J, i, 1, y); // df/dA
		gsl_matrix_set(J, i, 0, 1. / s); // df/dC
	}
	return GSL_SUCCESS;
}

int gauss_fdf(const gsl_vector *parms, void *data, gsl_vector *f, gsl_matrix *J){
	gauss_f(parms, data, f);
	gauss_df(parms, data, J);
	return GSL_SUCCESS;
}

// pseudorandom numbers generator initialisation
void urandom_ini(){
	double tt = dtime() * 1e6;
	double mx = (double)LONG_MAX;
	srand48((long)(tt - mx * floor(tt/mx)));
}


#define FIT(i) gsl_vector_get(s->x, i)
#define ERR(i) sqrt(gsl_matrix_get(covar,i,i))
#define FN(i) gsl_vector_get(s->f, i)
/*
 * Gaussian parameters calculation y=A/sqrt(2*pi*sigma^2) exp(-(x-x_0)^2/2/sigma^2),
 * which approximates the points set pts
 * Parameters A_, sigma_, x0_ may be NULL (if you don't need any of them)
 */
void gauss_fit(Points *pts, double *C_, double *A_, double *sigma_, double *x0_){
	// VVVV lower parameters may be formed as a structure to change as function argument
	double
		epsabs   = 1e-8,// absolute error
		epsrel   = 1e-5,// relative error
		chi_max  = 0.01;// max chi value for iterations criteria
	int max_iter = 300; // limit iterations number of gsl_multifit_fdfsolver
	size_t N_MIN = 10;	// minimum points for approximation
	double x_init[4];
	// AAAA upper parameters may be formed as a structure to change as function argument
/* x_init, the best approximations:
 * x0 - not far from real (the nearest is the better)
 * sigma - not far from real (the nearest is the better)
 * A - not large ~10 (it has a weak effect)
 */
	const gsl_multifit_fdfsolver_type *T;
	gsl_multifit_fdfsolver *s;
	int status;
	#ifdef EBUG
	int appNo = 0;
	#endif
	int iter;
	size_t i, j, n = pts->n, oldn;
	const size_t p = 4;
	gsl_matrix *covar = gsl_matrix_alloc (p, p);
	#ifdef EBUG
	double t0;
	#endif
	double *x, *y, *dy, chi, C, A, sigma, x0;
	if(n < 1) return;
	x = malloc(n * sizeof(double));
	y = malloc(n * sizeof(double));
	dy = malloc(n * sizeof(double));
	struct data d = {n, x, y, dy};
	gsl_multifit_function_fdf f;
	gsl_vector_view xx = gsl_vector_view_array(x_init, p);
	const gsl_rng_type *type;
	gsl_rng *r;

	gsl_rng_env_setup();
	type = gsl_rng_default;
	r = gsl_rng_alloc (type);
	f.f = &gauss_f;
	f.df = &gauss_df;
	f.fdf = &gauss_fdf;
	f.n = n;
	f.p = p;
	f.params = &d;
	// fill data structure. Don't forget Okkam's razor!!!
	{
		Point *pt = pts->data;
		double *px = x, *py = y, *pdy = dy, sum = 0.;
		for(i = 0; i < n; i++, pt++){
			*pdy++ = 1.; // I have no idea what is it, so init by 1
			*px++  = pt->x;
			*py++  = pt->y;
			sum   += pt->y;
			//DBG("point %d: (%g, %g)", i, pt->x, pt->y);
		}
		// fill x_init: x0, sigma, C, A (it can be a funtion parameter)
		x_init[3] = (*(--px) + *x) / 2.;
		x_init[2] = fabs((*x - *px) / 4.);
		x_init[0] = sum/(double)n;
		x_init[1] = sum;
		DBG("\nInitial parameters: x0=%.1f, sigma=%.1f, A=%.1f, C=%.1f",
			x_init[3], x_init[2], x_init[1], x_init[0]);
	}
	T = gsl_multifit_fdfsolver_lmder; // or also gsl_multifit_fdfsolver_lmsder
	s = gsl_multifit_fdfsolver_alloc(T, n, p);
	#ifdef EBUG
	t0 = dtime();
	#endif
	do{
		double dof, tres, c;
		DBG("\n************ Approximation %d ******************\n", appNo++);
		iter = 0;
		gsl_multifit_fdfsolver_set(s, &f, &xx.vector);
		do{
			iter++;
			status = gsl_multifit_fdfsolver_iterate(s);
			if(status)
				break;
			status = gsl_multifit_test_delta(s->dx, s->x, epsabs, epsrel);
		}while(status == GSL_CONTINUE && iter < max_iter);
		DBG("time=%g\n", dtime()-t0);
		gsl_multifit_covar(s->J, 0.0, covar);
		chi = gsl_blas_dnrm2(s->f);
		dof = n - p;
		tres = chi;
		c = chi / sqrt(dof); // GSL_MAX_DBL(1., chi / sqrt(dof));
		C = FIT(0), A = FIT(1), sigma = FIT(2), x0 = FIT(3);
		DBG("Number of iteratons = %d\n", iter);
		DBG("chi = %g, chi/dof = %g\n", chi, chi / sqrt(dof));
		DBG("C      = %.5f +/- %.5f\n", C, c*ERR(0));
		DBG("A      = %.5f +/- %.5f\n", A, c*ERR(1));
		DBG("sigma = %.5f +/- %.5f\n", sigma, c*ERR(2));
		DBG("x0     = %.5f +/- %.5f\n", x0, c*ERR(3));
		j = 0;
		oldn = n;
		if(c < chi_max) break;
		// throw out bad (by chi) data
		for(i = 0; i < n; i++){
			if(fabs(FN(i)) < tres){
				if(i != j){
					x[j] = x[i];
					y[j] = y[i];
					dy[j] = dy[i];
				}
				j++; continue;
			}
		}
		if(j != n){
			DBG("Chi tresholding %g, %zd points of %zd\n", tres, j, n);
			n = j;
			d.n = n;
		}
	}while(chi > chi_max && n != oldn && n > N_MIN);
	if(C_) *C_ = C;
	if(A_) *A_ = A;
	if(sigma_) *sigma_ = sigma;
	if(x0_) *x0_ = x0;
	//printf ("status = %s\n", gsl_strerror (status));
	gsl_multifit_fdfsolver_free(s);
	gsl_matrix_free(covar);
	gsl_rng_free(r);
	free(x); free(y); free(dy);
}
/*
int circle_f(const gsl_vector *parms, void *data, gsl_vector *f){
	size_t n = ((struct data *)data)->n;
	double *x = ((struct data *)data)->x;
	double *y = ((struct data *)data)->y;
	//double *dy = ((struct data *)data)->dy;
	double X = gsl_vector_get(parms, 0);
	double Y = gsl_vector_get(parms, 1);
	double R = gsl_vector_get(parms, 2);
	double sign, Yi, DX;
	size_t i;
	for (i = 0; i < n; i++){
		sign = (y[i] > Y) ? 1. : -1.;
		DX = x[i] - X;
		// Model Yi = Yc +- sqrt[R^2-(x-Xc)^2]
		Yi = Y + sign * sqrt(R*R - DX*DX);
		gsl_vector_set(f, i, (Yi - y[i]));
	}
	return GSL_SUCCESS;
}

int circle_df(const gsl_vector *parms, void *data, gsl_matrix *J){
	size_t n = ((struct data *)data)->n;
	double *x = ((struct data *)data)->x;
	double *y = ((struct data *)data)->y;
	//double *dy = ((struct data *)data)->dy;
	double X = gsl_vector_get(parms, 0);
	double Y = gsl_vector_get(parms, 1);
	double R = gsl_vector_get(parms, 2);
	double sign, Yi, DX;
	size_t i;
	for(i = 0; i < n; i++){
		// Jacobian matrix J(i,j) = dfi / dxj,
		// and the xj are the parameters (Xc, Yc, R)
		DX = x[i] - X;
		sign = (y[i] > Y) ? 1. : -1.;
		Yi = sign / sqrt(R*R - DX*DX);
		gsl_matrix_set(J, i, 2, R*Yi); // df/dR
		gsl_matrix_set(J, i, 1, 1.); // df/dYc
		gsl_matrix_set(J, i, 0, DX * Yi); // df/dXc
	}
	return GSL_SUCCESS;
}

int circle_fdf(const gsl_vector *parms, void *data, gsl_vector *f, gsl_matrix *J){
	circle_f(parms, data, f);
	circle_df(parms, data, J);
	return GSL_SUCCESS;
}

void circle_fit(struct data *d, double *Xcenter,
				double *Ycenter, double *Radius){
	double
		epsabs   = 1e-8,
		epsrel   = 1e-5;
	int max_iter = 30;
	double x_init[3], c, chi;
	const gsl_multifit_fdfsolver_type *T;
	gsl_multifit_fdfsolver *s;
	int status;
	int iter = 0;
	size_t n = d->n;
	const size_t p = 3;
	FNAME();
	gsl_matrix *covar = gsl_matrix_alloc (p, p);
	gsl_multifit_function_fdf f;
	gsl_vector_view xx = gsl_vector_view_array(x_init, p);
	const gsl_rng_type *type;
	gsl_rng *r;
	gsl_rng_env_setup();
	type = gsl_rng_default;
	r = gsl_rng_alloc(type);
	f.f = &circle_f;
	f.df = &circle_df;
	f.fdf = &circle_fdf;
	f.n = n;
	f.p = p;
	f.params = d;
	T = gsl_multifit_fdfsolver_lmder;
	s = gsl_multifit_fdfsolver_alloc(T, n, p);
	gsl_multifit_fdfsolver_set(s, &f, &xx.vector);
	x_init[0] = *Xcenter;
	x_init[1] = *Ycenter;
	x_init[2] = *Radius;
	do{
		iter++;
		status = gsl_multifit_fdfsolver_iterate(s);
		if(status) break;
		DBG("iter %d", iter);
		status = gsl_multifit_test_delta(s->dx, s->x, epsabs, epsrel);
	}while(status == GSL_CONTINUE && iter < max_iter);
	DBG("iter end");
	gsl_multifit_covar(s->J, 0.0, covar);
	chi = gsl_blas_dnrm2(s->f);
	c = chi / sqrt((double)(n - p));
	*Xcenter = FIT(0);
	*Ycenter = FIT(1);
	*Radius  = FIT(2);
	DBG("Number of iteratons = %d\n", iter);
	DBG("chi = %g, chi/dof = %g\n", chi, c);
	DBG("Xc = %.5f +/- %.5f\n", *Xcenter, c*ERR(0));
	DBG("Yc = %.5f +/- %.5f\n", *Ycenter, c*ERR(1));
	DBG("R  = %.5f +/- %.5f\n", *Radius,  c*ERR(2));
	gsl_multifit_fdfsolver_free(s);
	gsl_matrix_free(covar);
	gsl_rng_free(r);
}
*/

#else // GSL not found

//double gaussian_pt(double C, double A, double sigma, double x0, double x){
//	GSLERR; return A+C+sigma+x+x0;}
void gaussian_v(double C __attribute__((unused)), double A  __attribute__((unused)),
				double sigma __attribute__((unused)), double x0 __attribute__((unused)),
				Points *pts  __attribute__((unused))){GSLERR;}
void gauss_fit(Points *pts __attribute__((unused)), double *C_ __attribute__((unused)),
				double *A_ __attribute__((unused)), double *sigma_ __attribute__((unused)),
				double *x0_ __attribute__((unused))){	GSLERR;}
#endif // GSL_FOUND
