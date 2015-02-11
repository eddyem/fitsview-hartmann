//      NOCUDA.c - CPU-variants when there's no CUDA
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

#include "fitsview.h"
#include "gtk.h"
#include "CUtools.h"
#include <pthread.h>
#include <fftw3.h>

static float *Sin_d = NULL, *Cos_d = NULL;
static int sincosize = 0;

#ifndef THREAD_NUMBER
	#define THREAD_NUMBER 2		// default - 2 threads
#endif

typedef struct{
	float *image;
	float min;
	float max;
	float *hough;
	int w;
	int Y0;
	int Y1;
	int Rmax;
	int angles;
	pthread_mutex_t *mutex;
}Hough_kernel_args;

// sin/cos ini, [-90,+180)
void init_sincos(int angles){
	float anglestep = (float)angles / 270.; // step by an angle
	float conv = M_PI / 180.;
	float *cp, *sp;
	int k;
	if(!Sin_d || !Cos_d || angles != sincosize){ // first time initialize arrays
		free(Sin_d);
		free(Cos_d);
		Sin_d = malloc(angles * sizeof(float));
		Cos_d = malloc(angles * sizeof(float));
		cp = Cos_d; sp = Sin_d;
		for(k = 0; k < angles; k++){ // fill arrays from -90 to +90 degr in radians
			float theta = ((float)k / anglestep - 90.) * conv;
			*sp++ = sinf(theta);
			*cp++ = cosf(theta);
			//sincosf(theta, sp++, cp++);
		}
		sincosize = angles;
	}
}

// Line Hough transform kernel
void *fill_lin_hough_array(void *data){
	Hough_kernel_args *HD = (Hough_kernel_args*) data;
	float *ima = HD->image, *hough = HD->hough;
	float min = HD->min, wd = HD->max - min;
	pthread_mutex_t *mutex = HD->mutex;
	int Y0 = HD->Y0, Y1 = HD->Y1, imW = HD->w, Rmax = HD->Rmax, angles = HD->angles;
	int i, j, k, Y, R;
	for(j = Y0; j < Y1; j++){
		for(i = 0.; i < imW; i++, ima++){
			float imdata = (*ima - min) / wd;
			if(imdata > 0.1){
				Y = 0;
				pthread_mutex_lock(mutex);
				for(k = 0; k < angles; k++, Y+=Rmax){
					// R = x*cos(theta) + y*sin(theta)
					R = (int)(0.5 + i*Cos_d[k] + j*Sin_d[k]);
					if(R > 0 && R < Rmax) hough[R + Y] += imdata;
				}
				pthread_mutex_unlock(mutex);
			}
		}
	}
	return NULL;
}
/*
 * Build of a Hough transform to find lines
 * Input:
 *		ima - the image data
 *		min, max - range of the data in it
 *		imW, imH - image width and height
 *		Rmax - the maximum value of R
 *		angles - the array size of angles (the angle step is 180/angles degrees)
 * Output:
 *		hough - initialized by an external function array for Hough transform
 *		!!! array must be initialized with zeros before calling this function !!!
 * Output array is normalized to unity
 */
int fill_hough_lines(float *ima, float min, float max, int imW, int imH, int Rmax, int angles, float *hough){
	int i, Y0, Y1, dY;
	float *hptr = hough, hmax;
	Hough_kernel_args HD[THREAD_NUMBER];
	pthread_t threads[THREAD_NUMBER];
	pthread_mutex_t mutex;
	pthread_mutex_init(&mutex, NULL);
	init_sincos(angles);
	dY = (imH + THREAD_NUMBER/2) / THREAD_NUMBER;
	Y0 = 0; Y1 = dY;
	for(i = 0; i < THREAD_NUMBER; i++, Y0+=dY, Y1+=dY){
		if(Y0 >= imW) break;
		if(Y1 > imW) Y1 = imW;
		HD[i].min = min; HD[i].max = max;
		HD[i].image = ima; HD[i].hough = hough;
		HD[i].Y0 = Y0; HD[i].Y1 = Y1; HD[i].w = imW;
		HD[i].Rmax = Rmax; HD[i].angles = angles;
		HD[i].mutex = &mutex;
		pthread_create(&threads[i], NULL, fill_lin_hough_array, &HD[i]);
	}
	for(i = 0; i < THREAD_NUMBER; i++)
		pthread_join(threads[i], NULL);
	pthread_mutex_destroy(&mutex);
	hmax = *hough++;
	Y1 = Rmax * angles;
	for(i = 1; i < Y1; i++, hough++)
		if(*hough > hmax) hmax = *hough;
	hough = hptr;
	for(i = 0; i < Y1; i++, hough++)
		*hough /= hmax;
	return 1;
}




// from http://graphics.stanford.edu/%7Eseander/bithacks.html#RoundUpPowerOf2
int nextpow2(int i, int j){
	inline int p2oi(int i){
		unsigned int v = (unsigned int)i - 1;
		v |= v >> 1;
		v |= v >> 2;
		v |= v >> 4;
		v |= v >> 8;
		v |= v >> 16;
		v++;
		return (int) v;
	}
	int p1 = p2oi(i), p2 = p2oi(j);
	return MAX(p1, p2);
}

// FFT
fftw_plan fftimg, fftmask, ifft;
// FFT buffrers for image & filter
fftw_complex *Fmask = NULL,		// FFT of a filter
			*Fimg = NULL;		// picture FFT
static double	*mask = NULL,	// filter mask
				*img = NULL,	// picture
				*resm = NULL;

void fftshift(int size, double *m){
	int h = size/2, ss, ss1, i, j, k, l, p;
	double tmp;
	ss = (2*h+1)*h;
	ss1 = (2*h-1)*h;
	for(j = 0; j < h; j++){
		k = j * size;
		l = k + h;
		for(i = 0; i < h; i++, k++, l++){
			p = k + ss;
			tmp = m[k]; m[k] = m[p]; m[p] = tmp;
			p = l + ss1;
			tmp = m[l]; m[l] = m[p]; m[p] = tmp;
		}
	}
}

// Lapgauss mask building
void build_LG_filter(int size, Filter *f){
	int y, y0=0,y1=size, x0=0, x1=size;
	double sx2 = f->sx * f->sx, sy2 = f->sy * f->sy;
	double hw, hh;
	#ifdef EBUG
	double t0=dtime();
	#endif
	if(f->w < size && f->w > 0){
		x0 = (size - f->w + 1) / 2;
		x1 = x0 + f->w;
	}
	if(f->h < size && f->h > 0){
		y0 = (size - f->h + 1) / 2;
		y1 = y0 + f->h;
	}
	hh = -(double)size / 2.;
	hw = -(double)size / 2.+(double)x0;
	DBG("y0=%d, y1=%d, hw=%g, hh=%g",y0,y1,hw,hh);
	double ss = 3. / hh / hh / sqrt(-hh);
//	double ss = 1./sqrt(2*M_PI*f->sx*f->sy);
	#pragma omp parallel for
	for(y = y0; y < y1; y++){
		double X, Y, y2, x2, R;
		int str, x;
		X = hw;
		str = y * size;
		Y = ((double)y) + hh;
		y2 = Y*Y/sy2;
		for(x = x0; x < x1; x++, X+=1.){
			x2 = X*X/sx2;
			R = x2 + y2;
			mask[str + x] =ss * ((x2-1.)/sx2+(y2-1.)/sy2)*exp(-R/2.);
		}
	}
	DBG("time=%f\n", dtime()-t0);
}

// Gaussian mask building/
void build_G_filter(int size, Filter *f){
	int y, y0=0,y1=size, x0=0, x1=size;
	double sx2 = f->sx * f->sx, sy2 = f->sy * f->sy;
	double hw, hh;
	#ifdef EBUG
	double t0=dtime();
	#endif
	if(f->w < size && f->w > 0){
		x0 = (size - f->w + 1) / 2;
		x1 = x0 + f->w;
	}
	if(f->h < size && f->h > 0){
		y0 = (size - f->h + 1) / 2;
		y1 = y0 + f->h;
	}
	hh = -(double)size / 2.;
	hw = -(double)size / 2.+(double)x0;
	DBG("y0=%d, y1=%d, hw=%g, hh=%g",y0,y1,hw,hh);
	double ss = 1. / (M_PI*2.) / sx2 / sy2 / hh / hh;
	#pragma omp parallel for
	for(y = y0; y < y1; y++){
		double X, Y, y2, x2, R;
		int str, x;
		X = hw;
		str = y * size;
		Y = ((double)y) + hh;
		y2 = Y*Y/sy2;
		for(x = x0; x < x1; x++, X+=1.){
			x2 = X*X/sx2;
			R = x2 + y2;
			mask[str + x] =ss * exp(-R/2.);
		}
	}
	DBG("time=%f\n", dtime()-t0);
}

// Elementary filter mask building
void build_S_filter(int size, Filter *f){
	int y, a0, a1;
	double hh, Y, pt = 0.;
	a0 = (size - 2) / 2;
	a1 = a0 + 3;
	hh = -(double)(size / 2);
	double ss = 1. / (M_PI*2.) / hh / hh / sqrt(-hh);
	Y = -1.;
	for(y = a0; y < a1; y++, Y+=1.){
		double X = -1.;
		int str, x;
		str = y * size;
		for(x = a0; x < a1; x++, X+=1.){
			switch(f->FilterType){
				case SOBELH:
					pt = -X*(2.-fabs(Y));
				break;
				case SOBELV:
					pt = -Y*(2.-fabs(X));
				break;
				case PREWITTH:
					pt = X;
				break;
				case PREWITTV:
					pt = Y;
				break;
			/*	case :
					pt =
				break;
				case :
					pt =
				break;
				case :
					pt =
				break;
				case :
					pt =
				break;
				case :
					pt =
				break;*/
			}
			mask[str + x] = ss*pt;
		}
	}
	//mask[size*size/2]=0.;
}

/*
 * Filtering by convolution with a filter
 * Input:
 *		ima - a pointer to picture data, which should be filtered
 *		f - filter parameters
 * Output:
 *		result - allocated by this function memory area where
 *		the filtered image is placed
 * Returns TRUE, if the filter succeeded
 */
int DiffFilter(float *ima,
					float **result,
					Filter *f,
					int sizex,
					int sizey
					){
	int ssize, i, j, k, l;
	static int fftw_ini = 0;
	int size2;
	#ifdef EBUG
	double t0 = dtime();
	#endif
	size2 = nextpow2(sizex, sizey);
	ssize = size2 * size2; // FFT image size
	if(!fftw_ini)
		if(!(fftw_ini = fftw_init_threads())){
			//g_err(_("FFTW error");
			return FALSE;
		}
	void free_all(){
		_FREE(mask); _FREE(Fmask); _FREE(img);
		_FREE(Fimg); _FREE(resm);
	}
	free_all();
	DBG("img (%d x %d) -> (%d x %d), time=%f\n", sizex,sizey, size2,size2, dtime()-t0);
	// allocate memory for objects
	img = (double*)calloc(ssize, sizeof(double));
	Fimg = (fftw_complex*)calloc(ssize, sizeof(fftw_complex));
	if(!img || !Fimg){free_all(); return FALSE;}
	// define direct FFT
	fftw_plan_with_nthreads(THREAD_NUMBER);
	fftimg = fftw_plan_dft_r2c_2d(size2, size2, img, Fimg, FFTW_ESTIMATE);
	// copy ima -> img
	for(j = 0; j < sizey; j++){
		k = j * size2;
		l = j * sizex;
		for(i = 0; i < sizex; i++, k++, l++)
			img[k] = ima[l];
	}
	fftw_execute(fftimg);
	_FREE(img);
	// build filter
	DBG("build filter, time=%f\n", dtime()-t0);
	mask = (double*)calloc(ssize, sizeof(double));
	Fmask = (fftw_complex*)calloc(ssize, sizeof(fftw_complex));
	if(!mask || !Fmask){free_all(); return FALSE;}
	switch(f->FilterType){
		case LAPGAUSS:
			build_LG_filter(size2, f);
			break;
		case GAUSS:
			build_G_filter(size2, f);
			break;
		case SOBELH:
		case SOBELV:
		case PREWITTH:
		case PREWITTV:
			build_S_filter(size2, f);
			break;
	}
	// define filter FFT
	fftw_plan_with_nthreads(THREAD_NUMBER);
	fftmask = fftw_plan_dft_r2c_2d(size2, size2, mask, Fmask, FFTW_ESTIMATE);
	fftw_execute(fftmask);
	_FREE(mask);
	// filtered picture:
	DBG("filter image, time=%f\n", dtime()-t0);
	resm = (double*)calloc(ssize, sizeof(double));
	if(!resm){free_all(); return FALSE;}
	// define inverse FFT
	fftw_plan_with_nthreads(THREAD_NUMBER);
	ifft = fftw_plan_dft_c2r_2d(size2, size2, Fimg, resm, FFTW_ESTIMATE);
	// DON'T PARALLEL THIS, it will be slower
	for(i=0; i<ssize; i++){ // convolution by multiplication in Fourier space
	/**/
		double a,b,c,d;
		a=Fimg[i][0]; c=Fmask[i][0];
		b=Fimg[i][1]; d=Fmask[i][1];
		Fimg[i][0] = a*c - b*d;
		Fimg[i][1] = b*c + a*d;
	/**
		Fimg[i][0] = Fmask[i][0];
		Fimg[i][1] = Fmask[i][1];
	**/
	}
	_FREE(Fmask);
	fftw_execute(ifft);
	_FREE(Fimg);
	fftshift(size2, resm);
	*result = calloc(ssize, sizeof(float));
	if(!*result){free_all(); return FALSE;}
	float *tmp = *result;
	// copy inv FFT -> result
	for(j = 0; j < sizey; j++){
		k = j * size2;
		l = j * sizex;
		for(i = 0; i < sizex; i++, k++, l++)
			tmp[l] = resm[k];
	}
	_FREE(resm);
	fftw_destroy_plan(fftimg);
	fftw_destroy_plan(ifft);
	fftw_destroy_plan(fftmask);
	fftw_cleanup_threads();
	DBG("time=%f\n", dtime()-t0);
	return TRUE;
}

/*
 * Simple gradient filter based on two Sobel filters
 * output = sqrt(SobelH(input)^2+SobelV(input)^2)
 */
int GradFilterSimple(float *ima, float **result, Filter *f, int w, int h){
	float *dst1, *dst2;
	int res = FALSE, y;
	#ifdef EBUG
	double t0 = dtime();
	#endif
	f->FilterType = SOBELH;
	res = DiffFilter(ima, &dst1, f, w, h);
	if(!res) return FALSE;
	f->FilterType = SOBELV;
	res = DiffFilter(ima, &dst2, f, w, h);
	if(!res){free(dst1); return FALSE;}
	*result = dst1;
	#pragma omp parallel for
	for(y = 0; y < h; y++){
		int x;
		float *in, *out;
		in = dst2 + y*w;
		out = dst1 + y*w;
		for(x = 0; x < w; x++, in++, out++)
			*out = sqrtf((*in)*(*in) + (*out)*(*out));
	}
	free(dst2);
	DBG("time=%f\n", dtime()-t0);
	return TRUE;
}

/*
 * Quick median functions stolen from
 *
 *   Fast median search: an ANSI C implementation
 *    Nicolas Devillard - ndevilla AT free DOT fr
 *                  July 1998
 */
#define PIX_SORT(a,b) { if ((*a)>(*b)) ELEM_SWAP((a),(b)); }
#define ELEM_SWAP(a,b) { register float *t=(a);(a)=(b);(b)=t; }
float opt_med3(float **p, int n __attribute__((unused))){
	PIX_SORT(p[0],p[1]) ; PIX_SORT(p[1],p[2]) ; PIX_SORT(p[0],p[1]) ;
	return(*p[1]) ;
}float opt_med5(float **p, int n __attribute__((unused))){
	PIX_SORT(p[0],p[1]) ; PIX_SORT(p[3],p[4]) ; PIX_SORT(p[0],p[3]) ;
	PIX_SORT(p[1],p[4]) ; PIX_SORT(p[1],p[2]) ; PIX_SORT(p[2],p[3]) ;
	PIX_SORT(p[1],p[2]) ; return(*p[2]) ;
}float opt_med7(float **p, int n __attribute__((unused))){
	PIX_SORT(p[0], p[5]) ; PIX_SORT(p[0], p[3]) ; PIX_SORT(p[1], p[6]) ;
	PIX_SORT(p[2], p[4]) ; PIX_SORT(p[0], p[1]) ; PIX_SORT(p[3], p[5]) ;
	PIX_SORT(p[2], p[6]) ; PIX_SORT(p[2], p[3]) ; PIX_SORT(p[3], p[6]) ;
	PIX_SORT(p[4], p[5]) ; PIX_SORT(p[1], p[4]) ; PIX_SORT(p[1], p[3]) ;
	PIX_SORT(p[3], p[4]) ; return (*p[3]) ;
}float opt_med9(float **p, int n __attribute__((unused))){
	PIX_SORT(p[1], p[2]) ; PIX_SORT(p[4], p[5]) ; PIX_SORT(p[7], p[8]) ;
	PIX_SORT(p[0], p[1]) ; PIX_SORT(p[3], p[4]) ; PIX_SORT(p[6], p[7]) ;
	PIX_SORT(p[1], p[2]) ; PIX_SORT(p[4], p[5]) ; PIX_SORT(p[7], p[8]) ;
	PIX_SORT(p[0], p[3]) ; PIX_SORT(p[5], p[8]) ; PIX_SORT(p[4], p[7]) ;
	PIX_SORT(p[3], p[6]) ; PIX_SORT(p[1], p[4]) ; PIX_SORT(p[2], p[5]) ;
	PIX_SORT(p[4], p[7]) ; PIX_SORT(p[4], p[2]) ; PIX_SORT(p[6], p[4]) ;
	PIX_SORT(p[4], p[2]) ; return(*p[4]) ;
} float opt_med25(float **p, int n __attribute__((unused))){
	PIX_SORT(p[0], p[1])  ; PIX_SORT(p[3], p[4])  ; PIX_SORT(p[2], p[4]) ;
	PIX_SORT(p[2], p[3])  ; PIX_SORT(p[6], p[7])  ; PIX_SORT(p[5], p[7]) ;
	PIX_SORT(p[5], p[6])  ; PIX_SORT(p[9], p[10]) ; PIX_SORT(p[8], p[10]) ;
	PIX_SORT(p[8], p[9])  ; PIX_SORT(p[12], p[13]); PIX_SORT(p[11], p[13]) ;
	PIX_SORT(p[11], p[12]); PIX_SORT(p[15], p[16]); PIX_SORT(p[14], p[16]) ;
	PIX_SORT(p[14], p[15]); PIX_SORT(p[18], p[19]); PIX_SORT(p[17], p[19]) ;
	PIX_SORT(p[17], p[18]); PIX_SORT(p[21], p[22]); PIX_SORT(p[20], p[22]) ;
	PIX_SORT(p[20], p[21]); PIX_SORT(p[23], p[24]); PIX_SORT(p[2], p[5]) ;
	PIX_SORT(p[3], p[6])  ; PIX_SORT(p[0], p[6])  ; PIX_SORT(p[0], p[3]) ;
	PIX_SORT(p[4], p[7])  ; PIX_SORT(p[1], p[7])  ; PIX_SORT(p[1], p[4]) ;
	PIX_SORT(p[11], p[14]); PIX_SORT(p[8], p[14]) ; PIX_SORT(p[8], p[11]) ;
	PIX_SORT(p[12], p[15]); PIX_SORT(p[9], p[15]) ; PIX_SORT(p[9], p[12]) ;
	PIX_SORT(p[13], p[16]); PIX_SORT(p[10], p[16]); PIX_SORT(p[10], p[13]) ;
	PIX_SORT(p[20], p[23]); PIX_SORT(p[17], p[23]); PIX_SORT(p[17], p[20]) ;
	PIX_SORT(p[21], p[24]); PIX_SORT(p[18], p[24]); PIX_SORT(p[18], p[21]) ;
	PIX_SORT(p[19], p[22]); PIX_SORT(p[8], p[17]) ; PIX_SORT(p[9], p[18]) ;
	PIX_SORT(p[0], p[18]) ; PIX_SORT(p[0], p[9])  ; PIX_SORT(p[10], p[19]) ;
	PIX_SORT(p[1], p[19]) ; PIX_SORT(p[1], p[10]) ; PIX_SORT(p[11], p[20]) ;
	PIX_SORT(p[2], p[20]) ; PIX_SORT(p[2], p[11]) ; PIX_SORT(p[12], p[21]) ;
	PIX_SORT(p[3], p[21]) ; PIX_SORT(p[3], p[12]) ; PIX_SORT(p[13], p[22]) ;
	PIX_SORT(p[4], p[22]) ; PIX_SORT(p[4], p[13]) ; PIX_SORT(p[14], p[23]) ;
	PIX_SORT(p[5], p[23]) ; PIX_SORT(p[5], p[14]) ; PIX_SORT(p[15], p[24]) ;
	PIX_SORT(p[6], p[24]) ; PIX_SORT(p[6], p[15]) ; PIX_SORT(p[7], p[16]) ;
	PIX_SORT(p[7], p[19]) ; PIX_SORT(p[13], p[21]); PIX_SORT(p[15], p[23]) ;
	PIX_SORT(p[7], p[13]) ; PIX_SORT(p[7], p[15]) ; PIX_SORT(p[1], p[9]) ;
	PIX_SORT(p[3], p[11]) ; PIX_SORT(p[5], p[17]) ; PIX_SORT(p[11], p[17]) ;
	PIX_SORT(p[9], p[17]) ; PIX_SORT(p[4], p[10]) ; PIX_SORT(p[6], p[12]) ;
	PIX_SORT(p[7], p[14]) ; PIX_SORT(p[4], p[6])  ; PIX_SORT(p[4], p[7]) ;
	PIX_SORT(p[12], p[14]); PIX_SORT(p[10], p[14]); PIX_SORT(p[6], p[7]) ;
	PIX_SORT(p[10], p[12]); PIX_SORT(p[6], p[10]) ; PIX_SORT(p[6], p[17]) ;
	PIX_SORT(p[12], p[17]); PIX_SORT(p[7], p[17]) ; PIX_SORT(p[7], p[10]) ;
	PIX_SORT(p[12], p[18]); PIX_SORT(p[7], p[12]) ; PIX_SORT(p[10], p[18]) ;
	PIX_SORT(p[12], p[20]); PIX_SORT(p[10], p[20]); PIX_SORT(p[10], p[12]) ;
	return (*p[12]);
}
float quick_select(float **arr, int n){
	int low, high;
	int median;
	int middle, ll, hh;
	float ret;
	low = 0 ; high = n-1 ; median = (low + high) / 2;
	for(;;){
		if(high <= low) /* One element only */
			break;
		if(high == low + 1){ /* Two elements only */
			PIX_SORT(arr[low], arr[high]) ;
			break;
		}
		/* Find median of low, middle and high items; swap into position low */
		middle = (low + high) / 2;
		PIX_SORT(arr[middle], arr[high]) ;
		PIX_SORT(arr[low], arr[high]) ;
		PIX_SORT(arr[middle], arr[low]) ;
		/* Swap low item (now in position middle) into position (low+1) */
		ELEM_SWAP(arr[middle], arr[low+1]) ;
		/* Nibble from each end towards middle, swapping items when stuck */
		ll = low + 1;
		hh = high;
		for(;;){
			do ll++; while (*arr[low] > *arr[ll]);
			do hh--; while (*arr[hh] > *arr[low]);
			if(hh < ll) break;
			ELEM_SWAP(arr[ll], arr[hh]) ;
		}
		/* Swap middle item (in position low) back into correct position */
		ELEM_SWAP(arr[low], arr[hh]) ;
		/* Re-set active partition */
		if (hh <= median) low = ll;
		if (hh >= median) high = hh - 1;
	}
	ret = *arr[median];
	return ret;
}
#undef PIX_SORT
#undef ELEM_SWAP

/*
#define _0(x)	(x & 0x7FF)			// lower 11 bits
#define _1(x)	(x >> 11 & 0x7FF)	// middle 11 bits
#define _2(x)	(x >> 22 )			// upper 10 bits
float hist_select(float **arr, int n){
	inline uint32_t FloatFlip(uint32_t f){
		uint32_t mask = -((int32_t)(f >> 31)) | 0x80000000;
		return f ^ mask;
	}
	inline uint32_t IFloatFlip(uint32_t f){
		uint32_t mask = ((f >> 31) - 1) | 0x80000000;
		return f ^ mask;
	}
	uint32_t b0
}
#undef _0
#undef _1
#undef _2
*/

/*
 * Median filtering
 * Input:
 * 		ima - picture (free() should be in a caller)
 * 		f - filter
 * 		sizex, sizey - picture size
 * Output:
 * 		result - filtered picture, memory allocates in this routine.
 * 				Image borders that didn't pass through the filter (+- 1/2 of filter size)
 * 				filled zeros
 */
int MedFilter(float *ima,
					float **result,
					Filter *f,
					int sizex,
					int sizey
					){
	int xlow = f->w / 2, xhigh = f->w - xlow;// filter area borders [x0-xlow, x0+xhigh)
	int ylow = f->h / 2, yhigh = f->h - ylow; // [y0-ylow, y0+yhigh)
	int x, xm=sizex-xhigh, ym=sizey-yhigh; // xm,ym - upper boundaries if filtered picture
	int ssize=sizex*sizey, fsz=f->w*f->h, H = f->h - 1;
	#ifdef EBUG
	double t0 = dtime();
	#endif
	float (*medfn)(float **p, int n);
	*result = calloc(ssize, sizeof(float));
	if(!result) return FALSE;
	switch(f->w*f->h){
		case 3:
			medfn = opt_med3;
		break;
		case 5:
			medfn = opt_med5;
		break;
		case 7:
			medfn = opt_med7;
		break;
		case 9:
			medfn = opt_med9;
		break;
		case 25:
			medfn = opt_med25;
		break;
		default:
			medfn = quick_select;
		break;
	}
	/*
	 * Selection of picture elements on the square wxh into an array arr
	 * x0, y0 - coordinates of the center of the square
	 * cntr - current replacement string, or -1 if it is necessary to fill the entire arr
	 * If cntr!=-1 next line replaces a string cntr
	 */
	float selarr0(int x0, int y0, float *arr, float **sel){
		int x,y,tx,ty;
		float *tmp, *ptr;
		tx=x0+xhigh; ty=y0+yhigh;
		tmp = arr;
		for(y=y0-ylow; y<ty; y++){ // completely fill the array
			ptr = &ima[y*sizex+x0-xlow]; // FILL BY LINES!!!
			for(x=x0-xlow; x<tx; x++)
				*tmp++ = *ptr++;
		}
		return medfn(sel, fsz);
	}
	float selarr1(int x0, int y0, float *arr, float **sel, int *cntr){
		int x,tx,ty;
		float *tmp, *ptr;
		tx=x0+xhigh; ty=y0+yhigh;
		tmp = &arr[*cntr*f->w]; // pointer to a changed line
		ptr = &ima[(ty-1)*sizex+x0-xlow];
		for(x=x0-xlow; x<tx; x++) // change data in column
			*tmp++ = *ptr++;
		if(++(*cntr) > H) *cntr = 0;
		return medfn(sel, fsz);
	}
	#pragma omp parallel
	{
		float *arr = calloc(fsz, sizeof(float*)); // array to sample storage
		float **sel = calloc(fsz, sizeof(float*));// array to sorted sample storage
		for(x = 0; x < fsz; x++) sel[x] = arr + x;
		if(arr){
		#pragma omp for
		for(x = xlow; x < xm; x++){
			int y, cntr = 0;
			(*result)[ylow*sizex + x] = selarr0(x,ylow, arr, sel);
			for(y = ylow+1; y < ym; y++){
				(*result)[y*sizex + x] = selarr1(x,y, arr, sel, &cntr);
			}
		}
		free(arr);
		free(sel);
		}
	}
	DBG("time=%f\n", dtime()-t0);
	return TRUE;
}

/*
 * Threshold filtering ("posterization")
 * Input:
 *		ima - picture (free() must be executed in the caller)
 *		f - filter:
 *			f-> w - number of levels of posterization, [2.255]
 *			f-> h - type of posterization (0 - uniform)
 *		sizex, sizey - image size
 *		min, max - minimum and maximum intensity of the image
 * Output:
 *		result - filtered image, the memory is allocated in this procedure
 *		scale - the scale of intensities, the memory is allocated here (if the scale!=NULL)
 *
 * TODO: save the result in the char, not float; learn display function
 */
int StepFilter(float *ima, float **result, Filter *f, int sizex, int sizey,
				float min, float max, float **scale){
	if(f->w < 2 || f->w > 255) return FALSE;
	int y;
	float Nsteps = (float)f->w; // amount of intervals
	float step;
	float wd = max - min;
	if(fabs(wd) < FLT_EPSILON) return FALSE;
	float (*stepfn)(float in);
	float Funiform(float in){
		return floor((in-min)/step);
	}
	float Flog(float in){
		return floor(logf(in-min+1.)/step);
	}
	float Fexp(float in){
		return floor(expf((in-min)/wd)/step);
	}
	float Fsqrt(float in){
		return floor(sqrtf(in-min)/step);
	}
	float Fpow(float in){
		return floor((in-min)*(in-min)/step);
	}
	#ifdef EBUG
	double t0 = dtime();
	#endif
	switch(f->h){
		case LOG:
			stepfn = Flog;
			step = logf(max-min+1.)/Nsteps;
		break;
		case EXP:
			stepfn = Fexp;
			step = expf(1.)/Nsteps;
		break;
		case SQRT:
			stepfn = Fsqrt;
			step = sqrtf(wd)/Nsteps;
		break;
		case POW:
			stepfn = Fpow;
			step = wd*wd/Nsteps;
		break;
		default:
			stepfn = Funiform;
			step = wd/Nsteps;
	}
	*result = calloc(sizex*sizey, sizeof(float));
	if(!result) return FALSE;
	#pragma omp parallel for
	for(y = 0; y < sizey; y++){
		int x; float *out, *in;
		in = ima + y*sizex; out = *result + y*sizex;
		for(x = 0; x < sizex; x++, in++, out++){
			*out = stepfn(*in);
		}
	}
	if(scale && !fillIsoScale(f, scale, min, wd)) g_err(_("No memory left"));
	DBG("SF: %d sublevels, step=%g, time=%f\n", f->w, step, dtime()-t0);
	return TRUE;
}




