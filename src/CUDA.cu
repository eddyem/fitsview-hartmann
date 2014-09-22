/*
 * 		CUDA.cu - subroutines for GPU
 *
 *      Copyright 2011 Edward V. Emelianoff <eddy@sao.ru>
 *
 *      This program is free software; you can redistribute it and/or modify
 *      it under the terms of the GNU General Public License as published by
 *      the Free Software Foundation; either version 2 of the License, or
 *      (at your option) any later version.
 *
 *      This program is distributed in the hope that it will be useful,
 *      but WITHOUT ANY WARRANTY; without even the implied warranty of
 *      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *      GNU General Public License for more details.
 *
 *      You should have received a copy of the GNU General Public License
 *      along with this program; if not, write to the Free Software
 *      Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 *      MA 02110-1301, USA.
 */
#define _CUDA_CU_
#include "include/CUtools.h"

#include <cuda.h>
#include <cufft.h>
//#include <cuda_runtime_api.h>
//#include <crt/device_runtime.h>
//#include <device_functions.h>

const int SHMEMSZ = 16383; // default constants, changed runtime
const int QBLKSZ  = 16;		// QBLKSZ = sqrt(LBLKSZ)
const int LBLKSZ  = 512;

// static arrays for sines & cosines values
static float *Sin_d = NULL, *Cos_d = NULL;
// array size
static int sincosize = 0;

cudaError_t CUerr;
inline int CUERROR(char *str){
	if(CUerr != cudaSuccess){
		fprintf(stderr, "%s, %s\n", str, cudaGetErrorString(CUerr));
		return 1;
	}else return 0;
}
// error macro (by default - nothing)
#define RETMACRO return
// memory macros
#define CUALLOC(var, size)		do{				\
	CUerr = cudaMalloc((void**)&var, size);		\
	if(CUERROR("CUDA: can't allocate memory")){	\
		RETMACRO;								\
}}while(0)
#define CUMOV2DEV(dest, src, size) do{			\
	CUerr = cudaMemcpy(dest, src, size,			\
				cudaMemcpyHostToDevice);		\
	if(CUERROR("CUDA: can't copy data to device")){\
		RETMACRO;}								\
}while(0)
#define CUMOV2HOST(dest, src, size) do{			\
	CUerr = cudaMemcpy(dest, src, size,			\
				cudaMemcpyDeviceToHost);		\
	if(CUERROR("CUDA: can't copy data to host")){\
		RETMACRO;}								\
}while(0)
#define CUFREE(var) do{cudaFree(var); var = NULL; }while(0)
#define  CUFFTCALL(fn)		do{					\
	cufftResult fres = fn;						\
	if(CUFFT_SUCCESS != fres){					\
		fprintf(stderr, "CUDA fft error %d\n", fres);\
		RETMACRO;}								\
}while(0)

#ifdef EBUG
	#define FNAME() fprintf(stderr, "\n%s (%s, line %d)\n", __func__, __FILE__, __LINE__)
	#define DBG(...) do{fprintf(stderr, "%s (%s, line %d): ", __func__, __FILE__, __LINE__); \
					fprintf(stderr, __VA_ARGS__);			\
					fprintf(stderr, "\n");} while(0)
#else
	#define FNAME()	 do{}while(0)
	#define DBG(...) do{}while(0)
#endif //EBUG

// getting the videocard parameters
extern "C" void getprops(){
	cudaDeviceProp dP;
	CUdevice dev; CUcontext ctx;
	cudaGetDeviceProperties(&dP, 0);
	cuDeviceGet(&dev,0);
	cuCtxCreate(&ctx, 0, dev);
	printf("\nDevice: %s, totalMem=%zd, memPerBlk=%zd,\n", dP.name, dP.totalGlobalMem, dP.sharedMemPerBlock);
	printf("warpSZ=%d, TPB=%d, TBDim=%dx%dx%d\n", dP.warpSize, dP.maxThreadsPerBlock,
			dP.maxThreadsDim[0],dP.maxThreadsDim[1],dP.maxThreadsDim[2]);
	printf("GridSz=%dx%dx%d, MemovrLap=%d, GPUs=%d\n", dP.maxGridSize[0],
			dP.maxGridSize[1],dP.maxGridSize[2],
			dP.deviceOverlap, dP.multiProcessorCount);
	printf("canMAPhostMEM=%d\n", dP.canMapHostMemory);
	printf("compute capability %d.%d.\n\n", dP.major, dP.minor);
	if(dP.major > 1){
	//	SHMEMSZ = 49151; QBLKSZ = 32; LBLKSZ = 1024;
	}
	size_t theFree, theTotal;
	CUresult aaa = cuMemGetInfo( &theFree, &theTotal );
	printf("CARD returns(err=%d):  free mem:%zd,  total mem:%zd\n", aaa, theFree, theTotal);
	cuCtxDetach(ctx);
}

// normalisation of array arr with size arrsize
__global__ void normalize_vec(float *arr, int arrsize){
	__shared__ float max[LBLKSZ];
	int idx = threadIdx.x;
	int blksize = (arrsize + blockDim.x - 1) / blockDim.x;
	int b_beg = idx * blksize;
	if(b_beg >= arrsize) return;
	int b_end = b_beg + blksize;
	if(b_end > arrsize) b_end = arrsize;
	int i; float *ptr = &arr[b_beg];
	float mm = *ptr++;
	for(i = b_beg +1 ; i < b_end; i++, ptr++)
		if(mm < *ptr) mm = *ptr;
	max[idx] = mm;
	__syncthreads();
	if(idx == 0){
		mm = max[0];
		for(i = 1; i < LBLKSZ; i++)
			if(mm < max[i]) mm = max[i];
		max[0] = mm;
	}
	__syncthreads();
	ptr = &arr[b_beg];
	mm = max[0];
	if(mm != 0.f)
		for(i = b_beg ; i < b_end; i++, ptr++) *ptr /= mm;
}
/*
__global__ void fill_zeros(float *arr, int arrsize, int W, int H){
	int x = blockIdx.x * blockDim.x + threadIdx.x;
	int y = blockIdx.y * blockDim.y + threadIdx.y;
	if(x > W || y > H) return;
	arr[x + y*W] = 0.f;
}*/

// kernel of function for sin/cos array initialisation
__global__ void fill_sincos(int angles,
							float *Sin_d, float *Cos_d,
							float anglestep, float conv){
	int k = blockIdx.x * blockDim.x + threadIdx.x;
	if(k >= angles) return;
	float theta = ((float)k / anglestep - 90.f) * conv;
	sincosf(theta, &Sin_d[k], &Cos_d[k]);
}

// initialisation of sin/cos arrays
extern "C" int init_sincos(int angles){
	#undef RETMACRO
	#define RETMACRO return 0
	// the value reciprocal for angle step
	float anglestep = (float)angles / 270.f;
	float conv = M_PI / 180.f;
	int blks = (angles + QBLKSZ - 1) / QBLKSZ;
	int threads = LBLKSZ;
	// first time we initialize arrays
	if(!Sin_d || !Cos_d || angles != sincosize){
		CUFREE(Cos_d);
		CUFREE(Sin_d);
		CUALLOC(Sin_d, angles*sizeof(float));
		CUALLOC(Cos_d, angles*sizeof(float));
		fill_sincos<<<blks, threads>>>(angles, Sin_d, Cos_d, anglestep, conv);
		sincosize = angles;
	}
	cudaThreadSynchronize();
	return 1;
	#undef RETMACRO
	#define RETMACRO return
}

/*
 * Lines Hough transform kernel
 * ima_d - device array with the image
 * imW, imH, min, max - width, height of the image and extreme values of its histogram
 * Sin_d, Cos_d - device array with sines and cosines of angles (-90..180degr increments 270/angles)
 * Rmax - the maximum range for the desired lines
 * angles - the number of angles in the range -90 .. 180
 * treshold - lower threshold of intensity (in relative units: I=tres*(max-min)+min) for inclusion of point into array
 * hough_d - output array with Hough transform
 */
__global__ void fill_lin_hough_array(float *ima_d,
						int imW, int imH,
						float min, float max,
						float *Sin_d, float *Cos_d,
						int Rmax, int angles, float treshold,
						float *hough_d){
	int xi = blockIdx.x * blockDim.x + threadIdx.x;
	int yi = blockIdx.y * blockDim.y + threadIdx.y;
	int i = xi + imW * yi;
	float x = (float)xi;
	float y = (float)yi;
	int k, R;
	if(xi >= imW || yi >= imH) return;
	float wd = max-min; if(wd == 0.f) wd = 1.f;
	float ima = (ima_d[i]-min)/wd;
	if(ima > treshold){
		for(k = 0; k < angles; k++){
			// R = x*cos(theta) + y*sin(theta)
			R = (int)(0.5f + x * Cos_d[k] + y * Sin_d[k]);
			// THIS IS VERY BAD, BUT atomicAdd doesn't work in old devices
			if(R > 0 && R < Rmax) hough_d[R + Rmax*k] += ima;
			//if(R > 0 && R < Rmax) atomicAdd(&hough_d[R + Rmax*k], ima);
		}
	}
}

/*
 * Build Hough transform to find lines
 * Input:
 *		ima - the image data
 *		min, max - range of the data in it
 *		imW, imH - image width and height
 *		Rmax - the maximum value of R
 *		Angles - the array size of angles (the angle of pitch is 180/angles degrees)
 * Output:
 *		hough - array initialized by an external function,
 *		in which the Hough transform will be
 *		!!! array must be initialized with zeros before calling this function
 * Output array is normalized to unity
 */
extern "C" int fill_hough_lines(float *ima, float min, float max, int imW, int imH,
								int Rmax, int angles, float *hough){
	#undef RETMACRO
	#define RETMACRO do{ ret = 0; goto free_all; }while(0)
	int sz, ret = 1;
	int lblksz = LBLKSZ;
	float *ima_d = NULL, *hough_d = NULL;
	float treshold = 0.1f;
	sz = imW * imH;
	getprops();
	dim3 blkdim(QBLKSZ, QBLKSZ);
	dim3 griddim((imW+QBLKSZ-1)/QBLKSZ, (imH+QBLKSZ-1)/QBLKSZ);
//	dim3 hgriddim((Rmax+QBLKSZ-1)/QBLKSZ, (angles+QBLKSZ-1)/QBLKSZ);
	if(!init_sincos(angles)) RETMACRO;
	CUALLOC(ima_d, sz*sizeof(float));
	CUMOV2DEV(ima_d, ima, sz*sizeof(float));
	sz = Rmax * angles;
	CUALLOC(hough_d, sz*sizeof(float));
	cudaMemset(hough_d, 0, sz*sizeof(float));
//	fill_zeros<<<hgriddim, blkdim>>>(hough_d, sz, Rmax, angles);
	cudaThreadSynchronize();
	CUMOV2DEV(hough_d, hough, sz*sizeof(float));
	fill_lin_hough_array<<<griddim, blkdim>>>(ima_d, imW,imH, min,max, Sin_d,Cos_d,
									Rmax, angles, treshold, hough_d);
	cudaThreadSynchronize();
	normalize_vec<<<1, lblksz>>>(hough_d, sz);
	cudaThreadSynchronize();
	CUMOV2HOST(hough, hough_d, sz*sizeof(float));
free_all:
	CUFREE(hough_d);
	CUFREE(ima_d);
	return ret;
	#undef RETMACRO
	#define RETMACRO return
}

/*
 * Kernels of the threshold filtering
 * in, out - in and out
 * stepfn - a pointer to a function of conversion
 * sizex, sizey - image size
 * min - the minimum intensity
 * wd - range of the data
 * step - a step for stepfn
 */
// uniform intensity distribution
__global__ void Funiform(float *in, int sizex, int sizey, float min, float step){
	int xi = blockIdx.x * blockDim.x + threadIdx.x;
	int yi = blockIdx.y * blockDim.y + threadIdx.y;
	int i = xi + yi*sizex;
	if(xi >= sizex || yi >= sizey) return;
	in[i] = floor((in[i]-min)/step);
}
// logarithm distribution
__global__ void Flog(float *in, int sizex, int sizey, float min, float step){
	int xi = blockIdx.x * blockDim.x + threadIdx.x;
	int yi = blockIdx.y * blockDim.y + threadIdx.y;
	int i = xi + yi*sizex;
	if(xi >= sizex || yi >= sizey) return;
	in[i] = floor(logf(in[i]-min+1.f)/step);
}
// exponential distribution
__global__ void Fexp(float *in, int sizex, int sizey, float min, float wd, float step){
	int xi = blockIdx.x * blockDim.x + threadIdx.x;
	int yi = blockIdx.y * blockDim.y + threadIdx.y;
	int i = xi + yi*sizex;
	if(xi >= sizex || yi >= sizey) return;
	in[i] = floor(expf((in[i]-min)/wd)/step);
}
// distribution of a square root
__global__ void Fsqrt(float *in, int sizex, int sizey, float min, float step){
	int xi = blockIdx.x * blockDim.x + threadIdx.x;
	int yi = blockIdx.y * blockDim.y + threadIdx.y;
	int i = xi + yi*sizex;
	if(xi >= sizex || yi >= sizey) return;
	in[i] = floor(sqrtf(in[i]-min)/step);
}
// distribution of a x^2
__global__ void Fpow(float *in, int sizex, int sizey, float min, float step){
	int xi = blockIdx.x * blockDim.x + threadIdx.x;
	int yi = blockIdx.y * blockDim.y + threadIdx.y;
	int i = xi + yi*sizex;
	if(xi >= sizex || yi >= sizey) return;
	in[i] = floor((in[i]-min)*(in[i]-min)/step);
}

// functions for calculation of output scale
float Suniform(float in, float min, float wd, float step){
	return step*in + min;
}
float Slog(float in, float min, float wd, float step){
	return expf(in*step) + min - 1.f;
}
float Sexp(float in, float min, float wd, float step){
	return wd*logf(in*step) + min;
}
float Ssqrt(float in, float min, float wd, float step){
	return in*in*step*step + min;
}
float Spow(float in, float min, float wd, float step){
	return sqrtf(in*step) + min;
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
extern "C"  int StepFilter(float *ima, float **result,
							Filter *f, int sizex, int sizey,
							float min, float max, float **scale){
	#undef RETMACRO
	#define RETMACRO do{ ret = 0; goto free_all; }while(0)
	int ret = 1;
	float wd = max - min;
	int y;
	float Nsteps = (float)f->w; // number of intervals
	float step;
	float *in=NULL; // image and result array
	int sz = sizex*sizey*sizeof(float);
	dim3 blkdim(QBLKSZ, QBLKSZ);
	dim3 griddim((sizex+QBLKSZ-1)/QBLKSZ, (sizey+QBLKSZ-1)/QBLKSZ);
	*result = (float*)malloc(sz);
	if(!result) RETMACRO;
	float (*scalefn)(float,float,float,float);
	if(f->w < 2 || f->w > 255) return 0;
	if(wd == 0.f) return 0;
	CUALLOC(in, sz);
	CUMOV2DEV(in, ima, sz);
	switch(f->h){ // filter type
		case LOG: // logarithm
			scalefn = Slog;
			step = logf(max-min+1.f)/Nsteps;
			Flog<<<griddim, blkdim>>>(in, sizex, sizey, min, step);
		break;
		case EXP: // exponential
			scalefn = Sexp;
			step = expf(1.f)/Nsteps;
			Fexp<<<griddim, blkdim>>>(in, sizex, sizey, min, wd, step);
		break;
		case SQRT: // square root
			scalefn = Ssqrt;
			step = sqrtf(wd)/Nsteps;
			Fsqrt<<<griddim, blkdim>>>(in, sizex, sizey, min, step);
		break;
		case POW: // power of two
			scalefn = Spow;
			step = wd*wd/Nsteps;
			Fpow<<<griddim, blkdim>>>(in, sizex, sizey, min, step);
		break;
		default: // uniform
			scalefn = Suniform;
			step = wd/Nsteps;
			Funiform<<<griddim, blkdim>>>(in, sizex, sizey, min, step);
	}
	cudaThreadSynchronize();
	CUMOV2HOST(*result, in, sz);
	if(scale){
		int M = f->w;
		*scale = (float*)calloc(M, sizeof(float));
		if(*scale) for(y = 0; y < M; y++){
			(*scale)[y] = scalefn(y+1,min,wd,step);
		}
	}
free_all:
	CUFREE(in);
	return ret;
	#undef RETMACRO
	#define RETMACRO return
}


/*
 * A set of functions for constructing differential filters
 */
int p2oi(int i){
	unsigned int v = (unsigned int)i - 1;
	v |= v >> 1;
	v |= v >> 2;
	v |= v >> 4;
	v |= v >> 8;
	v |= v >> 16;
	v++;
	return (int) v;
}
int nextpow2(int i, int j){
	int p1 = p2oi(i), p2 = p2oi(j);
	return (p1 > p2)? p1 : p2;
}
// multiplication of two complex matrices with size x size
// result in the entry of the first matrix
__global__ void ComplexMul(cufftComplex *inout, cufftComplex *in, int size){
	int xi = blockIdx.x * blockDim.x + threadIdx.x;
	int yi = blockIdx.y * blockDim.y + threadIdx.y;
	if(xi >= size || yi >= size) return;
	int i = xi + yi*size;
	cufftComplex a = inout[i], b = in[i];
	inout[i].x = a.x * b.x - a.y * b.y;
	inout[i].y = a.x * b.y + a.y * b.x;
}

// restore coordinates of the Fourier transform
__global__ void fftshift(int size, float *m){
	int h = size/2;
	int xi = blockIdx.x * blockDim.x + threadIdx.x;
	int yi = blockIdx.y * blockDim.y + threadIdx.y;
	if(xi >= h || yi >= h) return;
	// k - point in left upper quadrant, k1 - in right upper
	int k = yi * size+xi, k1 = k + h;
	// p - point in right lower quadrant, p1 - in left lower
	int p = k + (size+1)*h, p1 = k1 + (size-1)*h;
	float tmp;
	tmp = m[k]; m[k] = m[p]; m[p] = tmp;
	tmp = m[k1]; m[k1] = m[p1]; m[p1] = tmp;
}
// data copying float->cufftReal (need because different sizes of picture and Fourier image)
__global__ void f2r(cufftReal *out, float *in, int sizex, int sizey, int size2){
	int xi = blockIdx.x * blockDim.x + threadIdx.x;
	int yi = blockIdx.y * blockDim.y + threadIdx.y;
	if(xi >= sizex || yi >= sizey) return;
	out[xi+yi*size2] = (cufftReal) in[xi+yi*sizex];
}
// data copying cufftReal->float
__global__ void r2f(float *out, cufftReal *in, int sizex, int sizey, int size2){
	int xi = blockIdx.x * blockDim.x + threadIdx.x;
	int yi = blockIdx.y * blockDim.y + threadIdx.y;
	if(xi >= sizex || yi >= sizey) return;
	out[xi+yi*sizex] = (float) in[xi+yi*size2];
}

/*
 * The kernel of the Laplacian of Gaussian
 * Output:
 *		mask - the filled array
 * Input:
 *		size - array size (size x size)
 *		x0, x1 - array bounds on x: [x0, x1)  (outside this array filled by zeros)
 *		y0, y1 - -//- on y
 *		half - half the size of the array
 *		ss - normalizing factor
 *		sx2, sy2 - the variance of the filter in x and y
 */
__global__ void LGf_kernel(cufftReal *mask, int size, int x0, int x1,
							int y0, int y1, float half, float ss,
							float sx2, float sy2){
	int xi = blockIdx.x * blockDim.x + threadIdx.x + x0;
	int yi = blockIdx.y * blockDim.y + threadIdx.y + y0;
	if(xi >= x1 || yi >= y1) return;
	int i = xi + yi*size;
	float x2 = (float)xi + half;
	float y2 = (float)yi + half;
	x2 = x2*x2 / sx2; y2 = y2*y2 / sy2;
	mask[i] = (cufftReal)(ss * ((x2-1.f)/sx2+(y2-1.f)/sy2)*expf(-(x2+y2)/2.f));
}
// The kernel of Gaussian filter
__global__ void Gf_kernel(cufftReal *mask, int size, int x0, int x1,
							int y0, int y1, float half, float ss,
							float sx2, float sy2){
	int xi = blockIdx.x * blockDim.x + threadIdx.x + x0;
	int yi = blockIdx.y * blockDim.y + threadIdx.y + y0;
	if(xi >= x1 || yi >= y1) return;
	int i = xi + yi*size;
	float x2 = (float)xi + half;
	float y2 = (float)yi + half;
	x2 = x2*x2 / sx2; y2 = y2*y2 / sy2;
	mask[i] = (cufftReal)(ss * expf(-(x2+y2)/2.f));
}
/*
 * Building mask of Gaussian or Laplasian of Gaussian
 * Output:
 *		mask - filter array
 * Input:
 * 		size - mask size (size x size)
 * 		f - filter parameters
 */
void build_GLG_filter(cufftReal *mask, int size, Filter *f){
	int y0=0,y1=size, x0=0, x1=size;
	float sx2 = f->sx * f->sx, sy2 = f->sy * f->sy;
	float half;
	dim3 blkdim(QBLKSZ, QBLKSZ);
	dim3 griddim((size+QBLKSZ-1)/QBLKSZ, (size+QBLKSZ-1)/QBLKSZ);
	if(f->w < size && f->w > 0){
		x0 = (size - f->w + 1) / 2;
		x1 = x0 + f->w;
	}
	if(f->h < size && f->h > 0){
		y0 = (size - f->h + 1) / 2;
		y1 = y0 + f->h;
	}
	half = -(float)size / 2.f;
	float ss = 3.f / half / half / sqrt(-half);
	switch(f->FilterType){
		case LAPGAUSS:
			LGf_kernel<<<griddim, blkdim>>>(mask, size, x0,x1, y0,y1, half, ss, sx2, sy2);
		break;
		case GAUSS:
			Gf_kernel<<<griddim, blkdim>>>(mask, size, x0,x1, y0,y1, half, ss, sx2, sy2);
		break;
		default:
			fprintf(stderr, "Error: bad filter\n");
	}
	cudaThreadSynchronize();
	DBG("size=%d, x0=%d,x1=%d, y0=%d,y1=%d, half=%g, ss=%g, sx2=%g, sy2=%g",
		size, x0,x1, y0,y1, half, ss, sx2, sy2);
}

/*
 * Building of elementary filter mask
 * Output:
 *		mask - filter array
 * Input:
 * 		size - mask size (size x size)
 * 		f - filter parameters
 */
void build_S_filter(cufftReal *mask, int size, Filter *f){
	int y, a0, a1;
	float hh, Y, pt = 0.f;
	a0 = (size - 2) / 2;
	a1 = a0 + 3;
	hh = -(float)(size / 2);
	float ss = 1.f / (M_PI*2.f) / hh / hh / sqrtf(-hh);
	Y = -1.f;
	for(y = a0; y < a1; y++, Y+=1.f){
		float X = -1.f;
		int str, x;
		str = y * size;
		for(x = a0; x < a1; x++, X+=1.f){
			switch(f->FilterType){
				case SOBELH:
					pt = -X*(2.f-fabs(Y));
				break;
				case SOBELV:
					pt = -Y*(2.f-fabs(X));
				break;
				case PREWITTH:
					pt = X;
				break;
				case PREWITTV:
					pt = Y;
				break;
			}
			cufftReal tmppar = (cufftReal)ss*pt;
			cudaMemcpy(&mask[str + x], &tmppar, sizeof(cufftReal), cudaMemcpyHostToDevice);
		}
	}
}
/*
 * Convolution filtering (convolution by FFT)
 * Input:
 * 		ima - picture, that need to be filtering
 * 		f - filter parameters
 * Output:
 * 		result - memory area, allocated by this function,
 * 				where the filtered picture to be store
 * return: TRUE if the filtering succeed
 */
extern "C"  int DiffFilter(float *ima, float **result,
					Filter *f, int sizex, int sizey){
	#undef RETMACRO
	#define RETMACRO do{ ret = 0; goto free_all; }while(0)
	int ssize, ret = 0, size2;
	float *tmp;
	size2 = nextpow2(sizex, sizey);
	ssize = size2 * size2; // Fourier image size
	dim3 blkdim(QBLKSZ, QBLKSZ);
	dim3 griddim((size2+QBLKSZ-1)/QBLKSZ, (size2+QBLKSZ-1)/QBLKSZ);
	dim3 halfgriddim((size2/2+QBLKSZ-1)/QBLKSZ, (size2/2+QBLKSZ-1)/QBLKSZ);
	dim3 imgriddim((sizex+QBLKSZ-1)/QBLKSZ, (sizey+QBLKSZ-1)/QBLKSZ);
	if(!result || !*result || !ima || !f){
		fprintf(stderr, "DiffFilter: bad parameters\n");
		return 0;
	}
	cufftHandle plan;
	cufftComplex *Fmask=NULL, *Fimg=NULL;
	cufftReal *mask=NULL, *img=NULL, *resm=NULL;
	#ifdef EBUG
	getprops();
	#endif
	// Allocate memory for new objects
	DBG("allocate");
	CUALLOC(img, ssize*sizeof(cufftReal));
	// fill it zeros
	cudaMemset(img, 0, ssize*sizeof(cufftReal));
	// copy ima -> img
	DBG("copy image to dev");
	CUALLOC(tmp, sizex*sizey*sizeof(float));
	CUMOV2DEV(tmp, ima, sizex*sizey*sizeof(float));
	f2r<<<imgriddim, blkdim>>>(img, tmp, sizex, sizey, size2);
	cudaThreadSynchronize();
	CUFREE(tmp);
	CUALLOC(Fimg, ssize*sizeof(cufftComplex));
	// make FFT
	DBG("doing image FFT");
	CUFFTCALL(cufftPlan2d(&plan, size2, size2, CUFFT_R2C));
	CUFFTCALL(cufftExecR2C(plan, img, Fimg));
	CUFREE(img);
	DBG("allocate");
	CUALLOC(mask, ssize*sizeof(cufftReal));
	cudaMemset(mask, 0, ssize*sizeof(cufftReal));
	CUALLOC(Fmask, ssize*sizeof(cufftComplex));
	switch(f->FilterType){
		case LAPGAUSS:
		case GAUSS:
			build_GLG_filter(mask, size2, f);
			break;
		case SOBELH:
		case SOBELV:
		case PREWITTH:
		case PREWITTV:
			build_S_filter(mask, size2, f);
			break;
		default:
			fprintf(stderr, "Error: bad filter\n");
			RETMACRO;
	}
	// swap filter quadrants
	fftshift<<<halfgriddim, blkdim>>>(size2, mask);
	cudaThreadSynchronize();
	// make FFT
	DBG("doing filter FFT");
	CUFFTCALL(cufftExecR2C(plan, mask, Fmask));
	CUFFTCALL(cufftDestroy(plan));
	CUFREE(mask);
	// make convolution in Fourier space
	DBG("multiplication");
	ComplexMul<<<griddim, blkdim>>>(Fimg, Fmask, size2);
	cudaThreadSynchronize();
	CUFREE(Fmask);
	// Inverse FFT
	DBG("doing inverse FFT");
	CUALLOC(resm, ssize*sizeof(cufftReal));
	CUFFTCALL(cufftPlan2d(&plan, size2, size2, CUFFT_C2R));
	CUFFTCALL(cufftExecC2R(plan, Fimg, resm));
	CUFFTCALL(cufftDestroy(plan));
	CUFREE(Fimg);
	DBG("allocate");
	CUALLOC(tmp, ssize*sizeof(float));
	*result = (float*)calloc(sizex*sizey, sizeof(float));
	if(!*result) RETMACRO;
	// copy iFFT -> res
	DBG("copy to host");
	r2f<<<imgriddim, blkdim>>>(tmp, resm, sizex, sizey, size2);
	cudaThreadSynchronize();
	CUMOV2HOST(*result, tmp, sizex*sizey*sizeof(float));
	ret = 1;
free_all:
	CUFREE(Fmask); CUFREE(Fimg); CUFREE(img);
	CUFREE(mask); CUFREE(resm); CUFREE(tmp);
	#ifdef EBUG
	getprops();
	#endif
	return ret;
	#undef RETMACRO
	#define RETMACRO return
}
extern "C"  int MedFilter(float *ima, float **result, Filter *f, int sizex, int sizey){return 0;}
extern "C"  int GradFilterSimple(float *ima, float **result, Filter *f, int sizex, int sizey){return 0;}
