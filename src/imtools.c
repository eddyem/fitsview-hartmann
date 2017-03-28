//      imtools.c - different image transforms
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
#include "fitsheaders.h"
#include "imtools.h"
#include "contours.h"
#include "opengl.h"
#include "gauss.h"
#include "CUtools.h"

// Histogram size, 256 by default
const int HIST_SIZE	= 256;
const int MIN_SPOT_NO = 200;	// minimum amount of spots on a hartmanogram
const int MAX_SPOT_NO = 258;	// max amount -//-

void destroy_image(IMAGE *image){
	if(!image) return;
	_FREE(image->data);
	g_free(image->filename);
	if(image->store) // clear parameters from FITS header
		gtk_list_store_clear(image->store);
	_FREE(image->cLevels); // isolines' intensity
	free_contours(image);  // isolines itself
	_FREE(image);
}

#ifdef LEPTONICA_FOUND
/*
 * Convert an image ima (after tresholding)
 * into a pixmap for pixConnComp from leptonica
 */
PIX *conv_image2pix(IMAGE *ima){
	PIX *pixs;
	l_uint32 byte, *pptr;
	GLfloat *iptr, *idata;
	int w, h, wpl, i, j, k, N;
	w = ima->width;
	h = ima->height;
	pixs = pixCreate(w, h, 1);
	wpl = pixGetWpl(pixs);
	pptr = pixGetData(pixs);
	idata = ima->data;
	for(i = h - 1; i > -1; i--){ // lines cycle
		N = 0;
		iptr = &idata[i * w];
		for(j = 0; j < wpl; j++){ // cycle by 32 columns
			byte = 0;
			for(k = 0; k < 32; k++){ // cycle by a 32-bit word bits
				byte <<= 1;
				// if string still doesn't end & intensity higher 0.1%
				if(N++ < w && *iptr++ > 0.001)
					byte |= 1; // set lesser bit
			}
			*pptr++ = byte;
		}
	}
	return pixs;
}

gboolean check_box(Window *window, BOX *box, int *X0, int *Y0, int *X, int *Y){
	if(!window->image->data) return FALSE;
	int WD = window->image->width, HT = window->image->height - 1;
	if(!(box && X0 && Y0 && X && Y)) return FALSE;
	*X0 = box->x;
	*Y0 = box->y;
	*X  = *X0 + box->w;
	*Y  = *Y0 + box->h; // coordinates of a box corners
	if(box->w < 5 || box->h < 5) return FALSE;
	if(*X0 < 0 || *X0 > WD || *X < 0 || *X > WD ||
	   *Y0 < 0 || *Y0 > HT || *Y < 0 || *Y > HT)
			return FALSE;
	return TRUE;
}

 //Calculation of min & max intensity by picture to show it
void compute_minmax(IMAGE *ima){
	float min, max, *dst = ima->data;
	int i, j, w=ima->width, h=ima->height;
	max = min = *dst;
	for(i=0; i<h; i++)
		for(j=0; j<w; j++, dst++){
			float tmp = *dst;
			if(tmp > max) max = tmp;
			else if(tmp < min) min = tmp;
		}
	ima->stat.max = max;
	ima->stat.min = min;
	//ima->width = w; ima->height = h;
	DBG("stat: max=%.2f, min=%.2f", max, min);
}

/*
 * Calculation of the center of box for image ima
 * by inscribing a Gaussian to a vector obtained
 * by addition of box columns (for X) and rows (for Y)
 */
gboolean get_gaussian_center(Window *window, BOX *box, Coordinates *crds){
	int WD = window->image->width, HT = window->image->height - 1;
	int x, y;
	GLfloat *data;
	Points row, col;
	Point *rptr, *cptr;
	double x0, sigma;
	int X0, Y0, X, Y;
	#ifndef GSL_FOUND
		return FALSE;
	#endif
	if(!check_box(window, box, &X0, &Y0, &X, &Y)) return FALSE;
	if(!crds) return FALSE;
	row.data = calloc(box->w, sizeof(Point));
	row.n = box->w;
	col.data = calloc(box->h, sizeof(Point));
	col.n = box->h;
	cptr = col.data;
	for(y = Y0; y < Y; y++, cptr++){
		data = &window->image->data[(HT - y) * WD + X0];
		cptr->x  = (double)y;
		rptr = row.data;
		for(x = X0; x < X; x++, data++, rptr++){
			double d = (double) *data;
			rptr->y += d;
			cptr->y += d;
		}
	}
	rptr = row.data;
	for(x = X0; x < X; x++, rptr++) rptr->x  = (double)x;
	// we inscribe Gaussian and change the the coordinates and half-widths of the spot
	gauss_fit(&row, NULL, NULL, &sigma, &x0);
	crds->x = x0; crds->rx = sigma;
	gauss_fit(&col, NULL, NULL, &sigma, &x0);
	crds->y = x0; crds->ry = sigma;
	_FREE(row.data); _FREE(col.data);
	return TRUE;
}

// count of bits amount in word i
// http://graphics.stanford.edu/~seander/bithacks.html#CountBitsSetParallel
// http://en.wikipedia.org/wiki/Hamming_weight
int count_bits(l_uint32 i){
	i = i - ((i >> 1) & 0x55555555);
	i = (i & 0x33333333) + ((i >> 2) & 0x33333333);
	return ((i + ((i >> 4) & 0xF0F0F0F)) * 0x1010101) >> 24;
}

gboolean get_circle_center(Window *window, BOX *box, PIX *pix, Coordinates *crds){
	int WD = window->image->width, HT = window->image->height - 1;
	int N = 0, X0, Y0, X, Y, n;
	l_uint32 *pixdata, *pixptr;
	l_int32 w, h, d, wpl, sz, i, x, y;
	double Xcenter, Ycenter, xx, yy, Intens, X2, Y2;
	GLfloat *data;
	if(!check_box(window, box, &X0, &Y0, &X, &Y)) return FALSE;
	if(!(pix && crds)) return FALSE;
	pixGetDimensions(pix, &w, &h, &d);
	if(d != 1) return FALSE;
	wpl = pixGetWpl(pix);
	pixdata = pixGetData(pix);
	pixptr = pixdata;
	sz = wpl * h;
	for(i = 0; i < sz; i++) N += count_bits(*pixptr++);
	if(N < 200000) return FALSE;
	n = 0;
	Xcenter = (double)(X0 + X) / 2.;
	Ycenter = (double)(Y0 + Y) / 2.;
	DBG("x0=%g, y0=%g, r0=%g, N=%d", Xcenter, Ycenter, (double)(w + h) / 4., N);
	Xcenter = Ycenter = X2 = Y2 = Intens = 0.;
	AxisX = 0.; AxisY = (double)HT;
	for(y = 0; y < h; y++){
		yy = (double)(HT - y - Y0);
		data = &window->image->data[(HT - y - Y0) * WD + X0];
		l_uint32 *line = pixdata + y * wpl;
		for(x = 0; x < w; x++, data++){
			if(GET_DATA_BIT(line, x)){
				double d = *data;
				xx= (X0 + x);
				n++;
				Intens += d;
				Xcenter += xx * d;
				Ycenter += yy * d;
				X2 += xx*xx*d;
				Y2 += yy*yy*d;
			}
		}
	}
	xx = Xcenter / Intens;
	yy = Ycenter / Intens;
	crds->x = xx;
	crds->y = yy;
	crds->rx= sqrt(X2 / Intens - xx*xx);
	crds->ry= sqrt(Y2 / Intens - yy*yy);
	DBG("xx: %g, yy: %g, rx:%g, ry: %g", xx, yy, crds->rx, crds->ry);
	return TRUE;
}

// Detection of spots' center, yet lousy through gaussians
void get_spots_coords(Window *window, IMAGE *ima){
	BOXA *boxa;
	BOX *box;
	PIX *pixs;
	Spot spot;
	Spots *spots;
	Coordinates scrds;
	int n, i, k;
	float xc, yc;
	if(!( /*(window->context->transformMask & TRANS_TRES) &&*/ window->image_transformed->data)){
		g_err(_("Define treshold limits first"));
		window->context->graphMode = GR_HISTOGRAM;
		show_histogram(window);
		return;
	}
	pixs = conv_image2pix(ima);
	boxa = pixConnComp(pixs, NULL, 8);
	n = boxaGetCount(boxa);
	red("Number of boxes: %d", n);
	spots_free(&window->spots);
	if(!(window->spots = spots_alloc(n))) return;
	spots = window->spots;
	#ifndef GSL_FOUND
		g_err(_("No GSL library found, don't calculate centroids"));
	#endif
	for(i = 0, k = 0; i < n; i++){
		if((box = boxaGetBox(boxa, i, L_CLONE)) == NULL){
			g_err(_("Box not found"));
			continue;
		}
		if(	box->w < Global->minSpotW ||
			box->h < Global->minSpotH ||
			box->w > Global->maxSpotW ||
			box->h > Global->maxSpotH){
			boxDestroy(&box);
			continue;
		}
		if(!get_gaussian_center(window, box, &scrds)){
			scrds.rx = ((double)box->w-1.)/2.;
			scrds.ry = ((double)box->h-1.)/2.;
			xc = box->x + scrds.rx;
			yc = box->y + scrds.ry;
			scrds.x = xc; scrds.y = yc;
		}
		spot.c = scrds; spot.id = k++;
		spot.box.x = box->x;
		spot.box.y = box->y;
		spot.box.w = box->w;
		spot.box.h = box->h;
		spots_add(spots, &spot, BY_COPY);
		boxDestroy(&box);
	}
	red("Number of spots: %d", k);
	boxaDestroy(&boxa);
	gdk_window_invalidate_rect(window->drawingArea->window,
		&window->drawingArea->allocation, FALSE);
	DBG("return");
}

void hartmann_spots(Window *window){ // 'r'
	Spots *spots = window->spots;
	if(!spots){
		g_err(_("Find spots first"));
		return;
	}
	int n = spots->n;
	if(n > MAX_SPOT_NO){
		g_err(_("Too many points"));
		return;
	}
	if(n < MIN_SPOT_NO){
		g_err(_("Not enough points"));
		return;
	}
	sort_spots(window);
	window->context->transformMask |= TRANS_HARTMANN; // points are sorted
	gdk_window_invalidate_rect(window->drawingArea->window,
		&window->drawingArea->allocation, FALSE);
}

void get_circles_params(Window *window, IMAGE *ima){
	BOXA *boxa;
	PIXA *pixa;
	BOX *box;
	PIX *pixs, *pix;
	Spot spot;
	Coordinates scrds;
	Spots *spots;
	int n, i;
	if(!( (window->context->transformMask & TRANS_TRES) && window->image_transformed->data)){
		g_err(_("Define treshold limits first"));
		window->context->graphMode = GR_HISTOGRAM;
		show_histogram(window);
		return;
	}
	pixs = conv_image2pix(ima);

	boxa = pixConnComp(pixs, &pixa, 8);
	n = boxaGetCount(boxa);

	red("Number of boxes: %d", n);
	spots_free(&window->spots);
	if(!(window->spots = spots_alloc(n))) return;
	spots = window->spots;
	for(i = 0; i < n; i++){
		if((box = boxaGetBox(boxa, i, L_CLONE)) == NULL){
			g_err(_("Box not found"));
			continue;
		}
		if((pix = pixaGetPix(pixa, i, L_CLONE)) == NULL){
			g_err(_("Pix not found"));
			continue;
		}
		if(get_circle_center(window, box, pix, &scrds)){
			spot.c = scrds; spot.id = i;
			spot.box.x = box->x;
			spot.box.y = box->y;
			spot.box.w = box->w;
			spot.box.h = box->h;
			spots_add(spots, &spot, BY_COPY);
		}
		boxDestroy(&box);
	}

	boxaDestroy(&boxa);
	pixaDestroy(&pixa);

	gdk_window_invalidate_rect(window->drawingArea->window,
		&window->drawingArea->allocation, FALSE);
	DBG("return");
}
#else // LEPTONICA_FOUND not defined
void get_spots_coords(Window *window __attribute__((unused)),
			IMAGE *ima __attribute__((unused))){
	LEPTERR;
}
gboolean get_gaussian_center(Window *window __attribute__((unused)),
			BOX *box __attribute__((unused)), Coordinates *crds __attribute__((unused))){
	LEPTERR;
	return FALSE;
}
PIX *conv_image2pix(IMAGE *ima __attribute__((unused))){
	LEPTERR;
	return NULL;
}
void get_circles_params(Window *window __attribute__((unused)),
			IMAGE *ima __attribute__((unused))){
	LEPTERR;
}
#endif  // LEPTONICA_FOUND


/*
 * Build the histogram of picture ima
 * The default is 256 shades of gray
 */
gboolean make_histogram(IMAGE *ima){
	FNAME();
	ssize_t i, sz;
	int j, pt, *hdata, max, min; // hmax
	GLfloat *ptr, hsize, datamin = ima->stat.min, datawd = ima->stat.max - datamin;
	#ifdef EBUG
		double t0 = dtime();
	#endif
	if(!ima->data){
		g_err(_("Image is empty"));
		return FALSE;
	}
	_FREE(ima->stat.histogram.data);
	// to change the size of the histogram replace HIST_SIZE by a variable value
	ima->stat.histogram.data = calloc(HIST_SIZE+1,  sizeof(int));
	if(!ima->stat.histogram.data){
		g_err(_("Can't allocate memory for histogram"));
		return FALSE;
	}
	hdata = ima->stat.histogram.data;
	ima->stat.histogram.size = HIST_SIZE;
	hsize = (GLfloat) HIST_SIZE / datawd;
	ima->stat.histogram.scale = 1. / hsize;
	//hmax = HIST_SIZE - 1;
	ptr = ima->data;
	sz = (ssize_t)(ima->width) * (ssize_t)(ima->height);
	DBG("datamin=%g, hsize=%g, sz=%zd, scale=%g", datamin, hsize, sz, ima->stat.histogram.scale);
	#pragma omp parallel for private(i, pt) shared(hdata)
	for(i = 0; i < sz; i++){
		pt = (int)((ptr[i] - datamin) * hsize);
		#pragma omp atomic
		hdata[pt]++;
	}
/*	for(i = 0; i < sz; i++, ptr++){
		pt = (int)(*ptr * hsize);
		if(pt < 0) pt = 0;
		else if(pt > hmax) pt = hmax;
		hdata[pt]++;
	}*/
	max = min = hdata[0];
	for(j = 1; j < HIST_SIZE; j++){
		pt = hdata[j];
		if(max < pt) max = pt;
		else if(min > pt) min = pt;
	}
	ima->stat.histogram.max = max;
	ima->stat.histogram.min = min;
	ima->stat.histogram.bot_tres = -1.;
	ima->stat.histogram.top_tres = -1.;
	DBG("max: %d, min: %d; time: %g", max, min, dtime() - t0);
	return TRUE;
}

/*
 * Fill the structure image_transformed by adjusted by the threshold values
 * Bottom - the lower bound, top - upper bound (in units of image intensity)
 * Brightness less than bottom, replaced by zero
 * Brightness larger than top, replaced by a unit
 * The rest - are calculated linearly (equalization histogram)
 */
void treshold_image(Window *window, GLfloat bottom, GLfloat top){
	FNAME();
	ssize_t i, sz;
	IMAGE *ima = window->image;
	IMAGE *i_new = window->image_transformed;
	GLfloat *p_old, *p_new, pt, wd;
	GLfloat scale = ima->stat.histogram.scale * (GLfloat)ima->stat.histogram.size;
	#ifdef EBUG
		double t0 = dtime();
	#endif
	if(!(p_old = ima->data)){
		g_err(_("Image is empty"));
		return;
	}
	_FREE(i_new->data);
	if((top - bottom) < 0.01) top = bottom + 0.01;
	sz = (ssize_t)(ima->width) * (ssize_t)(ima->height);
	i_new->data = malloc(sz * sizeof(GLfloat));
	i_new->width = ima->width; i_new->height = ima->height;
	i_new->bZero = 0.; i_new->bScale = 1.;
	i_new->maxVal = ima->maxVal;
	i_new->stat.max = 1.;//top * ima->stat.max;
	i_new->stat.min = 0.;
	p_new = i_new->data;
	DBG("top=%g, bottom=%g, scale=%g", top, bottom, ima->stat.histogram.scale);
	top = top * scale+ima->stat.min;
	bottom = bottom * scale+ima->stat.min;
	DBG("top=%g, bottom=%g", top, bottom);
	wd = top - bottom;
	for(i = 0; i < sz; i++, p_new++, p_old++){
		pt = *p_old;
		if(pt < bottom) pt = 0.;
		else if(pt > top) pt = 1.;
		else pt = (pt - bottom) / wd;
		*p_new = pt;
	}
	DBG("time: %g", dtime() - t0);
	window->context->transformMask |= TRANS_TRES;
	window->context->current_image = window->image_transformed;
	gen_texture(window->image_transformed, window, TRUE);
}

const int SINCOSSIZE = 540; // size of the array of sines & cosines
							// [-90 , +180) degrees

void hough_lines(Window *window){ // key 'i'
	FNAME();
	int w, h, Rmax, sz;
	double hd, wd;
	IMAGE *ima = window->image;
	Window *outwin = init_window(window, OPENGL_WINDOW);
	IMAGE *i_new = outwin->image;
	if(!ima->data){
		g_err(_("Image is empty"));
		return;
	}
	_FREE(i_new->data);
	w = ima->width; h = ima->height;
	wd = (double)w; hd = (double)h;
	Rmax = (int)(0.5 + sqrt(wd*wd + hd*hd));
	sz = SINCOSSIZE * Rmax;
	DBG("sz: %d", sz);
	i_new->data = calloc(sz, sizeof(GLfloat));
	if(!i_new->data){
		g_err(_("Can't allocate memory for Hough transform"));
		return;
	}
	i_new->width = Rmax; i_new->height = SINCOSSIZE;
	i_new->bZero = 0.; i_new->bScale = 1.;
	i_new->maxVal = 1.;
	i_new->dtype = FLOAT_IMG;
	gchar *str = malloc(256);
	snprintf(str, 255, "Hough transform of file %s", ima->filename);
	init_keylist(i_new);
	add_key(i_new, "TSTRING", "OBJECT", "Hough", str, TSTRING);
	_FREE(str);
#ifdef EBUG
	double t0 = dtime();
#endif
	if(!fill_hough_lines(ima->data, ima->stat.min, ima->stat.max, w, h, Rmax, SINCOSSIZE, i_new->data)){
		g_err(_("Error in Hough transform module"));
		destroy_window(outwin);
		return;
	}
	DBG("Hough transform for lines, time: %gs", dtime() - t0);
	//Context->transformMask |= TRANS_TRES;
	//DBG("HOUGH: maxR=%d, size=%d, max: %g at R=%d, theta=%g",
		//i_new->width, i_new->height, max, R%Rmax, ((double)(R/Rmax))/2.-90.);
	outwin->context->current_image = outwin->image;
	compute_minmax(outwin->image);
	gen_texture(outwin->image, outwin, FALSE);
	force_redraw(outwin->drawingArea);
	outwin->image->imagetype = HOUGH;
	CH_WIN_TITLE(outwin, "%s (Hough) %s", PRGNAME, window->image->filename);
}

/*
 * respect to the coordinates (x, y) of a click in the window of a Hough transform
 * we obtain the parameters of the line (R, phi) and display in the parent window
 * the corresponding line
 */
void get_houg_line( Window *window,
					double x, double y, // mouse click coordinates
					double *phi_, double *R_ // Output: parameters of line
					){
	double R, phi, ang_sep, sinphi, cosphi;
	double x1=0.,y1=0., x2=0.,y2=0., yl,yr, xu,xd, X1, X2, Y1, Y2;
	if(!window || !window->parent) return;
	double W = window->parent->image->width-1., H = window->parent->image->height-1.;
	int rdy = 0;
	IMAGE *image = window->image;
	conv_mouse_to_image_coords(x, y, &R, &phi, window);
	ang_sep = 270. / (double)image->height / 180. * M_PI;
	phi = phi*ang_sep - M_PI/2.;
	sinphi = sin(phi); cosphi = cos(phi);
	phi *= 180./M_PI; // convert rad => deg
	DBG("line: phi = %g, R = %g", phi, R);
	*phi_ = phi; *R_ = R;
	GdkDrawable *d = window->parent->drawingArea->window;
	static GdkGC *gc = NULL;
	if(!gc){
		gc = gdk_gc_new(d);
		gdk_gc_set_foreground(gc, &(window->parent->drawingArea->style->white));
		gdk_gc_set_function(gc, GDK_XOR);
	}
	yl = R/sinphi; yr = (R-W*cosphi)/sinphi; // left & right bounds of a picture
	xd = R/cosphi; xu = (R-H*sinphi)/cosphi; // down & up bounds
	DBG("yl=%g, yr=%g, xd=%g, xu=%g",yl,yr,xd,xu);
	if(yl < H && yl > -0.5){
		DBG("left selected");
		y1 = yl; x1 = 0.; rdy = 1;
	}
	if(yr < H && yr > -0.5){
		DBG("right selected");
		if(rdy == 0){
			y1 = yr; x1 = W; rdy = 1;
		}else{
			y2 = yr; x2 = W; rdy = 2; goto lines_selected;
		}
	}
	if(xd < W && xd > -0.5){
		DBG("down selected");
		if(rdy == 0){
			x1 = xd; y1 = 0.; rdy = 1;
		}else{
			x2 = xd; y2 = 0.; rdy = 2; goto lines_selected;
		}
	}
	if(xu < W && xu > -0.5){
		DBG("up selected");
		if(rdy == 0){
			x1 = xu; y1 = H; rdy = 1;
		}else{
			x2 = xu; y2 = H; rdy = 2;
		}
	}
lines_selected:
	if(rdy != 2) DBG("WTF? not two points");
	else{
		DBG("line in image CS: (%g, %g) - (%g, %g)", x1,y1,x2,y2);
		conv_image_to_mouse_coords(x1, y1, &X1, &Y1, window->parent);
		conv_image_to_mouse_coords(x2, y2, &X2, &Y2, window->parent);
		gdk_draw_line(d, gc, X1,Y1, X2,Y2);
		DBG("line in window CS: (%g, %g) - (%g, %g)", X1,Y1,X2,Y2);
	}
}

// Filtering of an picture window->image by filter f
void filter_image( Window *window,
					Filter *f
					){
	float *src, *dst;
	int w,h;
	gboolean res = FALSE;
	// If picture already transformed
	if(window->context->current_image == window->image_transformed){
		src = window->image_transformed->data;
		w = window->image_transformed->width;
		h = window->image_transformed->height;
	}else{
		src = window->image->data;
		w = window->image->width;
		h = window->image->height;
		if(window->image_transformed) _FREE(window->image_transformed->data);
//		window->image_transformed = COPY(window->image, IMAGE);
//		window->image_transformed->data = NULL;
	}
	switch(f->FilterType){
		case MEDIAN:
			res = MedFilter(src, &dst, f, w, h);
		break;
		case SIMPLEGRAD:
			res = GradFilterSimple(src, &dst, f, w, h);
		break;
		case STEP: // "posterisation"
			res = StepFilter(src, &dst, f, w, h, window->image->stat.min, window->image->stat.max, NULL);
		break;
		default: // differential filters
			res = DiffFilter(src, &dst, f, w, h);
	}
	if(!res){ g_err(_("Error in image filter module")); return; }
	window->context->transformMask |= TRANS_FILTER;
	_FREE(src);
	window->image->data = dst;
	window->context->current_image = window->image;
	compute_minmax(window->image);
	//gen_texture(window->image_transformed, window, TRUE);
	gen_texture(window->image, window, TRUE);
	force_redraw(window->drawingArea);
}


/*
 * Fill isolines' scale (an array)
 * Input:
 * 		f - filter for given method
 * 		min - minimum value of intensity
 * 		wd - max-min (dinamic range)
 * Output:
 * 		scale - a pointer to array (allocated in this function)
 */
int fillIsoScale(Filter *f, float **scale, float min, float wd){
	int M = f->w, y;
	float (*scalefn)(float in);
	float step, Nsteps = (float)f->w; // amount of intervals
	float Suniform(float in){
		return step*in + min;
	}
	float Slog(float in){
		return expf(in*step) + min - 1.;
	}
	float Sexp(float in){
		return wd*logf(in*step) + min;
	}
	float Ssqrt(float in){
		return in*in*step*step + min;
	}
	float Spow(float in){
		return sqrtf(in*step) + min;
	}
	if(!scale) return FALSE;

	*scale = calloc(M, sizeof(float));
	if(!*scale) return FALSE;
	switch(f->h){
		case LOG:
			scalefn = Slog; step = logf(wd+1.)/Nsteps;
		break;
		case EXP:
			scalefn = Sexp; step = expf(1.)/Nsteps;
		break;
		case SQRT:
			scalefn = Ssqrt; step = sqrtf(wd)/Nsteps;
		break;
		case POW:
			scalefn = Spow; step = wd*wd/Nsteps;
		break;
		default:
			scalefn = Suniform; step = wd/Nsteps;
	}
	for(y = 0; y < M; y++){
		(*scale)[y] = scalefn(y+1);
		DBG("level %d: I=%g", y, (*scale)[y]);
	}
	return TRUE;
}
