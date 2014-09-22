//      tracking.c - funcions for cutting tracking
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
#include "tracking.h"
#include "imtools.h"
#include "gauss.h"
#include "opengl.h"

const double LHWD = 5.; // length of the cross-mark bar
double x_1, y_1, x_2, y_2; // coordinates the selection frame or line
//double GW, GH;	// height and width of the graph

double gY0, gYm, gXm; // range of values of Y, the maximum value of X
double gYmax;  // limits to indicate the coordinates in the status bar
double xScale, yScale; // scale by axes
double bot_X = -1., top_X = -1.;  // values (in units of X) of the lower and upper bounds on the abscissa
gint  g_botX, g_topX; // the same, but in the coordinates of the window

/* the data for translation of the x coordinate from Pdata into coordinates x, y of image
 * X = x0 + x*cx, Y= y0 + x*cy (cx - the cosine of the angle of the cutting, cy - sine) */
struct{
	double x0;	// starting poing
	double y0;
	double cx;	// cosine of the angle of inclination
	double cy;	// sine -//-
} Pdata2image;
Points *Pdata = NULL; // array of data points
GdkPoint *gPoints = NULL;	// points in window coordinates
Points selection = {0, NULL}; // the same for the selected region
GdkPoint *selPoints = NULL;

gint Npts;					// number of data points

Point *alloc_Pdata(int N){
	if(Pdata){
		free(Pdata->data);
		Pdata->data = NULL;
	}
	else
		Pdata = malloc(sizeof(Pdata));
	if(!Pdata){
		g_err(_("Can't allocate memory for track points"));
		return NULL;
	}
	Pdata->data = malloc(N*sizeof(Point));
	if(!Pdata->data){
		g_err(_("Can't allocate memory for track points"));
		return NULL;
	}
	return Pdata->data;
}
// returns the coordinates xx and yy on the image (in CS of the window!)
// for the point with coordinate x in the graph
void get_image_crds(double x, double *xx, double *yy){
	*xx = Pdata2image.x0 + x * Pdata2image.cx;
	*yy = Pdata2image.y0 + x * Pdata2image.cy;
}

// reset the sample
void clear_selection(){
	free(selection.data);
	selection.data = NULL;
	selection.n = 0;
}
/*
 * filling an array Pdata based on the values of coordinates initial (X0, Y0)
 * and final (X, Y) points
 */
void get_data(double X0, double X, double Y0, double Y, Window *window){
	int n = 0, N; // new array size
	IMAGE *image = window->image;
	int x0, y0, yi, WD = image->width, HT = image->height;
	double x, y, A, xstep, ystep;
	Point *ptr;
	gboolean Vert = FALSE;
	double tgXY, cosXY;
	inline void minmax(double amp){
		if(gY0 > amp) gY0 = amp;
		else if(gYm < amp) gYm = amp;
	}
	inline double dist(double dx, double dy){
		return sqrtf(dx*dx + dy*dy);
	}
	gXm = dist(X-X0, Y-Y0);
	N = (int)(gXm + 1.);
	if(N < 2) return;
	if(!image->data)
		return;
	if( X0 < 0. || X0 > WD || X < 0. || X > WD ||
		Y0 < 0. || Y0 > HT || Y < 0. || Y > HT) return;
	if(!(ptr = alloc_Pdata(N))) return;
	Pdata2image.x0 = X0;
	Pdata2image.y0 = Y0; //HT - 1. - Y0;
	Pdata2image.cx = (X-X0)/gXm;
	Pdata2image.cy = (Y-Y0)/gXm;
	//DBG("N=%d, ", N);
	x0 = (int)(X0 + 0.5);
	y0 = (int)(Y0 + 0.5);
	if(fabs(X-X0) < fabs(Y-Y0)){ // closer to the vertical
		Vert = TRUE;
		tgXY = (X-X0) / (Y-Y0);
		cosXY = fabs(cos(atan(tgXY)));
	}
	else{
		tgXY = (Y-Y0) / (X-X0);
		cosXY = fabs(cos(atan(tgXY)));
	}
	yi = x0 + WD * y0;
	gY0 = gYm = image->data[yi];// * wd + min;
	if(Vert){
		ystep = (Y > Y0) ? cosXY : -cosXY;
		for(y = Y0; ; y+=ystep, n++, ptr++){
			if(n > N) break;
			if(ystep > 0.){
				if(y > Y) break;
			}
			else{
				if(y < Y) break;
			}
			ptr->x = (double)n;
			x = X0 + tgXY * (y - Y0);
			//yi = (int)(x+0.5) + WD * (HT - (int)(y+0.5));
			yi = (int)(x+0.5) + WD * ((int)(y+0.5));
			A = image->data[yi];// * wd + min;
			ptr->y = A;
			minmax(A);
		}
	}
	else{
		xstep = (X > X0)? cosXY : -cosXY;
		for(x = X0; ; x+=xstep, n++, ptr++){
			if(n > N) break;
			if(xstep > 0.){
				if(x > X) break;
			}
			else{
				if(x < X) break;
			}
			y = Y0 + tgXY * (x-X0);
			ptr->x = (double)n;
			//yi = (int)(x+0.5) + WD * (HT - (int)(y+0.5));
			yi = (int)(x+0.5) + WD * ((int)(y+0.5));
			A = image->data[yi];// * wd + min;
			ptr->y = A;
			minmax(A);
		}
	}
	//DBG("Ymin=%g, Ymax=%g\n", gY0, gYm);
	Pdata->n = n - 1;
}

void draw_cross(GdkDrawable *d, GdkGC *gc, double x, double y){
	gdk_draw_line(d, gc, x-LHWD,y, x+LHWD, y);
	gdk_draw_line(d, gc, x,y-LHWD, x, y+LHWD);
}
// draw rectangle (gdk_draw_rectangle don't match because it require x,y > x0,y0
void draw_rect(GdkDrawable *d, GdkGC *gc, double x0, double y0, double x, double y){
	gdk_draw_line(d, gc, x0,y0, x0,y);
	gdk_draw_line(d, gc, x0,y0, x,y0);
	gdk_draw_line(d, gc, x,y, x, y0);
	gdk_draw_line(d, gc, x,y, x0, y);
}

// get coordinates of selected region
void get_sel_region(double *x0, double *y0, double *x, double *y){
	*x0 = x_1; *y0 = y_1; *x = x_2; *y = y_2;
}
void do_tracking(double xx, double yy, gboolean start, Window *window){
	static double xold, yold, X0, Y0;
	double X, Y;
	double imW = (double)window->texture->w;
	double imH = (double)window->texture->h;
	Context *c = window->context;
	gboolean conv = FALSE;
	x_2 = xx; y_2 = yy;
	conv_mouse_to_image_coords(x_2, y_2, &X, &Y, window);
	//FNAME();
	if(!window->image->data) return;
	if(X > imW){ X = imW; conv = TRUE;}
	else if(X < 0){ X = 0; conv = TRUE;}
	if(Y > imH){ Y = imH; conv = TRUE;}
	else if(Y < 0){ Y = 0; conv = TRUE;}
	if(conv) conv_image_to_mouse_coords(X, Y, &x_2, &y_2, window);
	//DBG("X=%g, Y=%g, x=%g, y=%g", X, Y, x_2, y_2);
	if(start){
		clear_selection();
		X0 = X; Y0 = Y;
		x_1 = xold = x_2;
		y_1 = yold = y_2;
		return;
	}
	if(c->trackingMode == TRACK_VERT){
		X = X0; x_2 = x_1;
	}
	else if(c->trackingMode == TRACK_HORIZ){
		Y = Y0; y_2 = y_1;
	}
	else if(c->trackingMode == TRACK_DIAG){
		double s, adx, ady;
		gboolean changed = FALSE;
		adx = fabs(X-X0); ady = fabs(Y-Y0);
		void corrY(){
			s = (Y > Y0) ? 1. : -1.;
			Y = Y0 + s * adx;
			y_2 = y_1 + s * fabs(x_2-x_1);
			if(Y > imH || Y < 0.){
				if(Y > imH) Y = imH;
				else Y = 0.;
				conv_image_to_mouse_coords(X, Y, &x_2, &y_2, window);
				changed = TRUE;
			}
			/*
			if(Y > imH){ Y = imH; y = Y/Daspect + Dy; changed = TRUE;}
			else if(Y < 0){ Y = 0; y = Dy; changed = TRUE;}*/
			ady = fabs(Y-Y0);
			//DBG("Y=%g, Y0=%g", Y, Y0);
		}
		void corrX(){
			s = (X > X0) ? 1. : -1;
			X = X0 + s * ady;
			x_2 = x_1 + s * fabs(y_2-y_1);
			if(X > imW || X < 0.){
				if(X > imW) X = imW;
				else X = 0.;
				conv_image_to_mouse_coords(X, Y, &x_2, &y_2, window);
				changed = TRUE;
			}
			adx = fabs(X-X0);
			//DBG("X=%g, X0=%g", X, X0);
		}
		if(adx > ady){
			corrY();
			if(changed) corrX();
		}
		else{
			corrX();
			if(changed) corrY();
		}
	}

	GdkDrawable *d = window->drawingArea->window;
	static GdkGC *gc = NULL;
	if(!gc){
		gc = gdk_gc_new(d);
		gdk_gc_set_foreground(gc, &(window->drawingArea->style->white));
		gdk_gc_set_function(gc, GDK_XOR);
	}
	if(c->drawingMode == DRAW_QUAD){// draw rectangle
		draw_rect(d, gc, x_1,y_1, xold,yold);
		draw_rect(d, gc, x_1,y_1, x_2,y_2);
	}
	else{ // draw lint
		// first erase old one
		gdk_draw_line(d, gc, x_1,y_1, xold,yold);
		draw_cross(d, gc, xold, yold);
		// then draw new one
		gdk_draw_line(d, gc, x_1,y_1, x_2,y_2);
		draw_cross(d, gc, x_2, y_2);
		//DBG("Xima=%g, Yima=%g (X0=%g, Y0=%g)", X, Y, X0, Y0);
		get_data(X0, X, Y0, Y, window);
		bot_X = 0.; top_X = gXm;
		//DBG("X=%g, X0=%g, Y=%g, Y0=%g", X, X0, Y, Y0);
		gchar str[80];
		g_snprintf(str, 79, "length=%.1f, angle=%.1fdegr",
					sqrt((X-X0)*(X-X0)+(Y-Y0)*(Y-Y0)),
					atan2(Y-Y0,X-X0)/M_PI*180.);
		set_status_text(StatusText, str, window->parent);
		set_status_text(StatusText, "", window);
		Gconfigure(window->graphWindow);
		gdk_window_invalidate_rect(window->graphWindow->drawingArea->window,
			NULL, TRUE);
	}
	xold = x_2; yold = y_2;
}

gboolean Gexpose(Window *window){//, GdkEventExpose *event, gpointer userData){
	cairo_t *cr = NULL;
	int i;
	GtkWidget *Area = window->drawingArea;
	int H = (int)(Area->allocation.height + 0.5);
	GdkPoint *gptr = &gPoints[1];
	if(!cr){
		cr = gdk_cairo_create(Area->window);
		cairo_set_line_width (cr, 1.);
	}
	cairo_set_source_rgb(cr, 0., 0., 0.);
	if(!Pdata || !gPoints)
		return FALSE;
	if(!Pdata->data || Pdata->n < 2)
		return FALSE;
	cairo_move_to(cr, gPoints->x, gPoints->y);
	for(i = 1; i < Npts; i++, gptr++)
			cairo_line_to(cr, gptr->x, gptr->y);
	cairo_stroke(cr);
	cairo_set_source_rgb(cr, 0., 1., 0.);
	cairo_move_to(cr, g_botX, 0);
	cairo_line_to(cr, g_botX, H);
	cairo_stroke(cr);
	cairo_set_source_rgb(cr, 1., 0., 0.);
	cairo_move_to(cr, g_topX, 0);
	cairo_line_to(cr, g_topX, H);
	cairo_stroke(cr);
	if(selPoints && selection.n){
		int n = selection.n;
		cairo_set_source_rgb(cr, 0., 0., 1.);
		cairo_move_to(cr, selPoints->x, selPoints->y);
		gptr = &selPoints[1];
		//DBG("n=%d", n);
		for(i = 1; i < n; i++, gptr++)
			cairo_line_to(cr, gptr->x, gptr->y);
		cairo_stroke(cr);
	}
	cairo_destroy(cr);
	return FALSE;
}

// calculation of the chart points in window coordinates
gboolean Gconfigure(Window *window){//, GdkEventConfigure *event, gpointer none){
	int i, n;
	Point *ptr;
	//FNAME();
	double W, H, ym, y0;
	GdkPoint *gptr;
	gboolean log_scale = isYlog(window);
	if(!Pdata)
		return FALSE;
	if(!Pdata->data || Pdata->n < 2)
		return FALSE;
	GtkWidget *Area = window->drawingArea;
	free(gPoints);
	gPoints = NULL;
	Npts = Pdata->n;
	gPoints = malloc((Npts + 1)* sizeof(GdkPoint));
	if(!gPoints){
		g_err(_("Can't allocate memory for track points"));
		return FALSE;
	}
	ptr = Pdata->data;
	if(!ptr) return FALSE;
	gptr = gPoints;
	W = Area->allocation.width;
	H = Area->allocation.height;
	ym = gYm; y0 = gY0;
	if(log_scale){
		if(y0 < 0.) y0 = -log(-y0 + 1.);
		else y0 = log(y0 + 1.);
		if(ym < 0.) ym = -log(-ym + 1.);
		else ym = log(ym + 1.);
	}
	xScale = W / gXm;
	yScale = H / (ym - y0);
	gYmax = ym;
	for(i = 0; i < Npts; i++, ptr++, gptr++){
		gptr->x = (int)(ptr->x * xScale + 0.5);
		if(log_scale){
			double yy = ptr->y;
			if(yy < 0.)
				yy = -log(-yy + 1.);
			else
				yy = log(yy + 1.);
			gptr->y = (int)(H - (yy - y0) * yScale + 0.5);
		}
		else
			gptr->y = (int)(H - (ptr->y - y0) * yScale + 0.5);
	}
	g_botX = (int)(bot_X * xScale + 0.5);
	g_topX = (int)(top_X * xScale + 0.5);
	// Long live code monkeys of all time!
	if(selection.data) do{
		n = selection.n;
		free(selPoints);
		selPoints = malloc((n+1)* sizeof(GdkPoint));
		if(!selPoints) break;
		gptr = selPoints;
		ptr = selection.data;
		for(i = 0; i < n; i++, ptr++, gptr++){
			gptr->x = (int)(ptr->x * xScale + 0.5);
			if(log_scale){
				double yy = ptr->y;
				if(yy < 0.)
					yy = -log(-yy + 1.);
				else
					yy = log(yy + 1.);
				gptr->y = (int)(H - (yy - y0) * yScale + 0.5);
			}
			else
				gptr->y = (int)(H - (ptr->y - y0) * yScale + 0.5);
//			DBG("g: (%g, %g), w: (%d, %d)", ptr->x, ptr->y, gptr->x, gptr->y);
		}
	}while(0);
	/*
	if(ym > 1000.)
		set_Grulers(y0/1000, gXm, ym/1000);
	else*/
		set_Grulers(y0, gXm, ym, window);
	return FALSE;
}

// display mode of the histogram in the field of graph
void do_histogram(Window *window){
	FNAME();
	IMAGE *image = window->context->current_image;
	if(!image || !image->data) image = window->image;
	if(!image || !image->data) return;
	int i, hsize = image->stat.histogram.size;
	int *hptr = image->stat.histogram.data;
	Point *ptr;
	double histsiz = (double)image->stat.histogram.size;
	double scale = image->stat.histogram.scale * (double)image->stat.histogram.size;
	//double imin = image->stat.min;
	if(!(ptr = alloc_Pdata(hsize))) return;
	gXm = 1.;//scale+imin;
	scale /= histsiz;
	gYm = (double)(image->stat.histogram.max);
	gY0 = (double)(image->stat.histogram.min);
	for(i = 0; i < hsize; i++, ptr++, hptr++){
		ptr->x = (double) i / histsiz;
		//ptr->x = (double) i * scale + imin;
		ptr->y = (double) *hptr;
	}
	Pdata->n = hsize;
	DBG("bot: %g, top: %g", image->stat.histogram.bot_tres, image->stat.histogram.top_tres);
	if(image->stat.histogram.bot_tres > 0.)
		bot_X = image->stat.histogram.bot_tres;
	else
		bot_X = 0.;
	if(image->stat.histogram.top_tres > 0.)
		top_X = image->stat.histogram.top_tres;
	else
		top_X = gXm;
	Gconfigure(window->graphWindow);
	gdk_window_invalidate_rect(window->graphWindow->drawingArea->window, NULL, TRUE);
}

/*
 * will try to move boundaries when clicking near graph borders
 * when clicked near a curve - something else can be done ...
 * Window - graphics window
 * Return value: TRUE, if will need to handle mouse motion events further
 */
gboolean graph_mouse_btn(int x, Window *window){
	gboolean redraw = FALSE;
	gboolean ret = FALSE;
	double W, xScale, tmp;
	GtkWidget *graphArea = window->drawingArea;
	int dist = (g_topX - g_botX) / 2; // half of a distance between borders
	if(dist < PT_TRESHOLD) dist = PT_TRESHOLD;
	W = graphArea->allocation.width;
	xScale = W / gXm;
	if(abs(x - g_botX) < dist || x < g_botX){ // select down border
		tmp = (double) x / xScale;
		if(tmp < top_X) bot_X = tmp;
		if(bot_X < 0.) bot_X = 0.;
		ret = TRUE;
		redraw = TRUE;
	}
	else if(abs(x - g_topX) < dist || x > g_topX){ // select up border
		tmp = (double) x / xScale;
		if(tmp > bot_X) top_X = tmp;
		if(top_X > gXm) top_X = gXm;
		ret = TRUE;
		redraw = TRUE;
	}
	if(redraw){
		gchar msg[80];
		//DBG("move: botX=%g, topX=%g", bot_X, top_X);
		g_snprintf(msg, 79,
			_("Selected region from %.1f to %.1f, %d datapoints"),
			bot_X, top_X, (int)(((double)Npts)/gXm*(top_X-bot_X)+0.5));
		set_status_text(StatusText, msg, window);
		Gconfigure(window);
		gdk_window_invalidate_rect(graphArea->window, NULL, TRUE);
	}
	return ret;
}

void gen_tres_texture(Window *window){
	window->image->stat.histogram.bot_tres = bot_X;
	window->image->stat.histogram.top_tres = top_X;
	treshold_image(window, bot_X / gXm, top_X / gXm);
}

gboolean isYlog(Window *window){
	gboolean log_scale = FALSE;
	if(window->context->graphMode == GR_GRAPH){
		if(window->context->graphYaxis & Y_LOGGRAPH)
			log_scale = TRUE;
	} else if(window->context->graphYaxis & Y_LOGHIST)
			log_scale = TRUE;
	return log_scale;
}

gboolean fit_gaussian(Window *window){
	//FNAME();
	Point *pt, *pt0;
	double C, A, s, x0;
	double xx, yy;
	#ifndef GSL_FOUND
		GSLERR;
		return FALSE;
	#endif
	if(!window->parent || !window->parent->texture){
		g_err(_("No parent window or image"));
		return FALSE;
	}
	double imW = (double)window->parent->texture->w;
	double imH = (double)window->parent->texture->h;
	gchar msg[160];
	GtkWidget *graphArea = window->drawingArea;
	if(!Pdata){ // no points? (O_o, but such a thing be?)
		g_err(_("No data points"));
		return FALSE;
	}
	int n = (int)(((double)Npts)/gXm*(top_X-bot_X) + 0.5);
	if(n < 10 || n > (int)(sqrt(imW * imW + imH * imH))){ // few or many points
		g_err(_("Bad amount of data points"));
		return FALSE;
	}
	int n0 = (int)(((double)Npts)/gXm*bot_X); // the number of a first point
	if(n0 < 0 || (n0 + n) > Npts){ // do not get there?
		g_err(_("Invalid number of first point"));
		return FALSE;
	}
	clear_selection();
	selection.n = n;
	selection.data = malloc(n * sizeof(Point));
	if(!selection.data){
		g_err(_("Can't allocate memory"));
		return FALSE;
	}
	pt = selection.data;
	pt0 = &Pdata->data[n0];
	DBG("pt0: (%g, %g)", pt0->x, pt0->y);
	// Copy the desired sample into a separate array
	for(; n--; pt++, pt0++){
		pt->x = pt0->x;
		pt->y = pt0->y;
	}
	gauss_fit(&selection, &C, &A, &s, &x0);
	gaussian_v(C, A, s, x0, &selection);
	get_image_crds(x0, &xx, &yy);
	//DBG("move: botX=%g, topX=%g", bot_X, top_X);
	g_snprintf(msg, 159,
		_("Fit gaussian, x0=%.2f I(%.2f, %.2f), s=%.2f, A=%.1f, C=%.3f"),
		x0, xx, yy, s, A, C);
	set_status_text(StatusText, msg, window);
	Gconfigure(window);
	gdk_window_invalidate_rect(graphArea->window, NULL, TRUE);
	return TRUE;
}
