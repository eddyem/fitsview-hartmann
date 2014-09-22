#ifndef _TRACKING_H_
#define _TRACKING_H_

#include "fitsview.h"
#include "gtk.h"

typedef struct{
	double x;
	double y;
} Point;

typedef struct{
	int n;
	Point *data;
} Points;

extern double gYmax, gXmax;
extern double xScale, yScale;

void do_tracking(double x, double y, gboolean start, Window *window);
void do_histogram(Window *window);
gboolean graph_mouse_btn(int x, Window *window);
gboolean isYlog(Window *window);

void gen_tres_texture(Window *window);
gboolean Gexpose(Window *window);
gboolean Gconfigure(Window *window);
gboolean fit_gaussian(Window *window);
void get_image_crds(double x, double *xx, double *yy);
void get_sel_region(double *x1, double *y1, double *x2, double *y2);
#endif // _TRACKING_H_
