#ifndef _IMTOOLS_H_
#define _IMTOOLS_H_
#include "fitsview.h"
#include "gtk.h"
#include "spots.h"

void destroy_image(IMAGE *image);
void free_contours(IMAGE *image);
gboolean make_histogram(IMAGE *ima);
void treshold_image(Window *window, GLfloat bottom, GLfloat top);
void get_spots_coords(Window *window, IMAGE *ima);
void hartmann_spots(Window *window);
gboolean get_gaussian_center(Window *window, BOX *box, Coordinates *crds);
PIX *conv_image2pix(IMAGE *ima);
void get_circles_params(Window *window, IMAGE *ima);
void hough_lines(Window *window);
void get_houg_line(Window *window, double x, double y, double *phi_, double *R_);
void filter_image(Window *window, Filter *f);
#endif // _IMTOOLS_H_
