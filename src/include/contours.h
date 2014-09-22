#ifndef __CONTOURS_H__
#define __CONTOURS_H__
#include "fitsview.h"
#include "gtk.h"
#include "spots.h"
void find_contours(Window *window, int n, int type);
void free_contours(IMAGE *image);
int copy_contours(IMAGE *in, IMAGE *out, float dX, float dY);
#endif // __CONTOURS_H__
