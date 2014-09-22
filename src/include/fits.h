#ifndef _FITS_H_
#define _FITS_H_
#include "fitsview.h"

#define SCALE_MIN (0.01)		// scale smaller than this is wrong
gboolean cmplx_conv(gchar *value, double *re, double *im);
gboolean readfits(gchar *filename, IMAGE *data);
gboolean writefits(gchar *filename, IMAGE *data);
gboolean try_open_file(gchar *filename, IMAGE *data);
gboolean guess_fits(gchar *filename);
#endif // _FITS_H_
