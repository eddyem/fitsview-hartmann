#ifndef _FITSVIEW_H_
#define _FITSVIEW_H_

#ifndef _GNU_SOURCE
	#define _GNU_SOURCE
#endif

#include <unistd.h>
#include <err.h>
#include <stdlib.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <string.h>
#include <math.h>
#include <GL/gl.h>
#include <gtk/gtk.h>
#include <gtk/gtkgl.h>
#include <glib.h>
#include <glib/gprintf.h>
#include <libintl.h>  // for locale
#include <sys/time.h> // for dtime
#include <libgen.h>   // for basename
#include <fitsio.h>
#include <errno.h>
#include <limits.h>
#include <strings.h>
#include <stdarg.h> // for color text in tty


#include "CUtools.h"

#ifndef GETTEXT_PACKAGE
#define GETTEXT_PACKAGE "fitsview"
#endif
#define PRGNAME GETTEXT_PACKAGE
#ifndef LOCALEDIR
#define LOCALEDIR "/home/eddy/locale"
#endif

#define _(String)				gettext(String)
#define gettext_noop(String)	String
#define N_(String)				gettext_noop(String)

#ifndef THREAD_NUMBER
	#define THREAD_NUMBER 2
#endif
#ifndef OMP_NUM_THREADS
	#define OMP_NUM_THREADS THREAD_NUMBER
#endif

extern const int HIST_SIZE;
void *copy_struct(void *src, int size);
void *init_struct(int size);
#define INIT(var, type)		var = (type*)init_struct(sizeof(type))
#define COPY(src, type)		(type*)copy_struct((void*)src, sizeof(type))
#define _FREE(ptr)			do{free(ptr); ptr = NULL;}while(0)

/*
 * If the distance from the current position of the mouse to a certain point
 * less than this value, it is believed that this kind of chosen point
 */
#define PT_TRESHOLD		5

// debug mode, -DEBUG
#ifdef EBUG
	#define FNAME() fprintf(stderr, "\n%s (%s, line %d)\n", __func__, __FILE__, __LINE__)
	#define DBG(...) do{fprintf(stderr, "%s (%s, line %d): ", __func__, __FILE__, __LINE__); \
					fprintf(stderr, __VA_ARGS__);			\
					fprintf(stderr, "\n");} while(0)
#else
	#define FNAME()	 do{}while(0)
	#define DBG(...) do{}while(0)
#endif //EBUG

extern int (*red)(const char *fmt, ...);
extern int (*green)(const char *fmt, ...);

// contour poing
typedef struct _cPoint{
	float x;		// coordinates of the contour point in the CS picture
	float y;
	struct _cPoint *prev;// pointers to the previous and next points of the contour
	struct _cPoint *next;
} cPoint;
// contour structure
typedef struct _Contour{
	int N;				// amount of points
	cPoint *first;		// the first point of the contour (list)
	cPoint *last;		// last -//-
//	struct _Contour *prev; // pointers to the previous and following contour of the group
	struct _Contour *next;
	unsigned char closed;// == TRUE, if contour is closed
} Contour;
typedef struct{
	int L;	// level (at the lower boundary of 0, 1, etc)
	int N;	// number of units in the structure
	Contour *first; // pointer to the first contour in the list
	Contour *last;  // pointer to the last contour in the list
}cList;

// image histogram
typedef struct{
	int size;		// it is possible that the size of the histogram can be changed
	int max;		// maximum and minimum values of the histogram (the amplitude)
	int min;
	double bot_tres;	// upper and lower thresholds for equalization (normalized to 1)
	double top_tres;
	GLfloat scale;	// the scale - the range of intensities per unit of histogram
	int *data;
} ImageHistogram;
// image statistics
struct ImageStat{
	GLfloat max;
	GLfloat min;
	ImageHistogram histogram;
	/*GLfloat median;
	GLfloat mean;
	GLfloat std;*/
};

// image itself
typedef struct{
	int width;			// width
	int height;			// height
	float bZero;		// conditional zero point
	float bScale;		// scale
	float maxVal;		// the limiting value of the intensity (1<<BITTPIX)
	double val_f;		// the value of a focus in mm (if not defined - '-1e6')
	int dtype;			// data type
	struct ImageStat stat; // image statistica
	float *data;		// picture data
	GtkListStore *store;// list of options for each key
	gchar *filename;	// file name (basename) for the image
	int imagetype;		// type of image (for any changes)
	cList **contours;	// isolines array
	int Ncontours;		// dimension of the array of isolines
	float *cLevels;		// array of intensity levels corresponding to each contour
} IMAGE;
// imagetype
enum{
	 SIMPLE_IMAGE		// normal (not converted) image
	,HOUGH				// Hough transform
};

// the context
typedef struct{
	unsigned char drawingMode;	// drawing mode: cross section area, etc.
	unsigned char trackingMode;	// drawing mode sections
	unsigned char graphYaxis;	// form the Y-axis chart: logarithmic, or common
	unsigned char graphMode;	// chart mode (chart, histogram, etc.)
	unsigned char transformMask;// mask changes made to image_transformed
	unsigned char visualMode;	// display mode (additional areas, etc.)
	IMAGE *current_image;		// displayed image (window->image or window->image_transformed)
} Context;

typedef struct{
	unsigned char previewMode;	// preview settings
	gboolean add_all;			// try to add all the files (or only .fts / .fit [s])
	int minSpotW;				// the limit sizes of points
	int maxSpotW;
	int minSpotH;
	int maxSpotH;
} Global_context;

extern Global_context *Global;

/*
 * description fields of the structure of context
 * the first (zero) should go the default values
 */
// drawingMode:
enum{
	 DRAW_TRACK		// tracking mode
	,DRAW_SELECTION	// and so on
	,DRAW_QUAD
	,GAME_MODE
};
// trackingMode:
enum{
	 TRACK_ANY
	,TRACK_VERT
	,TRACK_HORIZ
	,TRACK_DIAG
};
// graphYaxis, may be masking: LSB - the graph, Sr. - histogram
enum{
	 Y_LINEAR  = 0
	,Y_LOGBOTH = 3
	,Y_LOGHIST = 2
	,Y_LOGGRAPH= 1
};
// graphMode - what is shown on graph (histogram or track)
enum{
	 GR_GRAPH
	,GR_HISTOGRAM
};
// transformMask - bit components
enum{
	 TRANS_NONE		= 0
	,TRANS_TRES		= 1
	,TRANS_ROT		= 2
	,TRANS_MOVE		= 4
	,TRANS_SCALE	= 8
	,TRANS_FILTER	= 16
	,TRANS_HARTMANN = 32
};
// visualMode - bit components
enum{
	 SHOW_ONLY_IMAGE	= 0
	,SHOW_POTBOXES		= 1
	,SHOW_POTSELECT		= 2
	,SHOW_ISOLINES		= 4
	,SHOW_ALL			= 255
};
// previewMode - bit components, masks
enum{
	 PREVIEW_DEFAULT			// (512x512, sqrt)
	,PREVIEW_128
	,PREVIEW_256
	,PREVIEW_512
	,PREVIEW_768
	,PREVIEW_1024
	,PREVIEW_SIZEMASK	= 7
	,PREVIEW_LINEAR		= 8
	,PREVIEW_LOG		= 16
	,PREVIEW_SQRT		= 24
	,PREVIEW_COLFNMASK	= 24
};

double dtime();

#endif // _FITSVIEW_H_
