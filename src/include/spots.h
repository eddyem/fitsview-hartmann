#ifndef _SPOTS_H_
#define _SPOTS_H_
#include "fitsview.h"

#ifndef LEPTONICA_FOUND
	#define LEPTERR g_err(_("Install Leptonica library and remake sources for this functions"))
	typedef struct{double x; double y; double w; double h;} BOX;
	typedef void* PIX;
#else
	#include <leptonica/allheaders.h>  // from leptonica
#endif // LEPTONICA_FOUND

/*
 * Coordinates of the center of area and its sizes was originally initiated
 * based on the box-structures, then adjusted to a coordinate system relative
 * to the optical axis, further refined by a Gaussian
 * (x and y - as the central parameters of Gaussians,
 * rx and ry - as the standard deviation of a Gaussian
 */
typedef struct{
	double x;
	double y;	// center coordinates
	double rx;
	double ry;	// half-sizes
} Coordinates;

// spot parameters
typedef struct{
	BOX	box;		// bounding box
	Coordinates c;	// dimensions and coordinates in CS relative to the optical axis, without rotation
	double r;		// module of the radius vector in the CS of the center of spots
	double phi;		// angular coordinate (from OX counter-clockwise)
	double xC;		// coordinates of spots in the CS of an optical axis, taking into account the rotation
	double yC;
	int id;			// number, id etc.
} Spot;

// spots array
struct Spots{
	int size;	// array size
	int n;		// amount of spots in array
	int memflag;// BY_COPY or BY_PTR - method of spots allocation
	Spot **spot;// spots themselves
};

typedef struct Spots Spots;

Spots* spots_alloc(int n);
Spot* spots_add(Spots *spots, Spot *spot, int copy);
enum {	// definition for the function spots_add
	BY_COPY,	// copy spot
	BY_PTR		// only insert pointer to a spot
};
void spots_free(Spots **spots);
struct Window;
void sort_spots(struct Window *window);
void spots_save(Spots *spots, gchar *filename);

extern double AxisX, AxisY;

#endif // _SPOTS_H_
