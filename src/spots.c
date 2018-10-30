//      spots.c - functions to work with spots (sorting, enumerating & so on)
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

/*
 * Since our image is mirrored, to establish the direction of the axis
 * of symmetry it is necessary to know the position angle of the camera, etc.
 * It is easier to work in a "mirror" coordinates. Similarly, the coordinates
 * of the mask should also be given the "mirror".
 */
#include "gtk.h"
#include "spots.h"
#include "gauss.h"

const double M_2PI = 2. * M_PI;
// about a quarter of 11.25 degrees - the threshold for the identification of radial points and frames
const double PHI_STEP = 0.05;
// angle between the markers (for identification), 56.25 degrees
const double DPHI_MARKERS = 0.98;
const double aStep = M_PI / 16.; // 11.25 degrees, the angle between adjacent rays
const int N_ANGLES = 16; // number of lines on hartmanogramm (half of number of rays)
const int N_RINGS = 8;	// number of rings on hartmanogramm

// angle of the Y-axis of the matrix relative to the Y-axis mask
double SlopeAngle;
// angular position of the marker (near the center) for a zero slope, -101.25 degrees
const double Slope0 = -1.767146;

/*
 * these constants are the position of the axis of rotation of the table.
 * In the coordinate system of the image (ie, the Y axis goes from top to bottom)
 * means that it coincides with the optical axis. Naturally, in the work
 * program these coordinates should be calculated
 */
double AxisX = 1523., AxisY = 1533.;


// declaration of the flags of errors leading to the exit from auto mode
enum{
	 TOO_MANY_MARKERS			// discovered more than two markers
	,NOT_ENOUGH_MARKERS			// not found the two markers
	,BAD_MARKERS_ANGLE			// angle between the markers over DPHI_MARKERS +- PHI_STEP
	,TOO_MANY_RINGS				// found more rings than it should be
};

// messages appropriate to flags
char *errmsg[] = {
	 N_("Too many markers")
	,N_("Not enough markers")
	,N_("Bad markers angle")
	,N_("Too many rings")
};

#ifdef LEPTONICA_FOUND
/*
 * Violated any restrictions on the automatic mode
 * In theory the function should exit the batch mode with the mapping errors
 */
void bad_circumstances(int errcode){
	red("\nerror, need to stop\n");
	g_err(_(errmsg[errcode]));
}

// memory allocation for n spots array
Spots* spots_alloc(int n){
	Spots *ptr = malloc(sizeof(Spots));
	if(!ptr){
		g_err(_("Can't allocate memory for new spots"));
		return NULL;
	}
	ptr->spot = malloc(n * sizeof(Spot*));
	if(!ptr->spot){
		g_err(_("Can't allocate memory for pointers in spots array"));
		free(ptr);
		return NULL;
	}
	ptr->size = n;
	ptr->n = 0;
	// originally put the flag BY_COPY, which is removed by the first copy BY_PTR,
	// order to do not try to free one area several times
	ptr->memflag = BY_COPY;
	return ptr;
}

void spots_free(Spots **spots){
	FNAME();
	int i, n;
	if(!(*spots)) return;
	n = (*spots)->n;
	Spot **spot = (*spots)->spot;
	if((*spots)->memflag != BY_PTR) // check whether spot-regions were allocated
		for(i = 0; i < n; i++, spot++)
			free(*spot);
	free((*spots)->spot);
	free(*spots);
	*spots = NULL;
}

/*
 * Add point spot into array spots by copy.
 * If array is full return NULL
 *
 */
Spot* spots_add(Spots *spots, Spot *spot, int copy){
	Spot *ptr = NULL;
	if(spots->size <= spots->n){
		g_err(_("Spots array is full, can't place new spot\n"));
		return NULL; // there's no enough space
	}
	if(!spot){
		g_err(_("Zero pointer to spot\n"));
		return NULL;
	}
	if(copy == BY_COPY){
		ptr = malloc(sizeof(Spot));
		if(!ptr){
			g_err(_("Can't allocate memory for new spot"));
			return NULL;
		}
		if(!memcpy(ptr, spot, sizeof(Spot))){
			free(ptr);
			g_err(_("Can't copy spot to new one"));
			return NULL;
		}
	}
	else if(copy == BY_PTR){
		spots->memflag = BY_PTR; // change flag to not do free()
		ptr = spot;
	}
	spots->spot[spots->n] = ptr; // jam into array
	spots->n++;	// increment points number
	return ptr;
}

/*
 * alignment of the coordinate system spot.c to CS relative to the optical axis,
 * approximate calculation of the angular and radial coordinates of each spot in
 * the coordinates with respect to the center of gravity of spots (for further identification)
 */
void spots_center(Spots *spots){
	int i, n;
	double X = 0., Y = 0., N, x0, y0;
	double *xx, *yy, *px, *py;
	Spot **spot;
	n = spots->n; N = (double) n;
	if(n == 0) return;
	xx = malloc(n * sizeof(double)); px = xx;
	yy = malloc(n * sizeof(double)); py = yy;
	spot = spots->spot;
	for(i = 0; i < n; i++, px++, py++, spot++){
				//x0 = (*spot)->c.x; y0 = (*spot)->c.y;
		x0 = (*spot)->c.x - AxisX;
		//y0 = (*spot)->c.y - AxisY;
		y0 = AxisY - (*spot)->c.y; // flip into right coordinates
		(*spot)->c.x = x0; (*spot)->c.y = y0;
		*px = x0;
		*py = y0;
		X += *px; Y += *py;
	}
	X /= N; Y /= N;
	green("\nspots center: %g, %g", X, Y);
	px = xx; py = yy;
	spot = spots->spot;
	for(i = 0; i < n; i++, px++, py++, spot++){
		x0 = *px - X; y0 = *py - Y;
		(*spot)->phi = atan2(y0, x0);
		(*spot)->r = sqrt(x0*x0 + y0*y0);
	}
	free(xx); free(yy);
}

int cmp_spots_r(const void *spot1, const void *spot2){
	return ((*(Spot**)spot1)->r > (*(Spot**)spot2)->r);
}
void sort_spots_by_r(Spots *spots){
	qsort(spots->spot, spots->n, sizeof(Spot*), cmp_spots_r);
}

int cmp_spots_phi(const void *spot1, const void *spot2){
	return ((*(Spot**)spot1)->phi > (*(Spot**)spot2)->phi);
}
void sort_spots_by_phi(Spots *spots){
	qsort(spots->spot, spots->n, sizeof(Spot*), cmp_spots_phi);
}

int cmp_spots_id(const void *spot1, const void *spot2){
	return ((*(Spot**)spot1)->id > (*(Spot**)spot2)->id);
}
void sort_spots_by_id(Spots *spots){
	qsort(spots->spot, spots->n, sizeof(Spot*), cmp_spots_id);
}


/*
 * Sorting the points in R and phi we find the markers.
 * After sorting, id of markers is set to <Beam number> + 100 * <ring number>.
 * The markers are numbered as 100 * <ring number> + <number of the next clockwise beam>
*/
void sort_spots(Window *window){
	FNAME();
	Spots *spots = window->spots;
	gboolean Prefocal = TRUE;
	int i, n;
	int	step = 0,	// a step - goto next group
		num = 0,	// number of ray or ring
		mZero = 0;	// id of zero ray
	Spot **spot = spots->spot;
	Spots *markers = spots_alloc(2); // allocate space for markers
	double *angles_pre, *angles_post; // angle coordinate of each ray in "normal" orientation
	const int N_RAY_MAX = 2 * N_ANGLES - 1; // max ray number
	const int N_RAYS = 2 * N_ANGLES;
	const int N_ANGLE_MAX = N_ANGLES - 1; // max line number
	double r_aver = 0.; // mid radius of spots
	n = spots->n;
	// fill the array angles_0 (for obtaining further angle of inclination of hartmanogramm)
	angles_pre = calloc(N_ANGLES, sizeof(double));
	angles_post = calloc(N_ANGLES, sizeof(double));
	{
		angles_post[0] = M_PI_2 - aStep/2.;
		angles_pre[0] = M_PI_2 + aStep/2.;
		for(i = 1; i < N_ANGLES; i++){
			angles_pre[i] = angles_pre[i-1] + aStep;
			angles_post[i] = angles_post[i-1] - aStep;
			green("angles[%d]: pre=%g, post=%g",
				i, angles_pre[i]*180/M_PI, angles_post[i]*180/M_PI);
		}
	}
	// find the center of gravity of points on hartmanogramm and determine
	// for each point the polar coordinates in the system of the center of gravity
	spots_center(spots);
	// sort points by phi to find markers and number by rays
	{
		sort_spots_by_phi(spots);
		double phi0 = (*spot)->phi, phi, dphi;
		green("\nSort by phi\n");
		for(i = 0; i < n; i++, spot++){
			phi = (*spot)->phi;
			dphi = fabs(phi-phi0);
			phi0 = phi;
			// If we two consecutive jump occurs, the previous point has to be a marker.
			// If we we find two more markers, generate an error
			if(dphi > PHI_STEP){
				num++;
				if(!step){
					step = 1;
				}else{
					num--;
					red("FOUND MARKER");
					(*(spot-1))->id = 99 - (*(spot-1))->id;
					DBG("spot[%d]: r=%g, phi=%g\n", (*(spot-1))->id,
						(*(spot-1))->r, (*(spot-1))->phi);
					if(!spots_add(markers, *(spot-1), BY_PTR)){
						bad_circumstances(TOO_MANY_MARKERS);
						goto ret;
					}
				}
			}
			else step = 0;
			(*spot)->id = num;
			//DBG("spot[%d]: r=%g, phi=%g, dphi=%g\n", (*spot)->id, (*spot)->r,
			//	phi, dphi);
		}
	}
	// check markers
	{
		if(markers->n != 2){
			if(markers->n == 0){
				bad_circumstances(NOT_ENOUGH_MARKERS);
				goto ret; // no marker found
			}
			spot = spots->spot;
			for(i = 0; i < n; i++, spot++)
				if((*spot)->id < N_RAYS) r_aver += (*spot)->r;
			r_aver /= (double)n;
			green("R_aver = %.2f, R=%.2f", r_aver, markers->spot[0]->r);
			if(window->image->val_f < 0.){
				// obtain the value of Prefocal from human as VAL_F is unknown
				get_prefocal(window, &Prefocal);
			}else Prefocal = (window->image->val_f < 100.) ? TRUE : FALSE;
			int m_id = markers->spot[0]->id;
			if(!Prefocal){
				if(markers->spot[0]->r < r_aver){ // nearest, zero marker
					mZero = 99 - m_id;
				}else{ // far, fifth marker
					mZero = 104 - m_id; // 5 rays clockwise
				}
			}else{
				if(markers->spot[0]->r < r_aver){
					mZero = 98 - m_id;
				}else{
					mZero = 93 - m_id;
				}
			}
			if(markers->spot[0]->r < r_aver)
				markers->spot[0]->id = 299; // finally change markers' numbers
			else
				markers->spot[0]->id = 794;
			if(mZero > N_RAY_MAX) mZero -= N_RAYS;
			else if(mZero < 0) mZero += N_RAYS;
		}else{
			sort_spots_by_r(markers);
			spot = markers->spot;
			for(i = 0; i < 2; i++)
				green("marker #%d: r=%g, phi=%g, tmp_id=%d, x=%g, y=%g", i,
					spot[i]->r, spot[i]->phi, spot[i]->id, spot[i]->c.x, spot[i]->c.y);
			double dphi = spot[0]->phi - spot[1]->phi;
			if(dphi < 0.) dphi += 2.*M_PI;
			DBG("Angle between markers: %g degr\n", dphi*180./M_PI);
			if(fabs(dphi - DPHI_MARKERS) > PHI_STEP){
				bad_circumstances(BAD_MARKERS_ANGLE);
				goto ret;
			}
			if(dphi > 0.){
				Prefocal = FALSE;
				mZero = 99 - spot[0]->id; // number of ray to be zero
			}else{
				mZero = 98 - spot[0]->id; // number of ray to be zero
				if(mZero > N_RAY_MAX) mZero -= N_RAYS;
			}
			spot[0]->id = 299;
			spot[1]->id = 794;
		}
	}
	// renumbers points in a clockwise direction, starting from the scheme
	// EEERRRRRRROOOOOOORRRRR!!!!
	// смены направления отсчетов нет: и на пред- и на зафокальном одинаково!
	// The numbering starts from zero!
	spot = spots->spot;
	green("\nmZero = %d\n", mZero);
	if(!Prefocal){ // postfocal image
		red("\nPostfocal image\n");
		for(i = 0; i < n; i++, spot++){
			if((*spot)->id > N_RAY_MAX) continue; // don't convert markers' numbers
			int sid = mZero - (*spot)->id; // change the numbering and its direction
			if(sid < 0) sid += N_RAYS;
			(*spot)->id = sid;
		}
	}else{ // prefocal image
		red("\nPrefocal image\n");
		for(i = 0; i < n; i++, spot++){
			if((*spot)->id > N_RAY_MAX) continue;
			int sid = (*spot)->id - mZero; // only change the numbering
			if(sid < 0) sid += N_RAYS;
			(*spot)->id = sid;
		}
	}
	/* by the method of least squares we are counting angle of inclination
	 * y=a+bx, b=tg(slope+alpha_i0)
	 * b = [\sum x_i \sum y_i - n\sum (x_i y_i)] /
	 *        [ (\sum x_i)^2 - n\sum (x_i)^2 ]
	 *
	 * Let S1 = \sum x_i, S2 = \sum y_i, S3 = \sum(x_i)^2,
	 * 		S4 = \sum(x_i y_i), then
	 * b = (S1*S2 - n*S4) / (S1*S1 - n*S3)
	 *
	 * Set up arrays four elements Sj of 16 (for each line by an item)
	 * In the loop over all spots (except markers) fill the array.
	 * Then we are counting the coefficients b and fill the matrix of slope angles
	 */
	{
		double	*S1, *S2, *S3, *S4, *N, *angles;
		double	x, y, *angles_0;
		int		idx;
		angles_0 = (Prefocal) ? angles_pre : angles_post;
		S1 = calloc(N_ANGLES, sizeof(double));
		S2 = calloc(N_ANGLES, sizeof(double));
		S3 = calloc(N_ANGLES, sizeof(double));
		S4 = calloc(N_ANGLES, sizeof(double));
		N  = calloc(N_ANGLES, sizeof(double));
		angles = calloc(N_ANGLES, sizeof(double));
		SlopeAngle = 0.;
		spot = spots->spot;
		for(i = 0; i < n; i++, spot++){
			idx = (*spot)->id;
			if(idx > N_RAY_MAX) continue;
			if(idx > N_ANGLE_MAX) idx -= N_ANGLES; // convert the ray number to the number of line
			x = (*spot)->c.x; y = (*spot)->c.y;
			S1[idx] += x;	// \sum x_i
			S2[idx] += y;	// \sum y_i
			S3[idx] += x*x;	// \sum (x_i)^2
			S4[idx] += x*y;	// \sum (x_i y_i)
			N[idx] += 1.;
		}
		for(i = 0; i < N_ANGLES; i++){
			double ang;
			// calculate an shift for i-th line
			angles[i] = atan2((S1[i]*S2[i] - N[i]*S4[i]) ,
				(S1[i]*S1[i] - N[i]*S3[i]));
			DBG("Angle: %g", angles[i]*180./M_PI);
			// check whether there's no jumps near +-pi
			if(i){
				ang = angles[i] - angles[i-1];
				DBG("Ang: %g", ang*180./M_PI);
				if(!Prefocal){
					if(ang > M_PI) angles[i] -= M_2PI;
					else if(ang > 0.) angles[i] -= M_PI;
					else if(ang < -M_2PI) angles[i] += M_2PI;
					else if(ang < -M_PI) angles[i] += M_PI;
				}else{
					if(ang < -M_PI) angles[i] += M_2PI;
					else if(ang < 0.) angles[i] += M_PI;
					else if(ang > M_2PI) angles[i] -= M_2PI;
					else if(ang > M_PI) angles[i] -= M_PI;
				}
			}
			ang = angles[i] - angles_0[i];
			if(ang > M_2PI) ang -= M_2PI;
			else if(ang < 0.) ang += M_2PI;
			SlopeAngle += ang;
			DBG("angle[%d] = %g; ang=%g; a0=%g\n", i, angles[i]*180./M_PI,
					ang*180./M_PI, angles_0[i]*180./M_PI);
		}
		SlopeAngle /= (double)N_ANGLES; // average
		red("Slope Angle: %g degr", SlopeAngle * 180. / M_PI);
		free(S1); free(S2); free(S3); free(S4); free(N); free(angles);
	}
	// sort spots by r to enumerate by rings
	{
		sort_spots_by_r(spots);
		spot = spots->spot;
		double r0 = (*spot)->r, r;
		// rtres == 1/7 of average distance between rings
		double rtres = (spot[n-1]->r - r0) / ((double)N_RINGS - 1.) / 7.;
		num = 0;
		green("\nSort by r\n");
		for(i = 0; i < n; i++, spot++){
			if((*spot)->id > N_RAY_MAX) continue; // we have already counts markers
			r = (*spot)->r;
			if(fabs(r-r0) > rtres){
				num+=100;
			}
			(*spot)->id += num;
		//	DBG("spot[%d]: r=%g, phi=%g, dr=%g\n", (*spot)->id, r,
		//		(*spot)->phi, r - r0);
			r0 = r;
		}
		if((num/100 + 1) > N_RINGS)
			bad_circumstances(TOO_MANY_RINGS);
	}
	/*
	 * At this stage it would be nice to display the picture boxes,
	 * numbered points, detected angle and center.
	 */

	/*
	 * Recalculate the the coordinates of centers points' for the translation
	 * to the coordinate system relative to the optical axis
	 * 			(rotation + reflection from the vertical axis).
	 * (\beta == SlopeAngle)
	 * x^0_i = -x_i\cos\beta - y_i\sin\beta;
	 * y^0_i = -x_i\sin\beta + y_i\cos\beta
	 */
	{
		double cosSA = cos(SlopeAngle);
		double sinSA = sin(SlopeAngle);
		double x,y;
		spot = spots->spot;
		if(!Prefocal){
			for(i = 0; i < n; i++, spot++){
				x = (*spot)->c.x; y = (*spot)->c.y;
				(*spot)->xC = -x*cosSA - y*sinSA;
				(*spot)->yC = -x*sinSA + y*cosSA;
			}
		}else{
			for(i = 0; i < n; i++, spot++){
				x = (*spot)->c.x; y = (*spot)->c.y;
				(*spot)->xC = x*cosSA + y*sinSA;
				(*spot)->yC = -x*sinSA + y*cosSA;
			}
		}
	}
ret:
	free(angles_pre); free(angles_post);
}

void spots_save(Spots *spots, gchar *filename){
	Spot **spot, *s;
	FILE *f;
	if(! spots || ! spots->spot || ! filename) return;
	int i, n = spots->n;
	sort_spots_by_id(spots);
	spot = spots->spot;
	f = fopen(filename, "w");
	if(!f){
		g_err(_("Can't open file"));
		return;
	}
	fprintf(f, "\t***** Bounding boxes *****\t** Coords w/o rot. comp. **\tCoords with rot. comp.\n");
	fprintf(f, "id\tx0\ty0\tw\th\txc\tyc\tw2\th2\txc\tyc\n");
	for(i = 0; i < n; i++, spot++){
		s = *spot;
		if(fprintf(f, "%03d\t%d\t%d\t%d\t%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n",
				s->id,
				s->box.x, s->box.y, s->box.w, s->box.h,
				s->c.x, s->c.y, s->c.rx, s->c.ry,
				s->xC, s->yC) <= 0){
			g_err(_("Can't write to file"));
			break;
		}
	}
	fclose(f);
}

#else // LEPTONICA_FOUND not set
Spots* spots_alloc(int n __attribute__((unused))){
	LEPTERR;
	return NULL;
}
Spot* spots_add(Spots *spots __attribute__((unused)),
			Spot *spot __attribute__((unused)), int copy __attribute__((unused))){
	LEPTERR;
	return NULL;
}
void spots_free(Spots **spots __attribute__((unused))){
	LEPTERR;
}
void sort_spots(struct Window *w __attribute__((unused))){
	LEPTERR;
}
#endif // LEPTONICA_FOUND
