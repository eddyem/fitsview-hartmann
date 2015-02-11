//      contours.c - find isophotos
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

#include "contours.h"
#include "opengl.h"
#include "imtools.h"

// contours' minimum size limits
const int MIN_CONTOUR_SIZE = 4;
const int MIN_CLOSED_CONTOUR_SIZE = 7;
// maximum amount of contour levels
const int MAX_CONTOUR_LEVELS = 255;

// frees contour's data
void free_contour(Contour **c){
	cPoint *pp, *p;
	p = (*c)->first;
	while(p){ // go through all points of current contour
		pp = p; p = p->next;
		_FREE(pp);
	}
	_FREE(*c);
}

// delete contours' list for the picture
void free_contours(IMAGE *image){
	FNAME();
	int i;
	cList *l;
	Contour *cp, *c;
	if(image->Ncontours < 1 || !image->contours) return;
	for(i = 0; i < image->Ncontours; i++){ // go through all levels
		l = image->contours[i];
		c = l->first;
		while(c){ // go throug all contours of current level
			cp = c; c = c->next;
			free_contour(&cp);
		}
	}
}

// create new element of contours' list
cList *new_clist(int lvl){
	DBG("level: %d", lvl);
	cList *c = calloc(1, sizeof(cList));
	c->L = lvl;
	return c;
}
// create new contour
Contour *new_contour(){
	Contour *c = calloc(1, sizeof(Contour));
	return c;
}
// create new point
cPoint *new_point(){
	cPoint *p = calloc(1, sizeof(cPoint));
	return p;
}
// add contour c to list
void cList_add(cList *cl, Contour *c){
	cl->N++;
	//DBG("contour #%d added to list with levl %d", cl->N, cl->L);
	if(!cl->first)
		cl->first = c;
	else
		cl->last->next = c;
	cl->last = c;
}
// add point p to contour c
// headflag == 0 - add to tail; ==1 - add to head
void c_add(Contour *c, cPoint *p, char headflag){
	c->N++;
	if(!c->first || !c->last){
		c->first = p;
		c->last = p;
	}else{
		if(headflag){ // add point to head
			c->first->prev = p;
			p->next = c->first;
			c->first = p;
		}
		else{
			c->last->next = p; // add point to tail
			p->prev = c->last;
			c->last = p;
		}
	}
}

/*
 * copy contours from image in to image out
 * with shift (dX, dY): Xnew=Xold+dX, Ynew=Yold+dY
 */
int copy_contours(IMAGE *in, IMAGE *out, float dX, float dY){
	int i, N = in->Ncontours;
	cList *l;
	cPoint *pCur, *p;
	Contour *cCur, *c;
	if(N < 1 || !N) return TRUE;
	if(!(out->contours = calloc(in->Ncontours, sizeof(Contour)))) return FALSE;
	for(i = 0; i < N; i++){ // go through all levels
		l = in->contours[i];
		if(!(out->contours[i] = new_clist(i))) goto badExit;
		c = l->first;
		while(c){ // go throug all contours of current level
			if(!(cCur = new_contour())) goto badExit;
			pCur = c->first;
			while(pCur){
				if(!(p = new_point())) goto badExit;
				p->x = pCur->x+dX; p->y = pCur->y+dY;
				c_add(cCur, p, 0);
				pCur = pCur->next;
			}
			cList_add(out->contours[i], cCur);
			c = c->next;
		}
	}
	out->Ncontours = in->Ncontours;
	return TRUE;
badExit:
	free_contours(out); return FALSE;
}

// direction bits
enum{
	 D_RIGHT = 1
	,D_LEFT  = 2
	,D_DOWN  = 4
	,D_UP    = 8
};

// directions of next moving for mask values
// spetial values 6 & 9 also checks separately
int directions[16] = {
	0, // no isoline
	D_RIGHT | D_DOWN,
	D_LEFT  | D_DOWN,
	D_RIGHT | D_LEFT,
	D_RIGHT | D_UP,
	D_UP    | D_DOWN,
	0,  // special point 6
	D_LEFT  | D_UP,
	D_LEFT  | D_UP,
	0,  // special point 9
	D_UP    | D_DOWN,
	D_RIGHT | D_UP,
	D_RIGHT | D_LEFT,
	D_LEFT  | D_DOWN,
	D_RIGHT | D_DOWN,
	0 // no isoline
};

/*
 * coordinate shifts  {dx, dy} according to direction value bits
 * their values are choosen so, that even if there would be some non-zero
 * bits in direction mask, new direction of isoline search will be choosen
 * as possible nearest to "right" or " down" direction
 */
int dxdy[9][2] = {
	{0 , 0}, // 0
	{1 , 0}, // D_RIGHT
	{-1, 0}, // D_LEFT
	{0 , 0},
	{0 , 1}, // D_DOWN
	{0 , 0},
	{0 , 0},
	{0 , 0},
	{0 ,-1}  // D_UP
};
// new directions in simplex case, 2nd dimension - head flag (if==1, rotate CCW)
int newdirs[16][2] = {
	{0, 0},     // 0
	{D_RIGHT,D_RIGHT},// D_RIGHT
	{D_LEFT, D_LEFT}, // D_LEFT
	{D_RIGHT, D_LEFT},// D_RIGHT | D_LEFT
	{D_DOWN, D_DOWN}, // D_DOWN
	{D_RIGHT, D_DOWN},// D_DOWN | D_RIGHT
	{D_DOWN, D_LEFT}, // D_DOWN | D_LEFT
	{D_RIGHT, D_LEFT},// D_DOWN | D_RIGHT | D_LEFT
	{D_UP, D_UP},     // D_UP
	{D_RIGHT, D_UP},  // D_UP | D_RIGHT
	{D_LEFT, D_UP},   // D_UP | D_LEFT
	{D_RIGHT, D_UP},  // D_UP | D_RIGHT | D_LEFT
	{D_DOWN, D_UP},   // D_UP | D_DOWN
	{D_RIGHT, D_UP},  // D_UP | D_DOWN | D_RIGHT
	{D_DOWN, D_UP},   // D_UP | D_DOWN | D_LEFT
	{D_RIGHT, D_UP}   // D_UP | D_DOWN | D_RIGHT | D_LEFT
};

/*
 * Functions that calculate the position of the contour (linear interpolation)
 * P - the point at which are recorded the coordinates (relative to the center of the upper left quadrant)
 * to the resulting coordinates p shold be added the coordinates of the square to obtain the coordinates in image SC
 * A, b, c, d - intensity in the UL, UR DL and DR corners of the square
 * Lvl - contours level
 */
// right & up borders
void getRU(cPoint *p, float a, float b, float c __attribute__((unused)), float d, float I){
	float x = 0., y = 0.;
	if(fabs(b-a) > FLT_EPSILON) x = (I-a)/(b-a);
	if(fabs(d-b) > FLT_EPSILON) y = (I-b)/(d-b);
	p->x = x/2.+0.5; p->y = y/2.;
}
// left & up
void getLU(cPoint *p, float a, float b, float c, float d __attribute__((unused)), float I){
	float x = 0., y = 0.;
	if(fabs(b-a) > FLT_EPSILON) x = (I-a)/(b-a);
	if(fabs(c-a) > FLT_EPSILON) y = (I-a)/(c-a);
	p->x = x/2.; p->y = y/2.;
}
// left & down
void getLD(cPoint *p, float a, float b __attribute__((unused)), float c, float d, float I){
	float x = 0., y = 0.;
	if(fabs(d-c) > FLT_EPSILON) x = (I-c)/(d-c);
	if(fabs(c-a) > FLT_EPSILON) y = (I-a)/(c-a);
	p->x = x/2.; p->y = y/2.+0.5;
}
// right & down
void getRD(cPoint *p, float a __attribute__((unused)), float b, float c, float d, float I){
	float x = 0., y = 0.;
	if(fabs(d-c) > FLT_EPSILON) x = (I-c)/(d-c);
	if(fabs(d-b) > FLT_EPSILON) y = (I-b)/(d-b);
	p->x = x/2.+0.5; p->y = y/2.+0.5;
}
// up & down
void getUD(cPoint *p, float a, float b, float c, float d, float I){
	float x1 = 0., x2 = 1.;
	if(fabs(d-c) > FLT_EPSILON) x1 = (I-c)/(d-c);
	if(fabs(b-a) > FLT_EPSILON) x2 = (I-a)/(b-a);
	p->x = (x1+x2)/2.; p->y = 0.5;
}
// right & left
void getRL(cPoint *p, float a, float b, float c, float d, float I){
	float y1 = 0., y2 = 1.;
	if(fabs(d-b) > FLT_EPSILON) y1 = (I-b)/(d-b);
	if(fabs(c-a) > FLT_EPSILON) y2 = (I-a)/(c-a);
	p->x = 0.5; p->y = (y1+y2)/2.;
}
// bung for a case
void getNull(cPoint *p __attribute__((unused)), float a __attribute__((unused)),
	float b __attribute__((unused)), float c __attribute__((unused)),
	float d __attribute__((unused)), float I __attribute__((unused))){
		p->x = 0.; p->y = 0.;
}

// array of a funtions to coordinates calculation
typedef void (*FnXY)(cPoint *p, float a, float b, float c, float d, float I);
FnXY getXY[16] = {
	getNull, // 0
	getRD,
	getLD,
	getRL,
	getRU,
	getUD,
	getNull,//6 must be getLU or get RD
	getLU,
	getLU,
	getNull,//9 must be getLD or get RU
	getUD,
	getRU,
	getRL,//12
	getLD,
	getRD,
	getNull
};

static float *isoscale = NULL;
static int w, w1, h, h1;
static cList **contours = NULL;

/*
 * filling of the contour line, beginning at the point x1, y1
 * flag==1 - add points not to the tail of the circuit but to the head
 * cCur - contour to add the point
 * imdata - the original picture
 * lvl - isoline level number
 * olddir - direction (D_...) to a previous point
 * mask - pointer to current mask
 */
int proc_contour(int x1, int y1, int flag,
				Contour *cCur, float *imdata,
				int lvl, int olddir,
				unsigned char *mask){
	int frst = flag;
	float Lvl = isoscale[lvl]; // intensity at isoline level lvl
	float *point;
	int newdir;
	do{
		int mid = y1*w1 + x1, iid = y1*w + x1;
		unsigned char m = mask[mid];
		if(m == 0 || m > 14) break;
		if(!frst){
			cPoint *p = new_point();
			if(!p) return FALSE;
			c_add(cCur, p, flag);
			// get coordinates p->x, p->y (by simplified linear approximation)
			point = imdata + iid; // LU pixel for a square
			// find coordinates
			switch(m){ // check special points 6 & 9, their values depends on olddir
				case 6: // m = 7/14, mask=14/7
					switch(olddir){
						case D_LEFT: // next - UP or LEFT
						case D_UP:
							m = 7; mask[mid] = 14; // clear only used point
						break;
						default: // D_DOWN/D_RIGHT, next - RIGHT or UP
							m = 14; mask[mid] = 7;
					}
				break;
				case 9: // m = 13/11, mask=11/13
					switch(olddir){
						case D_LEFT: // next - DOWN or LEFT
						case D_DOWN:
							m = 13; mask[mid] = 11; // clear only used point
						break;
						default: // D_UP/D_RIGHT, next - RIGHT or DOWN
							m = 11; mask[mid] = 13;
					}
				break;
				default:
					mask[mid] = 0; // just clear used mask's point
			}
			getXY[m](p,point[0],point[1],point[w],point[w+1],Lvl);
			if(p->x < -0.5 || p->x > 1.5) p->x = 0.5;
			if(p->y < -0.5 || p->y > 1.5) p->y = 0.5;
			// add .5, beacause of mask's shift
			p->x += (float)x1 + .5;
			p->y += (float)y1 + .5;
		}
		gboolean found = FALSE; // whether next isoline point was found
		int x2, y2;
		do{ // find next contour's point
			unsigned char pp;
			// check right new direction
			newdir = newdirs[directions[m] & (~olddir)][flag];
			//newdir = directions[m] & (~olddir);
			if(newdir > 8){DBG("WTF? m=%d, newdir=%d, olddir=%d",m,newdir, olddir); break;} // This can't be so, but WTF if?
			x2 = x1 + dxdy[newdir][0];
			y2 = y1 + dxdy[newdir][1];
			if(x1 == x2 && y1 == y2) break;
			// is next square out of an picture?
			if(x2 < 0 || x2 >= w1 || y2 < 0 || y2 >= h1) break;
			int ii = y2*w1+x2;
			pp = mask[ii];
			// has the next square a contour points?
			if(pp == 0 || pp > 14) break;
			found = TRUE;
		}while(0);
		x1 = x2; y1 = y2;
		if(!found){
			break; // end of isoline
		}
		switch(newdir){ // calculate old direction
			case D_LEFT: olddir = D_RIGHT; break;
			case D_RIGHT: olddir = D_LEFT; break;
			case D_UP: olddir = D_DOWN; break;
			case D_DOWN: olddir = D_UP; break;
			default: olddir = 0;
		}
		frst = FALSE;
	}while(1); // cycle through the contour
	return TRUE;
}
/*
 * Processing of the point with coordinates x, y
 * imdata - the original picture
 * lvl - number of contour's level
 * mask - current mask
 * return FALSE if failed
 */
int process_it_(int x, int y, int lvl, float *imdata, unsigned char *mask){
	Contour *cCur;
	unsigned char pt0 = mask[y*w1+x];
	if(pt0 == 0 || pt0 > 14) return TRUE;;
	do{
		if(!contours[lvl]){ // countour wasn't created - create it
			contours[lvl] = new_clist(lvl);
			if(!contours[lvl]){
				g_err(_("Can't allocate memory"));
				return FALSE;
			}
		}
		cCur = new_contour();
		if(!cCur){
			g_err(_("Can't allocate memory"));
			return FALSE;
		}
		// start processing contour, adding points to its tail
		// D_LEFT - beacause we check points from left to the right
		if(!proc_contour(x, y, 0, cCur, imdata, lvl, D_LEFT, mask))
			return FALSE;
		// continue processing contour, adding points to its head
		// D_UP - beacause there's no more contour points upper
	//	if((y + 1) < w1) // mask(x,y) == 0, so we check down point
	//		if(!proc_contour(x, y+1, 1, cCur, imdata, lvl, D_UP, mask))
	//			return FALSE;
		if(cCur->N < MIN_CONTOUR_SIZE){ // contour is too short
			free_contour(&cCur);
			break;
		}
		if( (fabs(cCur->first->x - cCur->last->x) +
			fabs(cCur->first->y - cCur->last->y)) < 2.){ // closed contour
				cCur->closed = TRUE;
				if(cCur->N < MIN_CLOSED_CONTOUR_SIZE){ // contour is too short
					free_contour(&cCur);
					break;
				}
			}
		cList_add(contours[lvl], cCur);
	}while(0);
	return TRUE;
}


/*
 * Find isophotos on picture window->image
 * 		nLevl - amount of isophotos
 * 		type  - type of isophotos segmentation (the same as f->h for STEP filter)
 */
void find_contours(Window *window,
					int nLevl,
					int type){
	IMAGE *ima = window->image;
	int y, level;
	float *image = ima->data;
	if(nLevl < 2 || nLevl > MAX_CONTOUR_LEVELS) return;
	unsigned char *mask; // array for marching squares matrix
	Filter f;
	#ifdef EBUG
	double t0 = dtime();
	#endif
	w = ima->width; w1=w-1; h = ima->height; h1=h-1;
	mask = calloc((w1)*(h1),1); // mask - one mask pixel is square of 2x2 pixels on an picture
	if(!mask){g_err(_("Can't allocate memory"));return;}
	// prepare isoscale
	f.FilterType = STEP; f.w = nLevl; f.h = type;
	if(!fillIsoScale(&f, &isoscale, ima->stat.min, ima->stat.max - ima->stat.min)){
		g_err(_("Can't allocate memory"));goto endOF;}
	// Go through all levels
	contours = calloc(nLevl, sizeof(Contour)); // all-picture contours array
	if(!contours){g_err(_("Can't allocate memory"));goto endOF;}
	for(level = 0; level < nLevl; ++level){
		/*
		 *	1. level segmentation & build mask
		 * Mask by square (if zero point in left upper corner)
		 *  _______
		 * | a | b |
		 * |---+---|
		 * | c | d |
		 *  -------
		 * bits order: 0000abcd
		 */
		int y;
		float lvl = isoscale[level]; // current isoline's intensity level
		#pragma omp parallel for
		for(y = 0; y < h1; y++){
			int x;
			unsigned char *out = &mask[y*w1];
			float *in = &image[y*w];
			for(x = 0; x < w1; x++, in++, out++)
				*out = (unsigned char)
					((in[0]>lvl)<<3)   |
					((in[1]>lvl)<<2) |
					((in[w]>lvl)<<1) |
					(in[w+1]>lvl);
		}
		// go through all mask's points
		for(y = 0; y < h1; y++){
			int x;
			for(x = 0; x < w1; x++)
				if(!process_it_(x, y, level, image, mask)) goto endOF;
		}
	}
endOF:
	_FREE(mask);
	if(contours){
		for(y = 0; y < nLevl; y++){
			g_print("Level: %d", y);
			if(contours[y]) g_print(", %d contours\n", contours[y]->N);
			else g_print("\n");
		}
		// fill contour structure
		ima->contours = contours; contours = NULL;
		ima->Ncontours = nLevl;
		if(isoscale){ima->cLevels = isoscale; isoscale = NULL;}
	}
	DBG("time: %g", dtime()-t0);
	window->context->visualMode |= SHOW_ISOLINES;
	gen_texture(ima, window, TRUE);
	force_redraw(window->drawingArea);
}
