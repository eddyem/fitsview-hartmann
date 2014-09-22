//      terrain.c - functions for 3D view
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
#include "opengl.h"
#include "terrain.h"
#include "contours.h"
#include "fits.h"
#include "imtools.h" // for destroy_image


#define MINSZ	5		// minimal size of a selection region

/*
 * Calculate the normals for the point (x, y) of the image.
 * First, fill arrays of normals for the the lower right and upper
 * left triangles with one vertex at a given point.
 * Then for each point the total normal is calculated
 *
 * 				Principle
 *  *---B---*---*		(rd is the lower right triangle, ul is upper left)
 *  | / | / | / |	Normal at the point A at the map of triangles (as shown)
 *  E---A---C---*	is composed of normal Aul + Bdr + Cul + Adr + Dul + Edr
 *  | / | / | / | 	The normals for individual triangles (X - intensity at a point):
 *  *---D---*---*		(normal = vertical vector x horizontal vector)
 * 			Not normalized		Aul = (E-A, A-B, 1)
 * 								Adr = (A-C, D-A, 1)
 */

/*
 * Computing of the normals UL & DR
 * REWRITE FOR CUDA
 */
GLfloat *mkTrNormals(IMAGE *image){
	FNAME();
	int W = image->width, H = image->height, x, y, i;
	/*
	 * Normal for the the lower right and upper left triangles coming out of a given point
	 * Component 0 - UL, 1 - DR
	 * To simplify further calculations, supplemented by a frame in one node (zeros)
	 */
	GLfloat *tr_norms = calloc((W+2)*(H+2), 6*sizeof(GLfloat));
	int Wstride = (W + 2) * 6;
	int ww = W - 1, hh = H - 1;
	GLfloat *heights = image->data;
	GLfloat *ptr;
	#ifdef EBUG
		double t0 = dtime();
	#endif
	inline void soft_normalize(GLfloat *v){
		GLfloat l,x,y;
		x = v[0];
		y = v[1];
		l = sqrtf((x*x) + (y*y) + 1.f);
		v[0] = x / l;
		v[1] = y / l;
		v[2] = 1.f/ l;
	}
	inline void calc_UL(GLfloat *height, GLfloat *ptr){
		GLfloat A = *height;
		if(A){};
		ptr[0] = height[-1] - A; // E - A
		ptr[1] = height[-W] - A; // B - A
		ptr[2] = 1.f;
		soft_normalize(ptr);
	}
	inline void calc_DR(GLfloat *height, GLfloat *ptr){
		GLfloat A = *height;
		if(A){};
		ptr += 3; // move to DR
		ptr[0] = A - height[1]; // A - C
		ptr[1] = A - height[W]; // A - D
		ptr[2] = 1.f;
		soft_normalize(ptr);
	}
	inline void calc_two_norms(GLfloat *height, GLfloat *ptr){
		calc_UL(height, ptr);
		calc_DR(height, ptr);
	}
	ptr = &tr_norms[Wstride + 1];
	#pragma omp parallel for
	for(x = 0; x < ww; x++){
		calc_DR(heights+x, ptr + 6*x); // upper border without right point, only DR
	}
	i = 0;
	#pragma omp parallel for private(y, ptr, x, i)
	for(y = 1; y < hh; y++){
		ptr = &tr_norms[Wstride*(y+1)+1];
		i = y * W;
		calc_DR(heights+(i++), ptr); // left border - only DR
		ptr += 6;
		for(x = 1; x < ww; x++, ptr += 6){
			calc_two_norms(heights+(i++), ptr);
		}
		calc_UL(heights+i, ptr); // right border - only UL
	}
	ptr = &tr_norms[Wstride * H + 1];
	heights = &image->data[hh * W];
	#pragma omp parallel for
	for(x = 1; x < W; x++)
		calc_UL(heights, ptr+6*x); // down border - onlyUL
	/*
	ptr = &tr_norms[Wstride + 1];
	for(x = 0; x < ww; x++, ptr += 6, heights++)
		calc_DR(heights, ptr);
	heights++;
	for(y = 1; y < hh; y++){
		ptr = &tr_norms[Wstride*(y+1)+1];
		calc_DR(heights++, ptr);
		ptr += 6;
		for(x = 1; x < ww; x++, ptr += 6, heights++){
			calc_two_norms(heights, ptr);
		}
		calc_UL(heights++, ptr);
	}
	heights++;
	ptr = &tr_norms[Wstride * H + 2];
	for(x = 1; x < W; x++, heights++, ptr += 6)
		calc_UL(heights, ptr);
	*/
	DBG("time: %g", dtime() - t0);
	return tr_norms;
}

void computeNormals(Window *window){
	FNAME();
	int x, y;
	IMAGE *verts = window->image;
	IMAGE *norms = window->image_transformed;
	int W = verts->width, H = verts->height;
	int Wstride = (W + 2) * 6;
	#ifdef EBUG
		double t0 = dtime();
	#endif
	GLfloat *tr_norms = mkTrNormals(verts);
	GLfloat *ptr;
	GLfloat *normal = norms->data;
	inline void addVectors(GLfloat *norm, GLfloat *dbl){
		GLfloat *Aul, *Adr, *Bdr, *Cul, *Dul, *Edr;
		Aul = dbl; Adr = dbl + 3;
		Bdr = dbl - Wstride + 3;
		Cul = dbl + 6;
		Dul = dbl + Wstride;
		Edr = dbl - 3;
		*norm++ = (*Aul++) + (*Bdr++) + (*Cul++) + (*Adr++) + (*Dul++) + (*Edr++);
		*norm++ = (*Aul++) + (*Bdr++) + (*Cul++) + (*Adr++) + (*Dul++) + (*Edr++);
		*norm   = (*Aul) + (*Bdr) + (*Cul) + (*Adr) + (*Dul) + (*Edr);
	}
	inline void normalize(GLfloat *v){
		GLfloat l,x,y,z;
		x = v[0];
		y = v[1];
		z = v[2];
		l = sqrtf((x*x) + (y*y) + (z*z));
		v[0] = x / l;
		v[1] = y / l;
		v[2] = z / l;
	}
	/*
	int WW = W*3;
	#pragma omp parallel for private(y, ptr, normal)
	for(y = 0; y < H; y++){
		ptr = &tr_norms[Wstride*(y+1)+1];
		normal = &norms->data[y*WW];
		for(x = 0; x < W; x++, normal+=3, ptr+=6){
			addVectors(normal, ptr);
			normalize(normal);
		}
	}*/
	for(y = 0; y < H; y++){
		ptr = &tr_norms[Wstride*(y+1)+1];
		for(x = 0; x < W; x++, normal+=3, ptr+=6){
			addVectors(normal, ptr);
			normalize(normal);
		}
	}
	free(tr_norms);
	DBG("time: %g", dtime() - t0);
}

/*
 * A simple 3D-display selected in window->parent field
 * Image corners coordinates (x0, y0), (x, y)
 * Coordinates in CS IMAGES, ANGLES are outermost, not known
 * window - a window with the 3D image
 * X0, y0, x, y - the image area of the parent window, from which information will be taken
 * from_file - TRUE to load the file into the same window (then the coordinates do not care)
 * 				FALSE for download from the [sub] the image of the parent window
 */
gboolean loadTerrain(Window *window,
					double x0, double y0, double x, double y,
					gboolean from_file){
	FNAME();
	int i, j, imW;
	IMAGE *verts, *norms; // vertexes & normals (pointers to w->im, w->im_tr)
	IMAGE *image = NULL; // pointer to the parent image or a new image
	int X0,Y0, X1,Y1; // here will be stored the coordinates of the upper left and lower right corners
	int W, H; // Here are the dimensions of the selected area
	GLfloat *in, *out; // Pointers to the points in the image
	GLfloat scale; // scale (1/8 of the perimeter of selected region)
	#ifdef EBUG
		double t0 = dtime();
	#endif
	if(!from_file){
		if(x0 < x){ X0 = (int)(x0+0.5); X1 = (int)(x+0.5);}
		else{ X0 = (int)(x+0.5); X1 = (int)(x0+0.5);}
		if(y0 < y){ Y0 = (int)(y0+0.5); Y1 = (int)(y+0.5);}
		else{ Y0 = (int)(y+0.5); Y1 = (int)(y0+0.5);}
		W = X1 - X0; H = Y1 - Y0;
		DBG("subimage size: %dx%d", W, H);
		DBG("subimage @ (%d,%d)-(%d,%d)", X0, Y0, X1, Y1);
		if(W < MINSZ || H < MINSZ){
			g_err(_("Selected area too small"));
			return FALSE;
		}
		image = window->parent->image;
		if(!copy_contours(image, window->image,
			(float)-X0, (float)-Y0))
				g_err(_("Error occured when tried to add contours"));
		window->image->cLevels = image->cLevels;
		
		window->image->stat.min = image->stat.min;
		window->image->stat.max = image->stat.max; 
	}else{ // load to the same window from file
		gchar *filename;
		INIT(image, IMAGE);
		if( !(filename = get_open_filename(window)) ||
			!readfits(filename, image)){
				g_free(filename);
				destroy_image(image);
				return FALSE;
		}
		g_free(filename);
		X1 = W = image->width;
		Y1 = H = image->height;
		X0 = Y0 = 0;
	}
	// check if there are shoals with the filling of the structure of points and normals
	free(window->image->data); window->image->data = NULL;
	free(window->image_transformed->data); window->image_transformed->data = NULL;
	if(!image || !image->data){
		g_err(_("No image in parent window"));
		goto abnormal_exit;
	}
	imW = image->width;
	verts = window->image;
	norms = window->image_transformed;
	verts->data = calloc(W*H, sizeof(GLfloat));
	if(!verts->data){
		g_err(_("Can't allocate memory"));
		goto abnormal_exit;
	}
	// normals for each point
	norms->data = calloc(W*H*3, sizeof(GLfloat));
	if(!norms->data){
		g_err(_("Can't allocate memory"));
		free(verts->data);
		verts->data = NULL;
		goto abnormal_exit;
	}
	verts->width = W; verts->height = H;
	norms->width = W; norms->height = H;
	scale = ((GLfloat)W + (GLfloat)H) / 4.f / (image->stat.max-image->stat.min);
	// Copy the selected piece of image with scaling
	#pragma omp parallel for
	for(j = Y0; j < Y1; j++){
		//int l = (j - Y0)*W;
		in = &image->data[j * imW + X0];
		out = &verts->data[(j - Y0)*W];
		for(i = 0; i < W; i++)
			(*out++) = (*in++) * scale;
	}
	if(from_file)
		destroy_image(image);
	DBG("time: %g", dtime() - t0);
	return TRUE;
abnormal_exit:
	if(from_file)
		destroy_image(image);
	return FALSE;
}

gboolean createList(Window *window){
	FNAME();
	TEXTURE *texture = window->texture; // number of list stores in texture->tex
	GLfloat *Verts = window->image->data;
	GLfloat *Norms = window->image_transformed->data;
	int W = window->image->width, H = window->image->height;
	int x, y, stride = 3 * W;
	GLfloat dX = -(GLfloat)W / 2.f, dY = -(GLfloat)H / 2.f;
	GtkWidget *Area = window->drawingArea;
	if(!Area){
		g_err("???");
		return FALSE;
	}
	GdkGLContext *glContext = gtk_widget_get_gl_context(Area);
	GdkGLDrawable *glDrawable = gtk_widget_get_gl_drawable(Area);
	if(!gdk_gl_drawable_gl_begin(glDrawable, glContext))
		g_assert_not_reached();
	#ifdef EBUG
		double t0 = dtime();
	#endif
	if(!texture){
		INIT(window->texture, TEXTURE);
		texture = window->texture;
		if(!texture){
			g_err(_("Can't allocate memory for texture"));
			return FALSE;
		}
		if(useVBO){
			glGenBuffersARB(1, &texture->tex);
			glGenBuffersARB(1, &texture->w);
		}else{
			texture->tex = glGenLists(1);
			texture->w = W;
		}
		if(texture->tex == 0){
			g_err(_("Can't generate VBO / GL list"));
			return FALSE;
		}
	}
	texture->h = H;
	if(useVBO){
		GLint sz = W*H*3*sizeof(GLfloat), szi = W*(H-1)*2*sizeof(GLuint);
		glBindBufferARB(GL_ARRAY_BUFFER_ARB, texture->tex);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, texture->w);
		glBufferDataARB(GL_ARRAY_BUFFER_ARB, sz*2,
						0, GL_STATIC_DRAW_ARB);
		GLfloat *vertexes = malloc(sz);
		GLuint *indexes = malloc(szi);
		GLfloat *ptr = vertexes;
		#pragma omp parallel for
		for(y = 0; y < H; y++){
			ptr = &vertexes[y*W*3];
			Verts = &window->image->data[y*W];
			for(x = 0; x < W; x++, Verts++){
				*ptr++ = dX + x;
				*ptr++ = dY + y;
				*ptr++ = *Verts;
			}
		}
		H--;
		#pragma omp parallel for
		for(y = 0; y < H; y++){
			GLuint *iptr = &indexes[y*W*2];
			GLuint idx = y*W;
			for(int x = 0; x < W; x++, idx++){
				*iptr++ = idx + W;
				*iptr++ = idx;
			}
		}
		glBufferDataARB(GL_ELEMENT_ARRAY_BUFFER, szi, 0, GL_STATIC_DRAW_ARB);
		glBufferSubDataARB(GL_ELEMENT_ARRAY_BUFFER, 0, szi, indexes);
		glBufferSubDataARB(GL_ARRAY_BUFFER_ARB, 0, sz, vertexes);
		glBufferSubDataARB(GL_ARRAY_BUFFER_ARB, sz,
							sz, window->image_transformed->data);
		free(vertexes);
		free(indexes);
		free(window->image->data);
		window->image->data = NULL;
		free(window->image_transformed->data);
		window->image_transformed->data = NULL;
		//~ free(window->image_transformed); window->image_transformed = NULL;

	}else{
		glNewList(texture->tex, GL_COMPILE);
		H--;
		for(y = 0 ; y < H; y++){
			glBegin(GL_TRIANGLE_STRIP);
			for(x = 0; x < W; x++, Verts++, Norms+=3){
				glNormal3fv(Norms);
				glVertex3f(dX + x, dY + y, *Verts);
				glNormal3fv(&Norms[stride]);
				glVertex3f(dX + x, dY + (y + 1), Verts[W]);
			}
			glEnd();
		}
		glEndList();
	}
	gdk_gl_drawable_gl_end(glDrawable);
	DBG("time: %g", dtime() - t0);
	return TRUE;
}

/*
 * External function for the preparation of [[section] [parental]] picture
 * to download into 3D-box
 * window - 3D-window
 * x0, y0, x, y - the corners of the selected area from picture of the parent window
 * from_file - TRUE to load from the called file selection window
 * 				FALSE for download from the parent window
 * If from_file == TRUE, the coordinates of the corners of subimage are ignored
 */
void terrain_3D(Window *window,
				double x0, double y0, double x, double y,
				gboolean from_file){
	FNAME();
	if(!window) return;
	if(!loadTerrain(window, x0, y0, x, y, from_file)){
		if(!from_file) destroy_window(window);
		return;
	}
	computeNormals(window);
	if(!createList(window)){
		destroy_window(window);
		return;
	}
	init3D(window); // initialisation of openGL window
}
