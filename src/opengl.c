//      opengl.c - functions to work with openGL
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
#include "spots.h"
#include "imtools.h"
#include "contours.h"
#include "open_dialog.h"

#define BUFSIZE 512 // select buffer sizw
GLdouble ObjX = 0., ObjY = 0.; // mouse click coordinates (in window CS)
gboolean forcesel = FALSE;
gboolean useVBO = FALSE; // use VBO or lists

// these functions are absent in headers
void glBlendEquation(GLenum);
void glBlendFuncSeparate(GLenum, GLenum, GLenum, GLenum);

#define BADOGL()	do{g_err(_("OpenGL not supported"));	\
						g_assert_not_reached();}while(0)

gboolean queryExtension(char *extName){
	char *p = (char *)glGetString(GL_EXTENSIONS);
	char *end = p + strlen(p);
	while(p < end){
		size_t n = strcspn(p, " ");
		if((strlen(extName)==n) && (strncmp(extName,p,n)==0))
			return TRUE;
		p += (n + 1);
	}
	return FALSE;
}

void getGLinfo(GtkWidget *Area){
	GLint a,b,c,d;
	GdkGLDrawable *glDrawable = gtk_widget_get_gl_drawable(Area);
	GdkGLContext *glContext = gtk_widget_get_gl_context(Area);
	if (!gdk_gl_drawable_gl_begin(glDrawable, glContext))
		g_assert_not_reached();
	g_print(_("\nOpenGL info:\n"));
	glGetIntegerv(GL_RED_BITS, &a);
	glGetIntegerv(GL_GREEN_BITS, &b);
	glGetIntegerv(GL_BLUE_BITS, &c);
	glGetIntegerv(GL_ALPHA_BITS, &d);
	g_print("%s (R, G, B, alpha) = (%d, %d, %d, %d)\n", _("Color bits"),
			a, b, c, d);
	glGetIntegerv(GL_DEPTH_BITS, &a);
	glGetIntegerv(GL_STENCIL_BITS, &b);
	glGetIntegerv(GL_MAX_LIGHTS, &c);
	g_print("%s = %d, %s = %d, %s = %d\n", _("Depth bits"), a,
			_("Stencil bits"), b, _("Max amount of lights"), c);
	glGetIntegerv(GL_MAX_TEXTURE_SIZE, &a);
	glGetIntegerv(GL_MAX_CLIP_PLANES, &b);
	g_print("%s = %d, %s = %d\n", _("Max texture size"), a,
			_("Max clip planes"), b);
	glGetIntegerv(GL_MAX_MODELVIEW_STACK_DEPTH, &a);
	glGetIntegerv(GL_MAX_PROJECTION_STACK_DEPTH, &b);
	glGetIntegerv(GL_MAX_ATTRIB_STACK_DEPTH, &c);
	glGetIntegerv(GL_MAX_TEXTURE_STACK_DEPTH, &d);
	g_print("%s: %s = %d, %s = %d, %s = %d, %s = %d\n", _("Max stack depths"),
			_("modelview"), a, _("projection"), b, _("attrib"), c,
			_("texture"), d);
	//g_print("%s: %s\n", _("All extensions"), (char *)glGetString(GL_EXTENSIONS));
	useVBO = queryExtension("GL_ARB_vertex_buffer_object");
	if(useVBO) g_print(_("VBO is supported\n"));
	gdk_gl_drawable_gl_end(glDrawable);
}

void freeGLmemory(Window *window){
	if(!window->texture) return;
	GLuint *tex = &window->texture->tex;
	if(!*tex) return;
	if(useVBO)
		glDeleteBuffersARB(1, tex);
	else
		glDeleteLists(*tex, 1);
}

/*
 * Coordinates transformation from CS of drawingArea into CS of picture
 * 		x,y - pointer coordinates
 * 		X,Y - coordinates of appropriate point at picture
 */
void conv_mouse_to_image_coords(double x, double y,
								double *X, double *Y,
								Window *window){
	double a = window->Daspect / window->Zoom;
	*X = x * a - window->mouse.x;
	*Y = window->mouse.y - y * a;
}

void conv_image_to_mouse_coords(double X, double Y,
								double *x, double *y,
								Window *window){
	double a = window->Zoom / window->Daspect;
	*x = (X + window->mouse.x) * a;
	*y = (window->mouse.y - Y) * a;
}

gboolean configure(GtkWidget *Area, GdkEventExpose *event, Window *window){
	//FNAME();
	double a, A, W, H, w, h;
	double Zoom = window->Zoom;
	if(event)
		if(event->count) return FALSE;
	TEXTURE *texture = (TEXTURE*)window->texture;
	GdkGLDrawable *glDrawable = gtk_widget_get_gl_drawable(Area);
	GdkGLContext *glContext = gtk_widget_get_gl_context(Area);
	if (!gdk_gl_drawable_gl_begin(glDrawable, glContext))
		g_assert_not_reached();
	gdk_gl_drawable_wait_gdk(glDrawable);
	// set position
	double W0, H0, x0,y0, xm,ym;
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	W = W0 = Area->allocation.width;
	H = H0 = Area->allocation.height;
	A = W / H;
	glViewport(0, 0, W, H);
	if(texture){
		// compute the right dimensions for rectangle with picture
		w = texture->w; h = texture->h;
		a = w / h; w /= 2.; h /= 2.;
		if(A > a){
			W = h * A; H = h;
			//window->mouse.rx = (W0 - H0*a)/2.; window->mouse.ry = 0.;
			window->Daspect = h / H0 * 2.;
		}
		else{
			H = w / A; W = w;
			//window->mouse.ry = (H0 - W0/a)/2.; window->mouse.rx = 0.;
			window->Daspect = w / W0 * 2.;
		}
		// recalculate limits for the rulers
		x0 = W/Zoom - w + window->move.x / Zoom;
		window->mouse.x = x0;
		y0 = H/Zoom - h + window->move.y / Zoom;
		xm = 2. * W / Zoom - x0;
		ym = 2. * H / Zoom - y0;
		window->mouse.y = ym;
		//DBG("W=%g, w=%g, Zoom=%g, x0=%g, xm=%g, Da=%g", W, w, Zoom, x0, xm, Daspect);
		set_Drulers(x0, y0, xm, ym, window);
		glOrtho(-W,W, -H,H, -1., 1.);
	}
	glMatrixMode(GL_MODELVIEW);
	glBlendEquation(GL_FUNC_ADD);
	glEnable(GL_BLEND);
	gdk_gl_drawable_gl_end(glDrawable);
	refresh_state(window);
	return FALSE;
}

/*
 * Ellipse with center at (x,y) and radii (rx,ry)
 * plots as polyline with stride of 10 degr
 */
void draw_ellipce(GLfloat x, GLfloat y, GLfloat rx, GLfloat ry){
	GLfloat pi2 = M_PI * 2., step = pi2 / 36., angle;
	glBegin(GL_LINE_LOOP);
	for(angle = 0.; angle < pi2; angle += step)
		glVertex2f(x+rx*cos(angle), y+ry*sin(angle));
	glEnd();
}

/*
 * Graphical indication of recognized spots
 * window - the window with picture
 * w, h - half-width & half-height of picture (texture)
 * CX0, CY0 - zero point coordinates in picture CS (in center)
 * box == TRUE - draw region, == FALSE - set selection region
 */
void draw_squares(Window *window, double w, double h, double CX0, double CY0, gboolean box){
	int i, n;
	double x0, y0, x1, y1;
	BOX *bbox; Spot **spot; Coordinates *crds;
	Spots *spots = window->spots;
	if(!spots) return;
	n = spots->n;
	spot = spots->spot;
	glPushMatrix();
	for(i = 0; i < n; i++){
		bbox = &spot[i]->box;
		crds = &spot[i]->c;
		x0 = bbox->x - w; y0 = h - bbox->y;
		x1 = x0 + bbox->w; y1 = y0 - bbox->h;
		if(box){
			glColor3f(0., 1., 0.); // green is box borders
			glBegin(GL_LINE_LOOP);
				glVertex2f(x0, y0); glVertex2f(x1, y0);
				glVertex2f(x1, y1); glVertex2f(x0, y1);
			glEnd();
			if(window->context->transformMask & TRANS_HARTMANN){
				x0 = crds->x + CX0;
				y0 = crds->y + CY0;
				glColor3f(1., 0., 0.); // red is spot half-widths
				draw_ellipce(x0, y0, crds->rx, crds->ry);
				//glColor3f(0., 1., 0.);
			}
		}else{ // box == FALSE - set frames for selection
			glLoadName(i);
			glRectf(x0, y0, x1, y1);
		}
	}
	if(box){ // draw cross - CS zero
		glColor3f(1., 0., 1.);
		glBegin(GL_LINES);
			glVertex2f(CX0-15., CY0-15.); glVertex2f(CX0+15., CY0+15.);
			glVertex2f(CX0-15., CY0+15.); glVertex2f(CX0+15., CY0-15.);
		glEnd();
	}
	glPopMatrix();
}

/*
 * graphical indication of isolines
 * window - the window
 * w, h - half-width & half-height of picture (texture)
 * CX0, CY0 - zero point coordinates in picture CS (in center)
 * _2D == TRUE for 2D images
 */
void draw_isolines(Window *window, double w, double h, gboolean _2D){
	//FNAME();
	IMAGE *ima = window->image;
	if(!ima) return;
	Contour *c;
	cPoint *p;
	cList **cl = window->image->contours;
	int i, n = window->image->Ncontours;
	if(!cl || n < 1) return;
	float m = ima->stat.min, wd = ima->stat.max - m;
	float *lvls = ima->cLevels, curLevel;
	guchar colr[3] = {255,0,0}; // isolines' colors
	glPushMatrix();
	for(i = 0; i < n; i++, cl++){ // cycle by isolines
		if(!*cl) continue; // given level doesn't have isolines
		c = (*cl)->first; // first contour for given level
		if(!c || (*cl)->N < 1) continue;
		gray2rgb((lvls[(*cl)->L]-m)/wd, colr); // convert level into color
		glColor3ubv(colr);
		if(_2D) curLevel = 0.;
		else curLevel = lvls[(*cl)->L] * ((GLfloat)w + (GLfloat)h) / 2./ wd;
		//DBG("curLevel=%g, gray=%g, color=(%d,%d,%d)", curLevel, (lvls[(*cl)->L]-m)/wd, colr[0],colr[1],colr[2]);
		while(c){ // cycle by contours
			p = c->first;
			if(!p || c->N < 1){ c = c->next; continue;}
			if(c->closed)
				glBegin(GL_LINE_LOOP);
			else
				glBegin(GL_LINE_STRIP);
			while(p){ // cycle by contour points
				glVertex3f(p->x-w, p->y-h, curLevel);
				p = p->next;
			}
			glEnd();
			c = c->next;
		}
	}
	glPopMatrix();
}

void processHits(Window *window, GLint hits, GLuint buffer[]){
	unsigned int i, j;
	GLuint names, *ptr;
	Spot *spot;
	//Coordinates crds;
	gchar buf[256];
	int l = 255;
	Spots *spots = window->spots;
	if(!spots) return;
	printf ("hits = %d\n", hits);
	ptr = (GLuint *) buffer;
	if(hits < 1){
		set_status_text(StatusText, _("Empty region"), window);
		return;
	}
	for (i = 0; i < (unsigned int)hits; i++){
		names = *ptr;
		DBG("number of names for this hit = %d\n", names);
		ptr += 3;
		g_print(_("Selected spot[s]:\n"));
		if(names == 1)
			g_sprintf(buf, _("Selected spot: "));
		else
			g_sprintf(buf, _("Selected spots: "));
		l -= strlen(buf);
		for (j = 0; j < names; j++){
			if(*ptr < (GLuint)spots->n){
				spot = spots->spot[*ptr];
				g_print("id: %d, xC=%.2f, yC=%.2f, r=%.2f, phi=%.2f, ",
						spot->id, spot->xC, spot->yC, spot->r, spot->phi);
				g_print("cn: (%.1f, %.1f), dr:(%.1f, %.1f)\n",
					spot->c.x, spot->c.y, spot->c.rx, spot->c.ry);
				g_print("box X: [%d, %d], box Y: [%d, %d], cn: [%d, %d]\n",
					spot->box.x, spot->box.x+spot->box.w,
					spot->box.y, spot->box.y+spot->box.h,
					spot->box.x+spot->box.w/2, spot->box.y+spot->box.h/2);
			/*	get_gaussian_center(&spot->box, &crds);
				g_print("\tCenter: (%.1f, %.1f), wd: (%.1f, %.1f)",
						crds.x, (double)image->height-crds.y-1., crds.rx, crds.ry);
			*/
				l -= g_snprintf(&buf[strlen(buf)], l, "%d ", spot->id);
			}
			ptr++;
		}
		g_print("\n");
		set_status_text(StatusText, buf, window);
	}
}

gboolean expose(GtkWidget *Area, GdkEventExpose *event, Window *window){
	TEXTURE *texture = (TEXTURE*)window->texture;
	double w, h;
	GLfloat tred[4] = {1.,0.,0.,1.};
	GLfloat twhite[4] = {1.,1.,1.,1.};
	if(event)
		if(event->count) return FALSE;
	//FNAME();
	GdkGLDrawable *glDrawable = gtk_widget_get_gl_drawable(Area);
	GdkGLContext *glContext = gtk_widget_get_gl_context(Area);
	if (!gdk_gl_drawable_gl_begin(glDrawable, glContext))
		g_assert_not_reached();
	glClearColor(.1, .1, .1, 1.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glColor3fv(tred);
	if(texture){
		w = texture->w / 2.; h = texture->h / 2.;
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
		glTranslatef(window->move.x, window->move.y, 0.);
		glScalef(window->Zoom, window->Zoom, 1.);
		//DBG("moveXY=%g, %g; Zoom=%g, texture: id=%d, w=%d",
		//	window->move.x, window->move.y, window->Zoom, texture->tex, texture->w);
		glEnable(GL_TEXTURE_2D);
		glBindTexture(GL_TEXTURE_2D, texture->tex);
		//DBG("texture: %d", texture->tex);
		//DBG("w=%g, h=%g\n", w, h);
		glTexEnvfv(GL_TEXTURE_ENV, GL_TEXTURE_ENV_COLOR, twhite);
		glDisable(GL_BLEND);
		glBegin(GL_QUADS);
			glTexCoord2f(0.0f, 0.0f); glVertex2f(-w, -h );
			glTexCoord2f(1.0f, 0.0f); glVertex2f( w, -h );
			glTexCoord2f(1.0f, 1.0f); glVertex2f( w,  h );
			glTexCoord2f(0.0f, 1.0f); glVertex2f(-w,  h );
		glEnd();
		glDisable(GL_TEXTURE_2D);
		/*glBegin(GL_LINE_LOOP);
			glTexCoord2f(0.0f, 0.0f); glVertex2f(-w, -h );
			glTexCoord2f(1.0f, 0.0f); glVertex2f( w, -h );
			glTexCoord2f(1.0f, 1.0f); glVertex2f( w,  h );
			glTexCoord2f(0.0f, 1.0f); glVertex2f(-w,  h );
		glEnd();*/
		glEnable(GL_BLEND);

		// draw addition information
		if(window->context->visualMode &
				(SHOW_POTBOXES|SHOW_POTSELECT|SHOW_ISOLINES)){
			GLuint selectBuf[BUFSIZE];
			GLint hits;
			double CX0 = AxisX - w, CY0 = h - AxisY;
			gboolean showiso = window->context->visualMode & SHOW_ISOLINES;
			gboolean showbox = window->context->visualMode & SHOW_POTBOXES;
			gboolean showsel = window->context->visualMode & SHOW_POTSELECT;
			GLint viewport[4];
			if(showsel && forcesel){ // get coordinates from selection buffer
				glPushMatrix();
				glGetIntegerv(GL_VIEWPORT, viewport);
				glSelectBuffer(BUFSIZE, selectBuf);
				glRenderMode(GL_SELECT);
				glInitNames();
				glPushName(0);
				glMatrixMode(GL_PROJECTION);
				glPushMatrix();
				glLoadIdentity();
				/*  create  pixel picking region near cursor location      */
				DBG("x=%g, y=%g", ObjX, ObjY);
				gluPickMatrix(ObjX, (viewport[3] - ObjY), 5.0, 5.0, viewport);
				//gluPickMatrix(ObjX, ObjY, 10.0, 10.0, viewport);
				double W = Area->allocation.width;
				double H = Area->allocation.height;
				double A = W / H;
				if(A > w / h){
					W = h * A; H = h;}
				else{
					H = w / A; W = w;}
				glOrtho(-W,W, -H,H, -1., 1.);
				glMatrixMode(GL_MODELVIEW);
				if(showbox) draw_squares(window, w, h, CX0, CY0, FALSE);
				if(showiso) draw_isolines(window, w, h, TRUE);
				glMatrixMode(GL_PROJECTION);
				glPopMatrix();
				hits = glRenderMode(GL_RENDER);
				processHits(window, hits, selectBuf);
				forcesel = FALSE;
				glPopMatrix();
			}
			glMatrixMode(GL_MODELVIEW);
			if(showbox) draw_squares(window, w, h, CX0, CY0, TRUE);
			if(showiso) draw_isolines(window, w, h, TRUE);
		}
/*
		// layers mix ===>
		glLoadIdentity();
		//glTranslatef(moveX,moveY,0);
		glScalef(window->Zoom, window->Zoom, 1.);
		glTexEnvfv(GL_TEXTURE_ENV, GL_TEXTURE_ENV_COLOR, tred);
		glBlendFunc(GL_ONE, GL_ONE);
		glEnable(GL_TEXTURE_2D);
			glBegin(GL_QUADS);
				glTexCoord2f(0.0f, 0.0f); glVertex2f(-w, -h );
				glTexCoord2f(1.0f, 0.0f); glVertex2f( w, -h );
				glTexCoord2f(1.0f, 1.0f); glVertex2f( w,  h );
				glTexCoord2f(0.0f, 1.0f); glVertex2f(-w,  h );
			glEnd();
		glDisable(GL_TEXTURE_2D);
		// <====
*/
	}
	else DBG("No texture");
	gdk_gl_drawable_swap_buffers(glDrawable);
	gdk_gl_drawable_wait_gl(glDrawable);
	gdk_gl_drawable_gl_end(glDrawable);
	return FALSE;
}

/*
gboolean idle(gpointer userData){
	GtkWidget *Area = GTK_WIDGET(userData);
	gdk_window_invalidate_rect(Area->window, &Area->allocation, FALSE);
	return TRUE;
}*/

void initGl(GtkWidget *Area){
	FNAME();
	GdkGLConfig *glConfig;
	GdkGLConfigMode mode =	GDK_GL_MODE_RGB   |
							GDK_GL_MODE_DEPTH |
							GDK_GL_MODE_ALPHA |
							GDK_GL_MODE_DOUBLE;
	if(!gdk_gl_query_extension())
        BADOGL();
	glConfig = gdk_gl_config_new_by_mode(mode);
	if(!glConfig)
		glConfig = gdk_gl_config_new_by_mode(mode & ~GDK_GL_MODE_ALPHA);
	if(!glConfig)
		BADOGL();
	if(!gtk_widget_set_gl_capability(Area, glConfig, NULL, TRUE,
									GDK_GL_RGBA_TYPE))
		BADOGL();
}

/*
 * generation of textures to show in Area widget of image
 * if redraw == TRUE texture substitutes previous otherwise generate the new one
 */
void gen_texture(IMAGE *image, Window *window, gboolean redraw){
	//~ FNAME();
	TEXTURE *texture = window->texture;
	GtkWidget *Area = window->drawingArea;
	GdkGLContext *glContext = gtk_widget_get_gl_context(Area);
	GdkGLDrawable *glDrawable = gtk_widget_get_gl_drawable(Area);
	if (!gdk_gl_drawable_gl_begin(glDrawable, glContext))
		g_assert_not_reached();
	glEnable(GL_TEXTURE_2D);
	if(!texture){
		INIT(window->texture, TEXTURE);
		texture = window->texture;
		if(!texture){
			g_err(_("Can't allocate memory for texture"));
			return;
		}
		//~ g_print(_("Generating texture\n"));
		glGenTextures(1, &texture->tex);
		//~ DBG("%s, tex=%d\n", gluErrorString(glGetError()), texture->tex);
	}
	//~ DBG("Image size: w=%d, h=%d", image->width, image->height);
	glBindTexture(GL_TEXTURE_2D, texture->tex);
	if(image && image->data){
		int i,j, k, w = image->width, h = image->height;
		GLfloat *tex = NULL;
		GLfloat min = image->stat.min, wd = image->stat.max - min;
		GLfloat *ptri = image->data, *ptro;
		DBG("min = %g, wd = %g", min, wd);
		if(image->stat.max > 1. || min < 0.){
			tex = malloc(w * h * sizeof(GLfloat));
			ptro = tex;
			#ifdef EBUG
				double t0 = dtime();
			#endif
			#pragma omp parallel for
			for(i=0; i<h; i++){
				k = i*w;
				for(j=0; j<w; j++, k++){//, ptri++, ptro++){
					//*ptro = (*ptri - min) / wd;
					ptro[k] = (ptri[k] - min) / wd;
				}
			}
			DBG("time: %g", dtime() - t0);
		}else ptro = image->data;
		if(redraw)
			glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, w, h,
					GL_LUMINANCE, GL_FLOAT, ptro); // s/image->data/tex/
		else
			glTexImage2D(GL_TEXTURE_2D, 0, GL_INTENSITY, w, h,
					0, GL_LUMINANCE, GL_FLOAT, ptro);
		//DBG("%s\n", gluErrorString(glGetError()));
		//DBG("W: %d, H:%d, im[100]=%f\n", image->width, image->height, image->data[100]);
		texture->w = image->width; texture->h = image->height;
		free(tex);
	}
	else{
		DBG("Empty texture");
		glTexImage2D(GL_TEXTURE_2D, 0, GL_INTENSITY, 10., 10.,
					0, GL_INTENSITY, GL_FLOAT, NULL);
	}
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_COMBINE);
	glTexEnvf(GL_TEXTURE_ENV, GL_SOURCE0_RGB, GL_TEXTURE);
	glTexEnvf(GL_TEXTURE_ENV, GL_SOURCE1_RGB, GL_CONSTANT);
	//glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_NEAREST);
	//glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_NEAREST);
	glDisable(GL_TEXTURE_2D);
	gdk_gl_drawable_gl_end(glDrawable);
/*
 * Events:
 * 		expose - redraw part of an image
 * 		configure - redraw window when size/position changed
 * if handler returns TRUE signal processing stops
 * otherwise it can be process by other handlers
 */
	g_signal_connect(Area, "expose-event",
					G_CALLBACK(expose), window);
	g_signal_connect(Area, "configure-event",
					G_CALLBACK(configure), window);
	//~ DBG("end");
	force_redraw(Area);
}

void force_redraw(GtkWidget *Area){
	gtk_widget_queue_resize(Area);
	gtk_widget_queue_draw(Area);
}

void pickRegion(double x, double y){
	FNAME();
	forcesel = TRUE;
	ObjX = x;
	ObjY = y;
}

gboolean configure3D(GtkWidget *Area, GdkEventExpose *event, Window *window){
	FNAME();
	GdkGLContext *glContext = gtk_widget_get_gl_context(Area);
	GdkGLDrawable *glDrawable = gtk_widget_get_gl_drawable(Area);
	double w, h;
	GLfloat ratio;
if(!window)return FALSE;
	//~ GLfloat H = (GLfloat)window->texture->h;
	if(event)
		if(event->count) return FALSE;
	if (!gdk_gl_drawable_gl_begin(glDrawable, glContext))
		g_assert_not_reached();
	w = Area->allocation.width;
	h = Area->allocation.height;
	ratio = w / h;
	glViewport(0, 0, w, h);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(45., ratio, 0.1, 1e10);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gdk_gl_drawable_gl_end(glDrawable);
	return FALSE;
}

gboolean expose3D(GtkWidget *Area, GdkEventExpose *event, Window *window){
	FNAME();
	TEXTURE *texture = window->texture;
	GdkGLContext *glContext = gtk_widget_get_gl_context(Area);
	GdkGLDrawable *glDrawable = gtk_widget_get_gl_drawable(Area);
	GLfloat H = (GLfloat)window->image->height;
	if(event)
		if(event->count) return FALSE;
	if (!gdk_gl_drawable_gl_begin(glDrawable, glContext))
		g_assert_not_reached();
	GLfloat spotDir[] = {0.,0.,-1.,1.};
	GLfloat pos[] = {0.,0.,0.,1.};
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glLoadIdentity();
	glLightfv(GL_LIGHT1,GL_POSITION,pos);
	glLightfv(GL_LIGHT1,GL_SPOT_DIRECTION, spotDir);
	if(window->context->drawingMode == GAME_MODE){ // mouse navigation
		GLfloat lx, ly, lz, za = window->Zangle, xa = window->Xangle;
		GLfloat x = window->move.x, y = window->move.y, z = H * window->Zoom;
		lx = cosf(za)*sinf(xa);
		ly = -sinf(za);
		lz = -cosf(za)*cosf(xa);
		gluLookAt(x, y, z,
		  x+lx, y+ly, z+lz,
		  0., 1., 0.);
			//~ GLfloat w = Area->allocation.width, h = Area->allocation.height;
		//GLfloat z = H * window->Zoom;
			//~ GLfloat yangle = 180. / M_PI * atan2f(window->Xangle - w/2., z);
			//~ GLfloat xangle = 180. / M_PI * atan2f(window->Zangle - h/2., z);
			//~ glRotatef(yangle, 0.,1.,0.);
			//~ glRotatef(xangle, 1.,0.,0.);
			//~ glTranslatef(window->move.x, 0., -z);
		//glRotatef(window->Xangle, 0.,1.,0.);
			//~ glRotatef(window->Zangle, 1.,0.,0.);
		//glRotatef(window->Zangle, cosf(window->Xangle*M_PI/180.),0.,sinf(window->Xangle*M_PI/180.));
		//glTranslatef(window->move.x, 0., -z);
	}else{ // key navigation
		glTranslatef(window->move.x, window->move.y, -H * window->Zoom);
		glRotatef(window->Xangle, 1.,0.,0.);
		glRotatef(window->Zangle, 0.,0.,1.);
	}
	glPushMatrix();
	glEnable(GL_LIGHTING);
	if(useVBO){
		GLint y, w = window->image->width, h = window->image->height;
		GLint pts = w*2;
		size_t sz = w * h * 3 * sizeof(GLfloat);
		glBindBufferARB(GL_ARRAY_BUFFER_ARB, texture->tex);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, texture->w);
		glEnableClientState(GL_NORMAL_ARRAY);
		glEnableClientState(GL_VERTEX_ARRAY);
		glEnableClientState(GL_INDEX_ARRAY);
		glNormalPointer(GL_FLOAT, 0, (GLfloat*)sz);
		glIndexPointer(GL_UNSIGNED_INT, 0, 0 );
		glVertexPointer(3, GL_FLOAT, 0, 0);
		h--;
		sz = 0; w = pts * sizeof(GLfloat);
		for(y = 0; y < h; y++, sz += w)
			glDrawElements(GL_TRIANGLE_STRIP, pts, GL_UNSIGNED_INT, (GLvoid*)(sz));
		glDisableClientState(GL_VERTEX_ARRAY);
		glDisableClientState(GL_NORMAL_ARRAY);
		glDisableClientState(GL_INDEX_ARRAY);
		glBindBufferARB(GL_ARRAY_BUFFER_ARB, 0);
	}else
		glCallList(texture->tex);
	glPushMatrix();
	glDisable(GL_LIGHTING);
	glLineWidth(3.0);
	draw_isolines(window, ((double)window->image->width)/2.+.5,
		((double)window->image->height)/2.+.5, FALSE);
	/*glColor3f(0.,0.,1.);
	glBegin(GL_QUADS);
		glVertex3f(-100.f, -100.f,  -12.f );
		glVertex3f(-100.f,  100.f,  -12.f );
		glVertex3f( 100.f,  100.f,  -12.f );
		glVertex3f( 100.f, -100.f,  -12.f );
	glEnd();
	glBegin(GL_LINES);
		glVertex3fv(pos);
		glColor3f(1.f,0.f,0.f);
		glVertex3fv(spotDir);
	glEnd();*/
	glPopMatrix();
	gdk_gl_drawable_swap_buffers(glDrawable);
	gdk_gl_drawable_wait_gl(glDrawable);
	gdk_gl_drawable_gl_end(glDrawable);
	return FALSE;
}

void init3D(Window *window){
	FNAME();
	GtkWidget *Area = window->drawingArea;
	GdkGLContext *glContext = gtk_widget_get_gl_context(Area);
	GdkGLDrawable *glDrawable = gtk_widget_get_gl_drawable(Area);
	GLfloat lAmbient[] = {0.3,0.3,0.3,1.0};
	GLfloat lDiffuse[] = {1.0,1.0,1.0,1.0};
	GLfloat lSpecular[]= {1.0,1.0,1.0,1.0};
	GLfloat l2Ambient[] = {0.0,0.0,0.1,1.0};
	GLfloat l2Diffuse[] = {0.0,0.0,0.3,1.0};
	GLfloat l2Specular[]= {0.0,0.0,0.5,1.0};
	GLfloat mShininess[] = {110.0};
	GLfloat mSpecular[] = {0.8,0.8,1.0,1.0};
	GLfloat mDiffuse[] = {0.8,0.8,1.0,1.0};
	GLfloat mAmbient[] = {0.2,0.2,0.2,1.0};
	if (!gdk_gl_drawable_gl_begin(glDrawable, glContext))
		g_assert_not_reached();
	glEnable(GL_RESCALE_NORMAL);
	glEnable(GL_DEPTH_TEST);
	//~ glEnable(GL_CULL_FACE);
	//~ glCullFace(GL_BACK);
	GLfloat spotDir[] = {0.,0.,-1.,1.};
	GLfloat pos[] = {0.,0.,0.,1.};
	glLightfv(GL_LIGHT1,GL_POSITION,pos);
	glLightfv(GL_LIGHT1,GL_SPOT_DIRECTION, spotDir);
	glLightf(GL_LIGHT1,GL_SPOT_CUTOFF, 90.);
	glLightf(GL_LIGHT1,GL_SPOT_EXPONENT, 20.);
	glLightfv(GL_LIGHT0,GL_AMBIENT,l2Ambient);
	glLightfv(GL_LIGHT0,GL_DIFFUSE,l2Diffuse);
	glLightfv(GL_LIGHT0,GL_SPECULAR,l2Specular);
	glLightfv(GL_LIGHT1,GL_AMBIENT,lAmbient);
	glLightfv(GL_LIGHT1,GL_DIFFUSE,lDiffuse);
	glLightfv(GL_LIGHT1,GL_SPECULAR,lSpecular);
	glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER,GL_TRUE);
	glMaterialfv(GL_FRONT, GL_SPECULAR, mSpecular);
	glMaterialfv(GL_FRONT, GL_SHININESS,mShininess);
	glMaterialfv(GL_FRONT, GL_AMBIENT, mAmbient);
	glMaterialfv(GL_FRONT, GL_DIFFUSE, mDiffuse);
	//~ glMaterialfv(GL_BACK, GL_SPECULAR, mSpecular);
	//~ glMaterialfv(GL_BACK, GL_SHININESS,mShininess);
	//~ glMaterialfv(GL_BACK, GL_AMBIENT, mGray);
	//~ glMaterialfv(GL_BACK, GL_DIFFUSE, mWhite);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHT1);
	gdk_gl_drawable_gl_end(glDrawable);
	g_signal_connect(Area, "expose-event",
					G_CALLBACK(expose3D), window);
	g_signal_connect(Area, "configure-event",
					G_CALLBACK(configure3D), window);
	force_redraw(Area);
}
