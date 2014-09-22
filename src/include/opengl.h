#ifndef _OPENGL_H_
#define _OPENGL_H_
#define GL_GLEXT_PROTOTYPES
#include "fitsview.h"
#include <GL/glu.h>
#include <GL/glext.h>
#include <GL/glut.h>
#include "gtk.h"

extern gboolean useVBO;

void initGl(GtkWidget *drawingArea);
void gen_texture(IMAGE *image, Window *window, gboolean redraw);
void force_redraw(GtkWidget *Area);
void pickRegion(double x, double y);
void conv_mouse_to_image_coords(double x, double y, double *X, double *Y, Window *window);
void conv_image_to_mouse_coords(double x, double y, double *X, double *Y, Window *window);
void init3D(Window *window);
void getGLinfo(GtkWidget *Area);
void freeGLmemory(Window *window);
#endif // _OPENGL_H_
