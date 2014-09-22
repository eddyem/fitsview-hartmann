#ifndef _GTK_H_
#define _GTK_H_
#include "fitsview.h"
#define EVENT_METHOD(i, x) GTK_WIDGET_GET_CLASS(i)->x

#define CH_WIN_TITLE(wndw, ...)	do{				\
	gchar *title = g_strdup_printf(__VA_ARGS__);	\
	gtk_window_set_title(GTK_WINDOW(wndw->window), title);\
	g_free(title);									\
	}while(0)

// Maximum number of cells lines the status window
#define SBAR_MAX 4
// Status bar fields
enum{
	 StatusText
	,StatusState
	,StatusCoords
	,StatusAdd
};

typedef struct{
	GLuint tex;		// texture itself
	GLuint w;		// its size
	GLuint h;
} TEXTURE;

typedef struct{
	double x;
	double y;
} XY;

// List of files in the current directory for this window
typedef struct{
	gint list_length;
	GList *files_list;
	GList *list_current;
	GList *list_end;
} FITSlist;

struct Spots;
// window description structure
struct Window{
	int id;						// identificator: MAIN_WINDOW, GRAPH_WINDOW, OPENGL_WINDOW
	GtkWidget *window;			// pointer to a window
	Context *context;			// context of a window
	struct Window *parent;		// parent window (NULL for the main)
	struct Window *graphWindow; // window's graph (while there's no - NULL)
	GtkWidget *drawingArea;		// drawing area (openGL)
	GtkWidget *hRule;			// rulers
	GtkWidget *vRule;
	GtkEntry *SBars[SBAR_MAX];	// an array of pointers to the cells of the status bar
	guint statusBlocks;			// status bar blocks amount
	IMAGE *image;				// image for a window
	IMAGE *image_transformed;	// and its transformation
	struct Spots *spots;		// spots array for a window
	TEXTURE *texture;			// texture or VBO buffers
	GtkRadioAction *LinLogMenu[2]; // scale lin/log
	double Zoom;				// zoom scale
	double Daspect;				// scale in image region (in window)
	XY move;					// image motions (in window)
	XY mouse;					// .x, .y - the beginning of the window coordinates in SC of a picture
	GLfloat Xangle;				// angles rotation relative to the X and Z axes
	GLfloat Zangle;
	FITSlist files;				// files in current directory for this window
};
typedef struct Window Window;
// window id
enum{
	MAIN_WINDOW,
	GRAPH_WINDOW,
	OPENGL_WINDOW,
	GL3D_WINDOW
};

extern Window *mainWindow;

void change_image(gchar *filename, Window *window);
void init_main_window(int *argc, char ***argv);
Window *init_window(Window *parent, int winId);
void destroy_window(Window* window);
void run_modal_window(GtkWindow *w, Window *parent);
gint run_modal_dialog(GtkDialog *dialog, Window *parent);
void g_err(gchar *text);
void set_Drulers(double x0, double y0, double xm, double ym, Window *window);
void set_Grulers(double y0, double xm, double ym, Window *window);
void show_histogram(Window *window);
void set_status_text(guint barName, gchar *text, Window *window);
void refresh_state(Window *window);
gchar *get_open_filename(Window *window);
void get_prefocal(Window *window, gboolean *prefocal);
#endif // _GTK_H_
