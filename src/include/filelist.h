#ifndef __FILELIST_H__
#define __FILELIST_H__

#include "fitsview.h"
#include "gtk.h"

enum{
	CURRENT,
	FIRST,
	LAST,
	PREVIOUS,
	NEXT
};

void free_filelist(Window *window);
gint fill_filelist(gchar *filename, Window *window);
gchar *get_filename(guchar type, Window *window);
gboolean get_ext(gchar *filename);
#endif // __FILELIST_H__
