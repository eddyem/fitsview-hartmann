#ifndef _OPEN_DIALOG_H_
#define _OPEN_DIALOG_H_
#include "gtk.h"

extern gchar *saved_path;
GtkFileChooser *open_fits_dialog();
void gray2rgb(float gray, guchar *rgb);
#endif
