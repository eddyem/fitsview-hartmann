//      fitsview.c - main project file
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
#include "fitsview.h"
#include "fits.h"
#include "gtk.h"
#include "opengl.h"
#include "open_dialog.h"
#include "filelist.h"

#define RED			"\033[1;31;40m"
#define GREEN		"\033[1;32;40m"
#define OLDCOLOR	"\033[0;0;0m"

int (*red)(const char *fmt, ...);
int (*green)(const char *fmt, ...);

int r_pr_(const char *fmt, ...){
	va_list ar; int i;
	printf(RED);
	va_start(ar, fmt);
	i = vprintf(fmt, ar);
	va_end(ar);
	printf(OLDCOLOR "\n");
	return i;
}

int g_pr_(const char *fmt, ...){
	va_list ar; int i;
	printf(GREEN);
	va_start(ar, fmt);
	i = vprintf(fmt, ar);
	va_end(ar);
	printf(OLDCOLOR "\n");
	return i;
}


Global_context *Global;

// init structures with zeros
void *init_struct(int size){
	int8_t *ptr = calloc(1, size);
	return (void*) ptr;
}
// copy structure
void *copy_struct(void *src, int size){
	int8_t *ptr = calloc(1, size);
	memcpy(ptr, src, size);
	return (void*)ptr;
}

// current time (for random number generator initialisation & for debug purposes)
double dtime(){
	double t;
	struct timeval tv;
	gettimeofday(&tv, NULL);
	t = tv.tv_sec + ((double)tv.tv_usec)/1e6;
	return t;
}

// if we'll save preferences global context should be read fom file
void init_contexts(){
	// limit spot size
	Global->maxSpotH = 100;
	Global->minSpotH = 5;
	Global->maxSpotW = 100;
	Global->minSpotW = 5;
}

// here will be a function to process command line arguments


int main(int argc, char *argv[]){
	gchar *filename = NULL, *str;
	if(isatty(STDOUT_FILENO)){ // make color output in tty
		red = r_pr_; green = g_pr_;
	}else{ // no colors in case of pipe
		red = printf; green = printf;
	}
#ifdef CUDA_FOUND
	DBG("CUDA found");
#else
	DBG("CUDA not found, using CPU instead");
#endif
	bindtextdomain(GETTEXT_PACKAGE, LOCALEDIR);
		bind_textdomain_codeset(GETTEXT_PACKAGE, "UTF-8");
	textdomain(GETTEXT_PACKAGE);
	INIT(Global, Global_context);
	init_contexts();
	init_main_window(&argc, &argv);
	getGLinfo(mainWindow->drawingArea);
	setlocale(LC_NUMERIC, "C"); // for right strto[df] working
	// Show all for begin
	mainWindow->context->visualMode = SHOW_ALL;
	mainWindow->image->data = NULL;
	if(argc > 1){
		str = basename(argv[1]);
		saved_path = get_current_dir_name();
		filename = g_strdup_printf("%s/%s", saved_path, str);
		fill_filelist(filename, mainWindow);
		change_image(filename, mainWindow);
		g_free(filename);
	}else
		gen_texture(mainWindow->image, mainWindow, FALSE);
	gtk_main();
	return 0;
}
