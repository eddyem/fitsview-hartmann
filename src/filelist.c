// filelist.c - functions to work with list of FITS files in current directory
//
// Copyright 2011 Edward V. Emelianoff <eddy@sao.ru>
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
// MA 02110-1301, USA.

#include "filelist.h"
#include "fits.h"
/*
 * receiving the file name of list, type can be one of the following:
 *		CURRENT		- current file
 *		FIRST		- first file
 *		LAST		- last file
 *		PREVIOUS	- previous file
 *		NEXT		- next file
 */
gchar *get_filename(guchar type, Window *window){
	DBG("length of list = %d", window->files.list_length);
	if(window->files.list_length < 1) return NULL;
	if(!window->files.list_current){
		DBG("No current file. WTF?");
		return NULL;
	}
	//~ g_list_foreach(window->files.files_list, (GFunc)prnt, NULL);
	switch(type){
		case FIRST:
			if(window->files.list_current == window->files.files_list){
				DBG("1st or only file");
				return NULL; // try to open current opened file or there's no files
			}
			window->files.list_current = window->files.files_list;
		break;
		case LAST:
			if(window->files.list_length < 2 ||
				window->files.list_current == window->files.list_end){
				DBG("last or only file");
				return NULL;
			}
			window->files.list_current = window->files.list_end;
		break;
		case PREVIOUS:
			if(!window->files.list_current->prev){
				DBG("1st file");
				return NULL;
			}
			window->files.list_current = window->files.list_current->prev;
		break;
		case NEXT:
			if(!window->files.list_current->next){
				DBG("last file");
				return NULL;
			}
			window->files.list_current = window->files.list_current->next;
		break;
		default: // current
		break;
	}
	if(!window->files.list_current || !window->files.list_current->data){
		DBG("no current file? WTF?");
		return NULL;
	}
	DBG("change filename: %s\n\n", (gchar*)window->files.list_current->data);
	return (gchar*)window->files.list_current->data;
}

// get filename suffix
gboolean get_ext(gchar *filename){
	gchar *basename, *ext;
	basename = strrchr(filename, '/');
	if(!basename) basename = filename;
	else basename++;
	ext = strrchr(basename, '.');
	if(!ext) return FALSE;
	else ext++;
	//DBG("basename = %s; ext = %s", basename, ext);
	if(strcasecmp(ext, "fts") == 0 || strcasecmp(ext, "fit") == 0 ||
		strcasecmp(ext, "fits") == 0) return TRUE;
	return FALSE;
}

// File test for loadable
gboolean is_loadable(gchar *filename){
    if(!Global->add_all){ // if we don't try to open any file
		if(!get_ext(filename)) // file suffix isn't fit/fits/fts
			return FALSE;
	}
    return guess_fits(filename);
}

gint add_file_to_list(gchar *filename, Window *window){
	if(is_loadable(filename)){
		DBG("File %s is FITS file\n", filename);
		window->files.files_list = g_list_prepend(window->files.files_list, (gpointer)filename);
		window->files.list_length++;
		if(window->files.list_length == 1)
			window->files.list_end = window->files.files_list;
		return 1;
	}
	return 0;
}

void free_filelist(Window *window){
	if(window->files.files_list){
		g_list_foreach(window->files.files_list, (GFunc)free, NULL);
		g_list_free(window->files.files_list);
		window->files.list_length = 0;
		window->files.files_list = NULL;
	}
}

/*
 * get a FITS files list in current directory
 *		filename - full filename (with path)
 * 		returns a files number in list
 */
gint fill_filelist(gchar *filename, Window *window){
	gint nb_inserted = 0;
	GDir *dir;
	GError *err = NULL;
	const gchar *dir_entry;
	gchar *dirname = g_path_get_dirname(filename);
	DBG("Open directory %s", dirname);
	dir = g_dir_open(dirname, 0, &err);
	if(dir == NULL){
		g_err(err->message);
		g_error_free(err);
		return 0;
	}
	free_filelist(window); // clear list if it was already created
	while((dir_entry = g_dir_read_name(dir))){
		gchar *full_path = g_build_filename(dirname, dir_entry, NULL);
		if(!g_file_test(full_path, G_FILE_TEST_IS_DIR)){
			nb_inserted += add_file_to_list(full_path, window);
		}
	}
	g_dir_close(dir);
	g_free(dirname);
	window->files.files_list = g_list_sort(window->files.files_list, (GCompareFunc)strcmp);
	DBG("\nSort %d items\n", nb_inserted);
	window->files.list_end = g_list_last(window->files.files_list);
	window->files.list_current = g_list_find_custom(window->files.files_list, filename,(GCompareFunc)strcmp);
//if(!window->files.list_current) window->files.list_current = g_list_first(window->files.files_list);
	//~ place_first = g_list_find_custom(window->files.files_list, filename, strcmp);
	//~ if(place_first && place_first->prev){
		//~ place_first->prev->next = place_first->next;
		//~ if(place_first->next)
			//~ place_first->next->prev = place_first->prev;
		//~ place_first->prev = NULL;
		//~ place_first->next = window->files.files_list;
		//~ window->files.files_list->prev = place_first;
		//~ window->files.files_list = place_first;
	//~ }
	return nb_inserted;
}

