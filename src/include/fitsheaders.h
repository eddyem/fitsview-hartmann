#ifndef _FITSHEADERS_H_
#define _FITSHEADERS_H_
#include "fitsview.h"
#include "gtk.h"

enum{			// key parameters field
	L_TYPENM,	// type name
	L_KEY,		// key name
	L_VAL,		// key value
	L_COMM,		// comment
	L_TYPE,		// key type (invisible column, we need you?)
	L_MAX		// amount of parameters
};

void edit_headers(Window *parent);
void add_key(	IMAGE *ima,
				gchar *keytype, gchar *keyname, gchar *keyval,
				gchar *keycomment, int keyfitstype);
void init_keylist(IMAGE *ima);
#endif
