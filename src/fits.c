//      fits.c - main functions to work with FITS files
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
#include "fits.h"
#include "gtk.h"
#include "fitsheaders.h"

/*
 * Macros for error processing when working with cfitsio functions
 */
#define TRYFITS(f, ...)						\
do{	gboolean status = FALSE;				\
	f(__VA_ARGS__, &status);				\
	if(status){								\
		fits_report_error(stderr, status);	\
		return FALSE;}						\
}while(0)
#define FITSFUN(f, ...)						\
do{	gboolean status = FALSE;				\
	int ret = f(__VA_ARGS__, &status);		\
	if(ret || status)						\
		fits_report_error(stderr, status);	\
}while(0)
#define WRITEKEY(...)							\
do{ gboolean status = FALSE;					\
	fits_write_key(__VA_ARGS__, &status);		\
	if(status) fits_report_error(stderr, status);\
}while(0)

/*
 * Check whether a number is a complex
 * If yes - returns true, no - false
 * If re and/or im!=NULL, in them stored values
 * of real and imaginary components
 */
gboolean cmplx_conv(gchar *value, double *re, double *im){
	double sign = 1.;
	double r,i;
	gchar *ptr1, *ptr2, *eptr;
	ptr1 = g_strstr_len(value, -1, "(");
	if(ptr1){ // form is (re, im) or (re im)
		ptr1++;
		ptr2 = g_strstr_len(ptr1, -1, ")");
		if(ptr2) *ptr2 = 0;
		ptr2 = g_strstr_len(ptr1, -1, ",");
		if(!ptr2)
			ptr2 = g_strstr_len(ptr1, -1, " ");
	}
	else{ // form is re + iim or re + jim
		ptr1 = value;
		ptr2 = g_strstr_len(ptr1, -1, "i");
		if(!ptr2)
			ptr2 = g_strstr_len(ptr1, -1, "j");
	}
	if(ptr2){
		*ptr2 = 0;
		ptr2++;
	}
	else
		return FALSE;
	if(g_strrstr_len(ptr1, 3, "-")) sign=-1.;
	r = strtod(ptr1, &eptr);
	if(eptr == ptr1) return FALSE;
	i = sign*strtod(ptr2, &eptr);
	if(eptr == ptr2) return FALSE;
	if(re) *re = r;
	if(im) *im = i;
	return TRUE;
}
/*
 * try to guess key type:
 * TLONGLONG, TDOUBLE, TDBLCOMPLEX or TSTRING
 */
int guess_type(gchar *keyval, gchar* keytype){
	char *eptr, *val;
	int dtype = TSTRING;
	val = strdup(keyval);
	if(strchr(val, '.')){
		errno = 0;
		strtod(val, &eptr);
		if(eptr != val && errno != ERANGE)
			dtype = TDOUBLE;
	}
	if(dtype == TSTRING){
		errno = 0;
		strtoll(val, &eptr, 0);
		if(eptr != val && errno != ERANGE) dtype = TLONGLONG;
	}
	if(dtype == TSTRING && cmplx_conv(val, NULL, NULL))
		dtype = TDBLCOMPLEX;
	if(keytype)	switch(dtype){
		case TDOUBLE:
			sprintf(keytype, "TDOUBLE");
		break;
		case TLONGLONG:
			sprintf(keytype, "TLONGLONG");
		break;
		case TDBLCOMPLEX:
			sprintf(keytype, "TDBLCOMPLEX");
		break;
		default:
			sprintf(keytype, "TSTRING");
	}
	//~ DBG("dtype: %d, keytype: %s", dtype, (keytype) ? keytype : "");
	free(val);
	return dtype;
}
/*
 * Save image-structure to a FITS file
 * (together whith keys)
 */
gboolean writefits(char *filename, IMAGE *image){
	FNAME();
	long naxes[2] = {image->width, image->height};
	int w = image->width, h = image->height;
	ssize_t sz = w * h;
	gchar *keytype, *keyname, *keyval, *keycomment; int keyfitstype;
	GtkTreeIter iter;
	fitsfile *fp;
	gchar *card;
	#ifdef EBUG
		int n = 0;
	#endif
	TRYFITS(fits_create_file, &fp, filename);
	TRYFITS(fits_create_img, fp, image->dtype, 2, naxes);
	if(image->store){
		GtkTreePath *path = gtk_tree_path_new_first();
		GtkTreeModel *tree_model = GTK_TREE_MODEL(image->store);
		if(gtk_tree_model_get_iter(tree_model, &iter, path)){ // there's keys
			DBG("filename: %s", filename);
			do{
				gtk_tree_model_get(tree_model, &iter,
					L_TYPENM,	&keytype,
					L_KEY,		&keyname,
					L_VAL,		&keyval,
					L_COMM,		&keycomment,
					L_TYPE,		&keyfitstype,
					-1);
				if(keyfitstype == TSTRING) // make sure, that string is in quotes
					if(*keyval != '\''){
						gchar *str = g_strdup_printf("'%s'", keyval);
						g_free(keyval);
						keyval = str;
					}
				card = g_strdup_printf("%-8s= %s / %s",
					keyname, keyval, keycomment);
				g_free(keytype); g_free(keyname); g_free(keyval);
				g_free(keycomment);
				DBG("key#%d: %s", ++n, card);
				FITSFUN(fits_write_record, fp, card);
				g_free(card);
			}while(gtk_tree_model_iter_next(tree_model, &iter));
		}
		gtk_tree_path_free(path);
	}
#ifdef EBUG
	GLfloat min, wd;
	min = image->stat.min;
	wd = image->stat.max - min;
	DBG("min = %f, wd = %f", min, wd);
#endif
/*	int i,j;
	GLfloat *newdata = malloc(sz*sizeof(GLfloat));
	GLfloat *ptro = newdata, *ptri = image->data;
	for(i=0; i<h; i++)
		for(j=0; j<w; j++, ptro++, ptri++){
			*ptro = min + (*ptri) * wd;
		}
	TRYFITS(fits_write_img, fp, TFLOAT, 1, sz, newdata);
	free(newdata);*/
	TRYFITS(fits_write_img, fp, TFLOAT, 1, sz, image->data);
	TRYFITS(fits_close_file, fp);
	return TRUE;
}

gboolean try_open_file(gchar *filename, IMAGE *ima){
	if(!readfits(filename, ima)){
		g_err(_("Can't read fits file"));
		ima->data = NULL;
		return FALSE;
	}
	return TRUE;
}
/*
char *del_spaces(char *str){
	char *beg, *end;
	beg = str;
	end = str + strlen(str);
	while(*beg == ' ' && beg < end) beg++;
	if(beg == end) return beg;
	end--;
	if(beg == end) return beg;
	while(*end == ' ' && end > beg)
		*end-- = 0;
	return beg;
}
gboolean strip_card(char *card, char *keyname, char *keyval, char *keycomment){
	char *ptr, *ptr1, *end;
	char *c = strdup(card);
	end = c + strlen(c);
	ptr = strchr(c, '=');
	if(!ptr) return FALSE;
	*ptr = 0; ptr++;
	ptr1 = del_spaces(c);
	strncpy(keyname, ptr1, FLEN_KEYWORD-1);
	ptr1 = strchr(ptr, '\'');
	if(ptr1){
		ptr1++;
		ptr = strchr(ptr1, '\'');
		if(!ptr){
			ptr = strchr(ptr1, '/');
			if(!ptr){
				free(c);
				return FALSE;
			}
			*ptr = 0; ptr++;
		}else{
			*ptr = 0; ptr = strchr(ptr+1, '/');
			if(!ptr){
				free(c);
				return FALSE;
			}
			*ptr = 0; ptr++;
		}
	}else{
		ptr1 = ptr;
		ptr = strchr(ptr1, '/');
		if(!ptr){
			free(c);
			return FALSE;
		}
		*ptr = 0; ptr++;
		ptr1 = del_spaces(ptr1);
	}
	if(ptr >= end) return FALSE;
	ptr1 = del_spaces(ptr1);
	strncpy(keyval, ptr1, FLEN_VALUE);
	ptr1 = del_spaces(ptr);
	strncpy(keycomment, ptr1, FLEN_COMMENT);
	free(c);
	return TRUE;
}*/

/*
 * Read FITS file filename to a structure image
 */
gboolean readfits(gchar *filename, IMAGE *image){
	FNAME();
	gboolean ret = TRUE;
	fitsfile *fp;
	GLfloat nullval = 0., imBits, bZero = 0., bScale = 1.;
	free(image->data);
	image->data = NULL;
	int i, j, hdunum, keynum, morekeys, dtype, stat;
	int naxis, w, h;
	double val_f = -1e6; // VAL_F value
	long naxes[4];
	ssize_t sz;
	//char card[FLEN_CARD];
	char keyname[FLEN_KEYWORD], keyval[FLEN_VALUE];
	char keycomment[FLEN_COMMENT], cdtype[32];
	TRYFITS(fits_open_file, &fp, filename, READWRITE);
/*	TRYFITS(fits_get_hdrspace, fp, &keynum, &morekeys);
	if(morekeys == -1){
		return FALSE;
	}
	g_print("Read header of %s, %d keys\n", filename, keynum);*/
	FITSFUN(fits_get_num_hdus, fp, &hdunum);
	//~ g_print(_("\nFile includes %d HDUs\n"), hdunum);
	if(hdunum < 1){
		g_err(_("Can't read HDU"));
		TRYFITS(fits_close_file, fp);
		return FALSE;
	}
	FITSFUN(fits_get_hdu_type, fp, &dtype);
	//~ if(dtype) g_print(_("Current HDU type: %d\n"), dtype);
	init_keylist(image);
	FITSFUN(fits_get_hdrpos ,fp, &keynum, &morekeys);

	for (j = 1; j <= keynum; j++){
		/*FITSFUN(fits_read_record, fp, j, card);
		if(!strip_card(card, keyname, keyval, keycomment)){
			g_print("Corrupted card: %s!\n", card);
			continue;
		}*/
		FITSFUN(fits_read_keyn, fp, j, keyname, keyval, keycomment);
		if(strcmp(keyname, "SIMPLE") == 0 || strcmp(keyname, "EXTEND") == 0) // key "file does conform ..."
			continue;
		else if(strcmp(keyname, "COMMENT") == 0) // comment of obligatory key in FITS head
			continue;
		else if(strstr(keyname, "NAXIS") == keyname || strcmp(keyname, "BITPIX") == 0) // NAXIS, NAXISxxx, BITPIX
			continue;
		else if(strcmp(keyname, "BZERO") == 0){
			bZero = atof(keyval);
			//continue;
		}
		else if(strcmp(keyname, "BSCALE") == 0){
			bScale = atof(keyval);
			if(bScale < SCALE_MIN) bScale = 1.;
			//continue;
		}
		else if(strcmp(keyname, "VAL_F") == 0){
			val_f = atof(keyval);
		}
		dtype = guess_type(keyval, cdtype);
		//DBG("key#%d; name=%s, dtype=%d (%s), value=%s, comment=%s\n",
		//		j, keyname, dtype, cdtype, keyval, keycomment);
		add_key(image, cdtype, keyname, keyval, keycomment, dtype);
	}
	if(hdunum > 1){
		for(i = 2; !(fits_movabs_hdu(fp, i, &dtype, &stat)); i++){
			FITSFUN(fits_get_hdrpos ,fp, &keynum, &morekeys);
			for (j = 1; j <= keynum; j++){
				FITSFUN(fits_read_keyn, fp, j, keyname, keyval, keycomment);
				dtype = guess_type(keyval, cdtype);
				add_key(image, cdtype, keyname, keyval, keycomment, dtype);
			}
		}
		if(stat != END_OF_FILE){
			fits_report_error(stderr, (stat));
			ret = FALSE;
		}
	}

	TRYFITS(fits_get_img_param, fp, 4, &dtype, &naxis, naxes);
	if(dtype > 0) imBits = (GLfloat)(1LL << dtype); // compute the max amplitude
	else imBits = -1.; // image has floating point type
	DBG("Image type: %d (%f bits)", dtype, imBits);
	image->maxVal = imBits;
	image->dtype = dtype;
	if(naxis != 2){
		g_err(_("Not an image? (dimensions != 2)"));
		return FALSE;}

	image->val_f = val_f;
	image->bZero = bZero;
	image->bScale = bScale;
	image->width = w = naxes[0];
	image->height = h = naxes[1];
	sz = (ssize_t)(w) * (ssize_t)(h);
	image->data = malloc(sz * sizeof(GLfloat));

	TRYFITS(fits_read_img, fp, TFLOAT, 1, sz, &nullval, image->data, &stat);
	GLfloat *ptr = image->data;
	GLfloat min, max;
	min = 1e6; max = 0.;
	for(i=0; i<h; i++)
		for(j=0; j<w; j++, ptr++){
			GLfloat tmp = *ptr;
			if(tmp > max) max = tmp;
			else if(tmp < min) min = tmp;
		}
	DBG("min = %f, wd = %f", min, max-min);
	image->stat.max = max;
	image->stat.min = min;
	// histogram isn't initialized yet
	image->stat.histogram.data = NULL;
	image->stat.histogram.size = 0;
	/*
	ptr = image->data;
	for(i=0; i<h; i++)
		for(j=0; j<w; j++, ptr++){
			// *ptr = (*ptr * bScale + bZero) / imBits;
			*ptr = (*ptr - min) / wd;
		}
	*/
	TRYFITS(fits_close_file, fp);
	return(ret);
}

// Try to open FITS and find an image in it
gboolean guess_fits(gchar *filename){
	fitsfile *fp;
	long naxes[4];
	int naxis, dtype;
	DBG("Try to open file %s\n", filename);
	TRYFITS(fits_open_file, &fp, filename, READONLY);
	TRYFITS(fits_get_img_param, fp, 4, &dtype, &naxis, naxes);
	if(naxis != 2) return FALSE;
	TRYFITS(fits_close_file, fp);
	return TRUE;
}
