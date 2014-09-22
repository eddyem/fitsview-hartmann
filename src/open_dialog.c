//      open_dialog.c - open file dialog with preview
//      Copyright:
//      Guillaume Chazarain <guichaz@gmail.com> (gliv project)
//      2011 Edward V. Emelianoff <eddy@sao.ru> (modification for fits open)
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
#include "open_dialog.h"

gchar *saved_path = NULL;

typedef struct{
	GtkLabel  *label;
	GtkWidget *image;
} fits_preview;

GdkPixbuf *getfits(char *filename, gchar **descr);

void update_preview(	GtkFileChooser *file_chooser,
								fits_preview *FP){
	GtkWidget *preview = FP->image;
	GtkLabel *ilabel = FP->label;
	char *filename;
	gchar *descr = NULL;
	GdkPixbuf *pixbuf = NULL;
	filename = gtk_file_chooser_get_preview_filename(file_chooser);
	if (filename == NULL) return;
	pixbuf = getfits(filename, &descr);
	if(pixbuf){
		gint w = gdk_pixbuf_get_width(pixbuf);
		w = w - w/4;
		if(w < 200) w = 200;
		gtk_widget_set_size_request(GTK_WIDGET(ilabel), w, -1);
		gtk_label_set_text(ilabel, descr);
	}
	g_free(filename);
	gtk_image_set_from_pixbuf(GTK_IMAGE(preview), pixbuf);
	if(pixbuf) g_object_unref(pixbuf);
	gtk_file_chooser_set_preview_widget_active(file_chooser, (pixbuf != NULL));
}

void cb_changed(GtkComboBox *combo, GtkObject *chooser){
	GtkTreeIter iter;
	GtkTreeModel *model;
	int val;
	gtk_combo_box_get_active_iter(combo, &iter);
	model = gtk_combo_box_get_model(combo);
	gtk_tree_model_get(model, &iter, 1, &val, -1);
	// clear needed config bits
	if(val & PREVIEW_SIZEMASK) Global->previewMode &= ~PREVIEW_SIZEMASK;
	else if(val & PREVIEW_COLFNMASK) Global->previewMode &= ~PREVIEW_COLFNMASK;
	Global->previewMode |= val; // set new
	gtk_signal_emit_by_name(chooser, "update-preview"); // redraw
}

// arguments of funtion which creates menu
enum{
	COMBO_SIZE,
	COMBO_COLORFN
};
GtkWidget *create_combo(unsigned char type, GtkObject *chooser){
	#define SIZES     5
	#define COLORFNS  3
	gchar *sizes[] = {"128 x 128", "256 x 256", "512 x 512",
						"768 x 768", "1024 x 1024"};
	gchar *colorfns[] = {N_("linear"), N_("log"), N_("square root")};
	gchar **ptr;
	GtkWidget *label = NULL, *vbox = NULL;
	gchar *label_text = NULL;
	int curpos, i, n;
	int startpos = 0; // bit position of appropriate numbers
	GtkListStore *store;
	GtkWidget *combo;
	GtkCellRenderer *cell;
	GtkTreeIter iter;
	if(type == COMBO_SIZE){
		ptr = sizes; n = SIZES;
		if((curpos = Global->previewMode & PREVIEW_SIZEMASK))
			curpos--;
		else
			curpos = PREVIEW_512 - 1;
		startpos = 0;
		label_text = g_strdup_printf(_("Preview size"));
	}
	else if(type == COMBO_COLORFN){
		ptr = colorfns; n = COLORFNS;
		if((curpos = (Global->previewMode & PREVIEW_COLFNMASK) >> 3))
			curpos--;
		else
			curpos = (PREVIEW_SQRT>>3) - 1;
		startpos = 3;
		label_text = g_strdup_printf(_("Colormap function"));
	}
	else return NULL;
	vbox = gtk_vbox_new(FALSE, 3);
	label = gtk_label_new(label_text);
	gtk_box_pack_start(GTK_BOX(vbox), label, FALSE, FALSE, 0);
	store = gtk_list_store_new(2, G_TYPE_STRING, G_TYPE_INT);
	for(i = 0; i < n; i++){
		gtk_list_store_append(store, &iter);
		gtk_list_store_set(store, &iter, 0, _(ptr[i]), 1, (i+1)<<startpos, -1);
	}
	combo = gtk_combo_box_new_with_model(GTK_TREE_MODEL(store));
	gtk_box_pack_start(GTK_BOX(vbox), combo, FALSE, FALSE, 0);
	cell = gtk_cell_renderer_text_new();
	gtk_cell_layout_pack_start(GTK_CELL_LAYOUT(combo), cell, TRUE);
	gtk_cell_layout_set_attributes(GTK_CELL_LAYOUT(combo), cell, "text", 0, NULL);
	gtk_combo_box_set_active(GTK_COMBO_BOX(combo), curpos);
	g_signal_connect(combo, "changed", G_CALLBACK(cb_changed), (gpointer)chooser);
	g_object_unref(store);
	free(label_text);
	return vbox;
}

//GtkFileChooser *open_fits_dialog(gboolean select_dir){
GtkFileChooser *open_fits_dialog(){
	GtkLabel *ilabel = NULL;
	GtkImage *preview = NULL;
	GtkWidget *vbox = NULL, *combo = NULL, *frame = NULL, *sbox = NULL;
	GtkFileChooserAction action;
	const gchar *label;
	GtkFileChooser *chooser;
	static fits_preview FP;
/*	if (select_dir) {
		action = GTK_FILE_CHOOSER_ACTION_SELECT_FOLDER;
		label = _("Select a folder to open");
	} else { */
		action = GTK_FILE_CHOOSER_ACTION_OPEN;
		label = _("Select fits file to open");
	//}
	chooser = GTK_FILE_CHOOSER(
		gtk_file_chooser_dialog_new(label, NULL,
									action,
									GTK_STOCK_CANCEL,
									GTK_RESPONSE_CANCEL,
									GTK_STOCK_OPEN,
									GTK_RESPONSE_ACCEPT,
									NULL)
								);
	// it CAN CAUSED MEMORY LEAKS: filter doesn't deletes
	GtkFileFilter *filter = gtk_file_filter_new();
	gtk_file_filter_set_name(filter, _("FITS files"));
	gtk_file_filter_add_pattern(filter, "*.fit");
	gtk_file_filter_add_pattern(filter, "*.FIT");
	gtk_file_filter_add_pattern(filter, "*.fts");
	gtk_file_filter_add_pattern(filter, "*.FTS");
	gtk_file_filter_add_pattern(filter, "*.fits");
	gtk_file_filter_add_pattern(filter, "*.FITS");
	gtk_file_chooser_add_filter(chooser, filter);
	filter = gtk_file_filter_new();
	gtk_file_filter_set_name(filter, _("All files"));
	gtk_file_filter_add_pattern(filter, "*");
	gtk_file_chooser_add_filter(chooser, filter);

	if(saved_path)
		gtk_file_chooser_set_current_folder(chooser, saved_path);
//	if (!select_dir) {
	vbox = gtk_vbox_new(FALSE, 8);
	ilabel = GTK_LABEL(gtk_label_new(NULL));
	gtk_label_set_line_wrap(ilabel, TRUE);
	gtk_label_set_max_width_chars(ilabel, 80);
	gtk_label_set_line_wrap_mode(ilabel, PANGO_WRAP_WORD);
	// either wrap or angle
	//gtk_label_set_angle(ilabel, 30.);
	preview = GTK_IMAGE(gtk_image_new());
	gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(preview),
					FALSE, FALSE, 0);
	gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(ilabel),
						FALSE, FALSE, 0);
	FP.label = ilabel;  FP.image = GTK_WIDGET(preview);
	//gtk_file_chooser_set_select_multiple(chooser, TRUE);
	gtk_widget_show_all(GTK_WIDGET(vbox));
	gtk_file_chooser_set_preview_widget(chooser, GTK_WIDGET(vbox));
	vbox = gtk_hbox_new(FALSE, 10);
	frame = gtk_frame_new(_("Preview settings"));
	gtk_box_pack_end(GTK_BOX(vbox), frame, FALSE, FALSE, 0);
	sbox = gtk_hbox_new(FALSE, 10);
	gtk_container_add(GTK_CONTAINER(frame), sbox);
	combo = create_combo(COMBO_SIZE, GTK_OBJECT(chooser));
	gtk_box_pack_start(GTK_BOX(sbox), combo, FALSE, FALSE, 0);
	combo = create_combo(COMBO_COLORFN, GTK_OBJECT(chooser));
	gtk_box_pack_start(GTK_BOX(sbox), combo, FALSE, FALSE, 0);
	gtk_widget_show_all(frame);
	gtk_file_chooser_set_extra_widget(chooser, vbox);
	g_signal_connect(chooser, "update-preview",
				G_CALLBACK(update_preview), (gpointer)&FP);
//	}
	return chooser;
}
/*
// gray - intensity from 0.f to 1.f
void gray2rgb(float gray, guchar *rgb){
	int i = (int)(gray * 4.);
	float x = gray - (float)i * .25;
	guchar r = 0, g = 0, b = 0;
	switch(i){
		case 0:
			b = (guchar)(255. * x);
		break;
		case 1:
			b = 255;
			g = (guchar)(255. * x);
		break;
		case 2:
			r = (guchar)(255. * x);
			g = 255;
			b = (guchar)(255. * (1. - x));
		break;
		case 3:
			r = 255;
			g = (guchar)(255. * (1. - x));
		break;
		default:
			if(i>3) r=255;
	}
	*rgb++ = r;
	*rgb++ = g;
	*rgb++ = b;
}*/
void gray2rgb(float gray, guchar *rgb){
	int i = (int)(gray * 4.);
	float x = gray - (float)i * .25;
	guchar r = 0, g = 0, b = 0;
	switch(i){
		case 0:
			g = (guchar)(255. * x);
			b = 255;
		break;
		case 1:
			g = 255;
			b = (guchar)(255. * (1. - x));
		break;
		case 2:
			r = (guchar)(255. * x);
			g = 255;
		break;
		case 3:
			r = 255;
			g = (guchar)(255. * (1. - x));
		break;
		default:
			r = 255;
	}
	*rgb++ = r;
	*rgb++ = g;
	*rgb++ = b;
}

GdkPixbufDestroyNotify free_preview_data(guchar *pixels, gpointer data){
	free(pixels);
	free(data);
	return FALSE;
}

// get preview from FITS file filename
GdkPixbuf *getfits(char *filename, gchar **descr){
	gboolean status;
	fitsfile *fp;
	double (*colorfun)(double);
	double linfun(double arg){ return arg; } // bung for PREVIEW_LINEAR
	double logfun(double arg){ return log(1.+arg); } // for PREVIEW_LOG
	gchar *description = NULL; // some keywords from FITS header
	#define TRYFITS(f, ...)						\
	do{	status = FALSE;				\
		f(__VA_ARGS__, &status);				\
		if(status){								\
			free(ima_data); free(pixbuf_data);	\
			free(description); free(pix);		\
			fits_close_file(fp, &status);		\
			return NULL;}						\
	}while(0)
	void add_keyw(char *keyw){
		char keyval[FLEN_VALUE], *ptr;
		if(VALUE_UNDEFINED == fits_read_key(fp,
				TSTRING, keyw, keyval, NULL, &status)) return;
		if(status) return;
		ptr = g_strdup_printf("%s;\t%s=%s", description, keyw, keyval);
		free(description);
		description = ptr;
	}
	int MAX_SIZE = 512; // max preview size
	unsigned char cntxt;
	GdkPixbuf *pixbuf = NULL;
	float nullval = 0.;
	int i, j, k, l, N, M, stat;
	int naxis, w, h, pixScale, Ws, Hs, dtype;
	int sz;
	cntxt = Global->previewMode & PREVIEW_SIZEMASK; // get preview size
	if(cntxt)
		switch(cntxt){
			case PREVIEW_128: MAX_SIZE = 128;
			break;
			case PREVIEW_256: MAX_SIZE = 256;
			break;
			case PREVIEW_512: MAX_SIZE = 512;
			break;
			case PREVIEW_768: MAX_SIZE = 768;
			break;
			case PREVIEW_1024:MAX_SIZE = 1024;
			break;
		}
	// array for preview picture line
	float *pix = malloc(MAX_SIZE * sizeof(float));
	long naxes[4];
	float *ima_data = NULL, *ptr, byte, n, m, max, min, wd, avr;
	guchar *pptr, *pixbuf_data = NULL;
	DBG("Try to open file %s\n", filename);
	TRYFITS(fits_open_file, &fp, filename, READONLY);
	TRYFITS(fits_get_img_param, fp, 4, &dtype, &naxis, naxes);
	if(naxis != 2) return NULL;
	w = naxes[0];
	h = naxes[1];
	sz = w * h;
	ima_data = malloc(sz * sizeof(float));
	pixbuf_data = malloc(3 * MAX_SIZE * MAX_SIZE * sizeof(guchar));

	TRYFITS(fits_read_img, fp, TFLOAT, 1, sz, &nullval, ima_data, &stat);
	ptr = ima_data;
	min = max = *ptr; avr = 0.;
	// get statistics:
	for(i=0; i<h; i++)
		for(j=0; j<w; j++, ptr++){
			float tmp = *ptr;
			if(tmp > max) max = tmp;
			else if(tmp < min) min = tmp;
			avr += tmp;
		}
	avr /= (float)sz;
	wd = max - min;
	i = (int)ceil((float)w / MAX_SIZE);
	j = (int)ceil((float)h / MAX_SIZE);
	DBG("i=%d, j=%d, ms=%d",i,j,MAX_SIZE);
	pixScale = (i > j) ? i : j;	// picture scale factor
	Ws = w / pixScale; 			// picture width in pixScale blocks
	Hs = h / pixScale; 			// -//- height pixScale
	DBG("w=%d, h=%d, Ws=%d, Hs=%d, pixScale=%d",w,h,Ws,Hs,pixScale);
	// prepare a comment to a prewiew:
	description = g_strdup_printf(
					"%s=%dx%d\n%s=%dx%d\n%s=%.1f%%;\tMax=%g,\tMin=%g,\tAvr=%g",
					_("Image size"), w, h,
					_("Preview size"), Ws, Hs,
					_("Preview scale"), 100./(double)pixScale,
					max, min, avr
					);
	add_keyw("BITPIX");
	add_keyw("IMAGETYP");
	add_keyw("OBJECT");
	add_keyw("EXPTIME"); add_keyw("EXP");
	add_keyw("AUTHOR");
	add_keyw("DATE");
	M = 0; // line number
	for(i = 0; i < Hs; i++){ // cycle through a blocks by lines
		//pptr = &pixbuf_data[i * Ws * 3];
		for(j = 0; j < MAX_SIZE; j++) pix[j] = 0;
		m = 0.; // amount of strings read in block
		for(l = 0; l < pixScale; l++, m++){ // cycle through a block lines
			ptr = &ima_data[M * w];
			N = 0; // number of column
			for(j = 0; j < Ws; j++){ // cycle through a blocks by columns
				n = 0.;	// amount of columns read in block
				byte = 0.; // average intensity in block
				for(k = 0; k < pixScale; k++, n++){ // cycle through block pixels
					if(N++ < w) // row didn't end
						byte += *ptr++; // sum[(pix-min)/wd]/n = [sum(pix)/n-min]/wd
					else break;
				}
				pix[j] += byte / n;//(byte / n - min)/wd;
			}
			if(++M >= h) break;
		}
		// fill unused picture pixels
		ptr = &ima_data[i*Ws];
		for(l = 0; l < Ws; l++)
			*ptr++ = pix[l] / m;
	}
	ptr = ima_data;
	sz = Ws * Hs;
	max = min = *ptr;
	avr = 0;
	for(i=0; i < sz; i++, ptr++){
		float tmp = *ptr;
		if(tmp > max) max = tmp;
		else if(tmp < min) min = tmp;
		avr += tmp;
	}
	avr /= (float)sz;
	wd = max - min;
	avr = (avr - min) / wd;	// normal average by preview
	avr = -log(avr);		// scale factor
	if(avr > 1.) wd /= avr;
	ptr = ima_data;
	colorfun = sqrt;
	if((Global->previewMode & PREVIEW_COLFNMASK) == PREVIEW_LINEAR)
		colorfun = linfun;
	else if((Global->previewMode & PREVIEW_COLFNMASK) == PREVIEW_LOG)
		colorfun = logfun;
	for(i = Hs - 1; i > -1; i--){// fill pixbuf mirroring image by vertical
		pptr = &pixbuf_data[Ws * i * 3];
		for(j = 0; j < Ws; j++){
			gray2rgb(colorfun((*ptr++ - min) / wd), pptr);
			pptr += 3;
		}
	}
	fits_close_file(fp, &status);
	pixbuf = gdk_pixbuf_new_from_data(
				pixbuf_data,		// guchar* data
				GDK_COLORSPACE_RGB,	// only this supported
				FALSE,				// no alpha
				8,					// number of bits in byte (WTF? who is this idiot?)
				Ws, Hs,				// size
				Ws * 3,				// line length in bytes
				(GdkPixbufDestroyNotify)free_preview_data, // function (*GdkPixbufDestroyNotify) (guchar *pixels, gpointer data);
				(gpointer)description	// pointer data
			);
	free(ima_data);
	free(pix);
	*descr = description;
	DBG("OK, preview is ready, previewMode = %d", Global->previewMode);
	return pixbuf;
}
