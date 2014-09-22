//      gtk.c - main funtions to work with windows
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
#include "gtk.h"
#include "tracking.h"
#include "opengl.h"
#include "ui.h"
#include "imtools.h"
#include "contours.h"
#include "fits.h"
#include "open_dialog.h"
#include "fitsheaders.h"
#include "terrain.h"
#include "filelist.h"
#include <gdk/gdkkeysyms.h>

Window *mainWindow;

double X0, Y0, lX0, lY0; // started coordinates of click, layer coordinates at click
// the longer key is pressed, the larger value of key_acceleration (3D moving):
double key_acceleration = 1.;
// multipliers for zoom in/out
const double ZOOM_NEAR = 1.1, ZOOM_FAR = 1./1.1;

void run_modal_window(GtkWindow *w, Window *parent){
	gtk_window_set_transient_for(w, GTK_WINDOW(parent->window));
	gtk_window_set_modal(w, TRUE);
	gtk_window_set_position(w, GTK_WIN_POS_MOUSE);
	gtk_widget_show_all(GTK_WIDGET(w));
}

gint run_modal_dialog(GtkDialog *dialog, Window *parent){
	run_modal_window(GTK_WINDOW(dialog), parent);
	return gtk_dialog_run(GTK_DIALOG(dialog));
}

/*
 * Error messages output
 * text - message text
 */
void g_err(gchar *text){
	GtkMessageDialog *dialog;
	g_printerr(_("Error: %s\n"), text);
	dialog = GTK_MESSAGE_DIALOG(gtk_message_dialog_new(NULL,
								GTK_DIALOG_MODAL,
								GTK_MESSAGE_ERROR,
								GTK_BUTTONS_CLOSE,
								"%s", text));
	run_modal_dialog(GTK_DIALOG(dialog), mainWindow);
	gtk_widget_destroy(GTK_WIDGET(dialog));
}

void show_about_dialog(){
	GtkAboutDialog *about = GTK_ABOUT_DIALOG(gtk_about_dialog_new());
	const gchar *author[] = { "Edward V. Emelianov <eddy@sao.ru>, <edward.emelianoof@gmail.com>", NULL };
	gtk_about_dialog_set_program_name(about, "FITS view & Hartman data reduction");
	gtk_about_dialog_set_version(about, PACKAGE_VERSION);
	gtk_about_dialog_set_copyright(about, "Copyright 2012 Edward V. Emelianov");
	gtk_about_dialog_set_website(about, "http://eddyem.narod.ru");
	gtk_about_dialog_set_authors(about, author);
	run_modal_dialog(GTK_DIALOG(about), mainWindow);
	gtk_widget_destroy(GTK_WIDGET(about));
}

/*
 * exit
 */
void chk_quit(){
	gtk_main_quit();
	// uncomment next, when autosaving of settings will be released
	/*GtkMessageDialog *dialog;
	dialog = GTK_MESSAGE_DIALOG(gtk_message_dialog_new(NULL,
								GTK_DIALOG_MODAL,
								GTK_MESSAGE_QUESTION,
								GTK_BUTTONS_YES_NO,
								_("You sure, you want exit?")));
	gint ans = run_modal_dialog(GTK_DIALOG(dialog), mainWindow);
	gtk_widget_destroy(GTK_WIDGET(dialog));
	if(ans == GTK_RESPONSE_YES){
		gtk_main_quit();
	}*/
}

void change_window_filename(gchar *filename, Window *window){
	gchar *title, *bsnm;
	if(!filename) return;
	g_free(window->image->filename);
	bsnm = g_path_get_basename(filename);
	if(window == mainWindow)
		title = g_strdup_printf("%s (main) - %s", PRGNAME, bsnm);
	else
		title = g_strdup_printf("%s - %s", PRGNAME, bsnm);
	gtk_window_set_title(GTK_WINDOW(window->window), title);
	window->image->filename = bsnm;
	g_free(title);
}

void change_image(gchar *filename, Window *window){
	if(!filename) return;
	if(try_open_file(filename, window->image)){
		change_window_filename(filename, window);
	}else
		change_window_filename("", window);
	window->context->current_image = window->image;
	gen_texture(window->image, window, FALSE);
	window->context->transformMask = TRANS_NONE;
}

// "save as" filename choosing dialog
gchar *get_save_filename(Window *window){
	FNAME();
	GtkFileChooser *dialog;
	gchar *filename = NULL, *retfilename = NULL;
	gint response;
	dialog = GTK_FILE_CHOOSER(
				gtk_file_chooser_dialog_new("Save File", NULL,
						GTK_FILE_CHOOSER_ACTION_SAVE,
						GTK_STOCK_CANCEL, GTK_RESPONSE_CANCEL,
						GTK_STOCK_SAVE, GTK_RESPONSE_ACCEPT,
						NULL));
	if(saved_path)
		gtk_file_chooser_set_current_folder(dialog, saved_path);
	gtk_file_chooser_set_current_name(dialog, "new_fits_file.fit");
	gtk_file_chooser_set_do_overwrite_confirmation(dialog, TRUE);
	response = run_modal_dialog(GTK_DIALOG(dialog), window);
	if(response == GTK_RESPONSE_ACCEPT){
		filename = gtk_file_chooser_get_filename(dialog);
	}
	gtk_widget_destroy(GTK_WIDGET(dialog));
	if(filename){
		if(!get_ext(filename)) // There wasn't .fit suffix
			retfilename = g_strdup_printf("!%s.fit", filename); // to rewrite file, if needed
		else
			retfilename = g_strdup_printf("!%s", filename);
		g_free(filename);
	}
	return retfilename;
}

// get choosed filename to open and fills filelist of current directory
gchar *get_open_filename(Window *window){
	FNAME();
	GtkFileChooser *dialog;
	gchar *filename = NULL;
	gint response;
	dialog = open_fits_dialog();
	response = run_modal_dialog(GTK_DIALOG(dialog), window);
	g_free(saved_path);
	saved_path = gtk_file_chooser_get_current_folder(dialog);
	if(response == GTK_RESPONSE_ACCEPT){
		filename = gtk_file_chooser_get_filename(dialog);
		DBG("open file %s", filename);
		gint nfiles = fill_filelist(filename, window);
		DBG("found %d files", nfiles);
		g_free(filename); filename = NULL;
		if(nfiles > 0 && window->files.list_current)
			filename = window->files.list_current->data;
		else
			g_err(_("Can't read fits file"));
	}
	gtk_widget_destroy(GTK_WIDGET(dialog));
	return filename;
}

void save_as(Window *window){
	gchar *filename = NULL;
	if((filename = get_save_filename(window))){
		change_window_filename(filename, window);
		if(!writefits(filename, window->image))
			g_err(_("Can't save file"));
		g_free(filename);
	}
}
void open_file(Window *window){
	gchar *filename;
	if((filename = get_open_filename(window)))
		change_image(filename, window);
}
void open_in_new(Window *window){
	gchar *filename;
	Window *newwin = init_window(window, OPENGL_WINDOW);
	if(!(filename = get_open_filename(newwin)))
		destroy_window(newwin);
	else
		change_image(filename, newwin);
}
void open_in_3D(Window *window){
	Window *terra = init_window(window, GL3D_WINDOW);
	terrain_3D(terra, 0.,0.,0.,0., TRUE);
}

void close_subwindow(Window *w){
	gtk_widget_hide(w->window);
}

// Refreshing State fields of status bar
void refresh_state(Window *window){
	//FNAME();
	Context *c = window->context;
	gchar buf[32];
	gchar draw[] = {'T', 'S', 'Q'};
	gchar track[] = {'A', 'V', 'H', 'D'};
	gchar graph[] = {'G', 'H'};
	if(window->id != GRAPH_WINDOW){
		g_sprintf(buf, "dX=%.1f, dY=%.1f", window->move.x, window->move.y);
		set_status_text(StatusAdd, buf, window);
		g_sprintf(buf, "%c %c %s (%.1f%%) Z=%.1f",
			draw[c->drawingMode], track[c->trackingMode],
			(c->transformMask ? "tr" : ""),
			100./window->Daspect*window->Zoom, window->Zoom);
		set_status_text(StatusState, buf, window);
	}else{
		g_sprintf(buf, "%c %s",
			graph[c->graphMode], (isYlog(window) ? "log" : "lin"));
		set_status_text(StatusState, buf, window);
	}
}

void set_status_text(guint barName, gchar *text, Window *window){
	if(!window || barName >= window->statusBlocks) return;
	gtk_entry_set_text(window->SBars[barName], text);
}

/*
 * This three short functions need to have ability to change scale
 * of axe Y from menu
 */
gboolean redraw_D(Window *window){
	Gconfigure(window);
	gdk_window_invalidate_rect(window->drawingArea->window, NULL, FALSE);
	refresh_state(window);
	return TRUE;
}
gboolean log_scale(Window *window){ // press 'e' in graph window
	window->context->graphYaxis |=
		(window->context->graphMode == GR_GRAPH) ? Y_LOGGRAPH : Y_LOGHIST;
	redraw_D(window);
	return TRUE;
}
gboolean lin_scale(Window *window){ // press 'l' in graph window
	window->context->graphYaxis &= (window->context->graphMode == GR_GRAPH) ?
						~Y_LOGGRAPH : ~Y_LOGHIST;
	redraw_D(window);
	return TRUE;
}

// Show histogram in chart window
void show_histogram(Window *window){
	IMAGE *ima;
	if(!window->context->current_image || !window->context->current_image->data){
		window->context->current_image = window->image;
		if(!window->image || !window->image->data) return;
	}
	ima = window->context->current_image;
	if(!ima->stat.histogram.data)
		if(!make_histogram(ima))
			return;
	gtk_toggle_action_set_active(GTK_TOGGLE_ACTION(window->graphWindow->LinLogMenu[isYlog(window)]), TRUE);
	gtk_widget_show_all(window->graphWindow->window);
	do_histogram(window);
}

void correct_zoom_by_move(Window *window, double Zoom_factor, double *x, double *y){
	double xx, yy, Xc, Yc, ZZ;
	if(!x) xx = (double)window->drawingArea->allocation.width / 2.;
	else xx = *x;
	if(!y) yy = (double)window->drawingArea->allocation.height / 2.;
	else yy = *y;
	//DBG("x=%g, y=%g",xx,yy);
	conv_mouse_to_image_coords(xx, yy, &Xc, &Yc, window);
	Xc -= window->texture->w / 2.; Yc -= window->texture->h / 2.;
	//DBG("Xc=%g, Yc=%g, ZF=%g",Xc, Yc, Zoom_factor);
	// It's amazing (if U don't know what float is), but  Zoom*(Zoom_factor - 1) isn't the same
	ZZ = window->Zoom*Zoom_factor - window->Zoom;
	window->move.x -= Xc * ZZ;
	window->move.y -= Yc * ZZ;
	//DBG("moveX=%g, moveY=%g\n", window->move.x, window->move.y);
	window->Zoom *= Zoom_factor;
}

void zoom_near(Window *window){ // run for '+' key
	correct_zoom_by_move(window, ZOOM_NEAR, NULL, NULL);
	force_redraw(window->drawingArea);
}
void zoom_far(Window *window){ // run for '-' key
	correct_zoom_by_move(window, ZOOM_FAR, NULL, NULL);
	force_redraw(window->drawingArea);
}

void zoom_selection(Window *window){ // function called after frame selection
	double x0,y0, x,y, zX, zY, ZF;
	get_sel_region(&x0, &y0, &x, &y);
	zX = (double)window->drawingArea->allocation.width / abs(x-x0);
	zY = (double)window->drawingArea->allocation.height / abs(y-y0);
	ZF = (zX < zY) ? zX : zY;
	x = (x + x0) / 2.;
	y = (y + y0) / 2.;
	correct_zoom_by_move(window, ZF, &x, &y);
	force_redraw(window->drawingArea);
}

void show_3D_terrain(Window *window){
	double imX0, imY0, imX, imY;
	double x0,y0, x,y;
	get_sel_region(&imX0, &imY0, &imX, &imY);
	conv_mouse_to_image_coords(imX0, imY0, &x0, &y0, window);
	conv_mouse_to_image_coords(imX, imY, &x, &y, window);
	Window *terra = init_window(window, GL3D_WINDOW);
	terrain_3D(terra, x0,y0, x,y, FALSE);
}

void full_3D_terrain(Window *window){
	Window *terra = init_window(window, GL3D_WINDOW);
	terrain_3D(terra, 0,0, window->image->width,window->image->height, FALSE);
}

gboolean is_tracking = FALSE; // TRUE is cutting mode
gboolean zoom_frame = FALSE;  // TRUE is zoom frame mode
gboolean select_3D = FALSE;   // TRUE is 3D frame selection mode
/*
 * Processing of mouse button press event
 * 1. Main window
 */
gboolean press_main_mbtn(Window *window, GdkEventButton *event){
	Context *c = window->context;
	switch(event->button){
		case 1: // Left button? Start cutting
			// Redraw window to delete old garbage
			if(!window->image->data) break;
			switch(c->drawingMode){
			case DRAW_TRACK:
				gtk_toggle_action_set_active(GTK_TOGGLE_ACTION(window->graphWindow->LinLogMenu[isYlog(window)]), TRUE);
				gtk_widget_show_all(window->graphWindow->window);
				c->graphMode = GR_GRAPH;
				if(event->state & GDK_CONTROL_MASK)
						c->trackingMode = TRACK_DIAG;
				else if(event->state & GDK_SHIFT_MASK)
						c->trackingMode = TRACK_VERT;
				else if(event->state & GDK_MOD1_MASK)
						c->trackingMode = TRACK_HORIZ;
				else
						c->trackingMode = TRACK_ANY;
			case DRAW_QUAD: // There's no break here and it's not an error!
				is_tracking = TRUE;
				do_tracking(event->x, event->y, TRUE, window);
			break;
			case DRAW_SELECTION:
				if(window->image->imagetype == HOUGH){
					double R, phi;
					get_houg_line(window, event->x, event->y, &phi, &R);
					gchar *text = g_strdup_printf("%s: R=%.0f, phi=%.1f",
						_("Selected line"), R, phi);
					set_status_text(StatusText, text, window);
					g_free(text);
				}
				else
					pickRegion(event->x, event->y);
			break;
			}
			gdk_window_invalidate_rect(window->drawingArea->window,
					&window->drawingArea->allocation, FALSE);
		break;
		case 2: // Middle button moves an picture
			X0 = event->x; Y0 = event->y;
			lX0 = window->move.x; lY0 = window->move.y;
		break;
		case 3:
		break;
	}
	refresh_state(window);
	return FALSE;
}
// 2. Auxiliary  window
gboolean press_sub_mbtn(Window *window, GdkEventButton * event){
	switch(event->button){
		case 1:
			is_tracking = graph_mouse_btn((int)event->x, window);
		break;
		case 2:
		break;
		case 3:
		break;
	}
	refresh_state(window);
	return FALSE;
}

/*
 * Reaction for button release
 * 1. Main window
 */
gboolean release_main_mbtn(Window *window, GdkEventButton * event){
	switch(event->button){
		case 1:
			// there was something selected by frame
			if(is_tracking && window->context->drawingMode == DRAW_QUAD){
				if(zoom_frame){ // It was a zoom selection
					zoom_selection(window);
				}else if(select_3D){
					show_3D_terrain(window);
				}
			}
		break;
		case 2:
		break;
		case 3:
		break;
	}
	return FALSE;
}

// 2. Auxiliary  window
gboolean release_sub_mbtn(Window *window, GdkEventButton * event){
	switch(event->button){
		case 1:
			switch(window->context->graphMode){
				case GR_GRAPH:
				break;
				case GR_HISTOGRAM:
					gen_tres_texture(window->parent);
				break;
			}
		break;
	}
	return FALSE;
}

/*
 * Functions of menus processing
 */
void show_hist(Window *window){ // 'h' in main window
	window->context->graphMode = GR_HISTOGRAM;
	show_histogram(window);
}
void spots_coords(Window *window){ // 'c'
	get_spots_coords(window, window->image_transformed);
}

void set_spots_parameters(Window *window){ // 'alt+s'
	GtkBuilder *builder;
	GError *err = NULL;
	GtkWidget *dialog, *btn;
	void getval(GtkSpinButton *spin_button, int *val){
		*val = gtk_spin_button_get_value_as_int(spin_button);
		DBG("val = %d", *val);
	}
	void chk_tres_vals(GtkWidget *wnd){
		GtkMessageDialog *dialog;
		gchar *text;
		gboolean W, H;
		H = (Global->maxSpotH > Global->minSpotH);
		W = (Global->maxSpotW > Global->minSpotW);
		if(H && W){ // OK
			gtk_widget_destroy(wnd);
			return;
		}
		if(!(H || W))
			text = _("Both max values are less than min");
		else{
			if(!H)
				text = _("Max height value less than min");
			else
				text = _("Max width value less than min");
		}
		dialog = GTK_MESSAGE_DIALOG(gtk_message_dialog_new(NULL,
									GTK_DIALOG_MODAL,
									GTK_MESSAGE_ERROR,
									GTK_BUTTONS_CLOSE,
									text));
		gtk_window_set_transient_for(GTK_WINDOW(dialog), GTK_WINDOW(wnd));
		gtk_dialog_run(GTK_DIALOG(dialog));
		gtk_widget_destroy(GTK_WIDGET(dialog));
	}
	void conn_sig(gchar *object, int *val){
		GtkAdjustment *adj;
		GtkWidget *spin = GTK_WIDGET(gtk_builder_get_object(builder, object));
		g_signal_connect(spin, "value-changed", G_CALLBACK(getval), val);
		adj = (GtkAdjustment*) gtk_adjustment_new(*val, 1.0, 4000.0, 1.0, 0.0, 0.0);
		gtk_spin_button_set_adjustment(GTK_SPIN_BUTTON(spin), adj);
	}
	builder = gtk_builder_new();
	if(!gtk_builder_add_from_string(builder, ui, -1, &err)){
		g_warning("%d, %s\n", err->code, err->message);
		exit(-1);
	}
	dialog = GTK_WIDGET(gtk_builder_get_object(builder, "TresWindow"));
	conn_sig("MinWidth",  &Global->minSpotW);
	conn_sig("MaxWidth",  &Global->maxSpotW);
	conn_sig("MinHeight", &Global->minSpotH);
	conn_sig("MaxHeight", &Global->maxSpotH);
	btn = GTK_WIDGET(gtk_builder_get_object(builder, "TresCancel"));
	g_signal_connect_swapped(btn, "clicked",
						G_CALLBACK(gtk_widget_destroy), dialog);
	btn = GTK_WIDGET(gtk_builder_get_object(builder, "TresOK"));
	g_signal_connect_swapped(btn, "clicked",
						G_CALLBACK(chk_tres_vals), dialog);
	run_modal_window(GTK_WINDOW(dialog), window);
}

void circle_params(Window *window){ // 'x'
	get_circles_params(window, window->image_transformed);
}
void spot_ident(Window *window){ // 's'
	window->context->visualMode |= SHOW_POTSELECT;
	window->context->drawingMode = DRAW_SELECTION;
	force_redraw(window->drawingArea);
}
void save_spots(Window *window){
	if(!window->spots){
		g_err(_("No recognized spots"));
		return;
	}
	GtkFileChooser *dialog;
	gchar *filename = NULL, *tmp = NULL;
	gint response;
	dialog = GTK_FILE_CHOOSER(
				gtk_file_chooser_dialog_new("Save File", NULL,
						GTK_FILE_CHOOSER_ACTION_SAVE,
						GTK_STOCK_CANCEL, GTK_RESPONSE_CANCEL,
						GTK_STOCK_SAVE, GTK_RESPONSE_ACCEPT,
						NULL));
	if(saved_path)
		gtk_file_chooser_set_current_folder(dialog, saved_path);
	filename = g_strdup_printf("%s", window->image->filename);
	DBG("saved path: %s; filename: %s", saved_path, filename);
	if( (tmp = strrchr(filename, '.')) ) *tmp = 0;
	tmp = g_strdup_printf("%s.spotlist", filename);
	gtk_file_chooser_set_current_name(dialog, tmp);
	g_free(tmp);
	g_free(filename); filename = NULL;
	gtk_file_chooser_set_do_overwrite_confirmation(dialog, TRUE);
	response = run_modal_dialog(GTK_DIALOG(dialog), window);
	if(response == GTK_RESPONSE_ACCEPT){
		filename = gtk_file_chooser_get_filename(dialog);
	}
	gtk_widget_destroy(GTK_WIDGET(dialog));
	if(!filename) return;
	spots_save(window->spots, filename);
	g_print("%s saved", filename);
	tmp = g_path_get_basename(filename);
	g_free(filename);
	filename = g_strdup_printf("%s saved", tmp);
	set_status_text(StatusText, filename, window);
	g_free(filename); g_free(tmp);
}
void draw_tracks(Window *window){ // 't'
	window->context->drawingMode = DRAW_TRACK;
	window->context->visualMode &= ~SHOW_POTSELECT;
}
void zoom_restore(Window *window){ // '0'
	window->Zoom = 1.;
	window->move.x = 0.;
	window->move.y = 0.;
	window->Zangle = 0.;
	window->Xangle = 0.;
	force_redraw(window->drawingArea);
}
void zooom_qframe(Window *window){ // 'z'
	window->context->drawingMode = DRAW_QUAD;
	zoom_frame = TRUE;
}
void select_terrain(Window *window){ // '3'
	window->context->drawingMode = DRAW_QUAD;
	select_3D = TRUE;
}
void filter(Window *window){ // 'f'
/*	Filter *f;
	INIT(f, Filter);*/
/*	f->FilterType = MEDIAN;
	f->sx=f->sy=9.;
	f->w = f->h = 3;
	filter_image(window, f);*/
/*	f->FilterType = SIMPLEGRAD;
	filter_image(window, f);
*/
	find_contours(window, 16, LOG);
/*	f->FilterType = STEP;
	f->w = 16;
	f->h = LOG;
	filter_image(window, f);*/
	/*f->FilterType = LAPGAUSS;
	f->w = f->h = 256;
	f->sx = f->sy = 5.f;
	filter_image(window, f);*/
}

/*
 * Processing of a key pressing events
 * 1. Main window
 */
gboolean press_main_key(Window *window, GdkEventKey * event){
	//FNAME();
	gboolean ret = TRUE;
	if(event->keyval == GDK_Escape){
		if(zoom_frame || select_3D){
			window->context->drawingMode = DRAW_TRACK;
			zoom_frame = FALSE; // ESC cansels zoom
			select_3D = FALSE;
		}
		return TRUE;
	}
	if(event->state & (GDK_CONTROL_MASK)){ // with ctrl
		switch(event->keyval){
			default: ret = FALSE;
		}
	}
	else{
		switch(event->keyval){
			case GDK_Page_Up: case GDK_Up: case GDK_Left:
				if(window->files.list_length > 1)
					change_image(get_filename(PREVIOUS, window), window);
			break;
			case GDK_Page_Down: case GDK_Down: case GDK_Right:
				if(window->files.list_length > 1)
					change_image(get_filename(NEXT, window), window);
			break;
			case GDK_Home:
				if(window->files.list_length > 1)
					change_image(get_filename(FIRST, window), window);
			break;
			case GDK_End:
				if(window->files.list_length > 1)
					change_image(get_filename(LAST, window), window);
			break;
			default: ret = FALSE;
		}
	}
	//refresh_state(window);
	return ret;
}
// 2. Auxiliary  window
/*
gboolean press_sub_key(Window *window, GdkEventKey * event){
	//FNAME();
	gboolean ret = TRUE;
	switch(event->keyval){
		default: ret = FALSE;
	}
	return ret;
}
*/

/*
gboolean release_key(GtkWidget * none, GdkEventKey * event){
	//DBG("Release key %d\n", event->keyval);
	return FALSE;
}*/

// Scroll
gboolean main_scroll(Window *window, GdkEventScroll *event){
	gboolean ret = FALSE;
	if(event->state & (GDK_CONTROL_MASK)){
		switch(event->direction){
			case GDK_SCROLL_UP:
				zoom_far(window);
				ret = TRUE;
			break;
			case GDK_SCROLL_DOWN:
				zoom_near(window);
				ret = TRUE;
			break;
			default: break;
		}
	}
	return ret;
}

/*
 * Mouse move events
 * 1. Main window
 */
gboolean main_motion(Window *window, GdkEventMotion * event){
	double x,y;
	gchar buf[40];
	conv_mouse_to_image_coords(event->x, event->y, &x, &y, window);
	g_snprintf(buf, 39, "(%.1f, %.1f)", x, y);
	set_status_text(StatusCoords, buf, window);
	if(event->state & GDK_BUTTON1_MASK){
		switch(window->context->drawingMode){
		case DRAW_TRACK: // In tracking mode show track
		case DRAW_QUAD:  // the same in frame selection mode
			if(is_tracking) do_tracking(event->x, event->y, FALSE, window);
		break;
		}
	}
	else if(event->state & GDK_BUTTON2_MASK){
		window->move.x = lX0 + (event->x - X0) * window->Daspect;
		window->move.y = lY0 + (Y0 - event->y) * window->Daspect;
		force_redraw(window->drawingArea);
	}
	else if(event->state & GDK_BUTTON3_MASK){
	}
	return TRUE;
}

// 2. Auxiliary  window
gboolean sub_motion(Window *window, GdkEventMotion * event){
	double x = event->x / xScale;
	double y = gYmax - event->y/yScale;
	double xx, yy;
	gchar buf[40];
	get_image_crds(x, &xx, &yy);
	if(isYlog(window)) y = exp(y);
	g_snprintf(buf, 39, "(%.3f, %.3f) I(%.1f, %.1f)", x, y, xx, yy);
	set_status_text(StatusCoords, buf, window);
	if(event->state & GDK_BUTTON1_MASK){
		if(is_tracking) graph_mouse_btn((int)event->x, window);
	}
	return TRUE;
}

// 3. OpenGL  window
void open3D(Window *window){
	terrain_3D(window, 0.,0.,0.,0., TRUE); // Show open dialog & load new picture
}
void rotXp(Window *window){
	window->Xangle += key_acceleration;
	if(window->Xangle > 360.) window->Xangle -= 360.;
	force_redraw(window->drawingArea);
}
void rotXm(Window *window){
	window->Xangle -= key_acceleration;
	if(window->Xangle < 0.) window->Xangle += 360.;
	force_redraw(window->drawingArea);
}
void rotZp(Window *window){
	window->Zangle += key_acceleration;
	if(window->Zangle > 360.) window->Zangle -= 360.;
	force_redraw(window->drawingArea);
}
void rotZm(Window *window){
	window->Zangle -= key_acceleration;
	if(window->Zangle < 0.) window->Zangle += 360.;
	force_redraw(window->drawingArea);
}
void moveXp(Window *window){
	window->move.x += key_acceleration;
	force_redraw(window->drawingArea);
}
void moveXm(Window *window){
	window->move.x -= key_acceleration;
	force_redraw(window->drawingArea);
}
void moveYp(Window *window){
	window->move.y += key_acceleration;
	force_redraw(window->drawingArea);
}
void moveYm(Window *window){
	window->move.y -= key_acceleration;
	force_redraw(window->drawingArea);
}
void game_move(Window *window, double where){
	double lx, ly, lz, za = window->Zangle, xa = window->Xangle;
	lx = cos(za)*sin(xa);
	ly = -sin(za);
	lz = -cos(za)*cos(xa);
	double delta = where * key_acceleration;
	window->move.x += delta*lx;
	window->move.y += delta*ly;
	window->Zoom += delta*lz/window->image->height;
}
void moveZp(Window *window){
	if(window->context->drawingMode == GAME_MODE)
		game_move(window, 1.);
	else
		window->Zoom *= 1.1;
	force_redraw(window->drawingArea);
}
void moveZm(Window *window){
	if(window->context->drawingMode == GAME_MODE)
		game_move(window, -1.);
	else
		window->Zoom /= 1.1;
	force_redraw(window->drawingArea);
}
void game_mode(Window *window){
	if(window->context->drawingMode == GAME_MODE){
		window->context->drawingMode = DRAW_TRACK;
		gdk_pointer_ungrab(GDK_CURRENT_TIME);
		gtk_window_unfullscreen(GTK_WINDOW(window->window));
		gdk_window_set_cursor(GDK_WINDOW(window->drawingArea->window), NULL);
	}
	else{
		window->context->drawingMode = GAME_MODE;
		gdk_pointer_grab(window->drawingArea->window, TRUE, 0,
						window->drawingArea->window, NULL,
						GDK_CURRENT_TIME);
		gtk_window_fullscreen(GTK_WINDOW(window->window));
		GdkCursor *cursor = gdk_cursor_new(GDK_BLANK_CURSOR);
		gdk_window_set_cursor(GDK_WINDOW(window->drawingArea->window), cursor);
	}
	zoom_restore(window);
}
gfloat oldXbtn, oldYbtn;
gint oldXscreen, oldYscreen;
gboolean CursorWasMoved = TRUE;
gboolean press_ogl_mbtn(Window *window, GdkEventButton *event){
	if(window->context->drawingMode == GAME_MODE &&
			event->button == 1){
		GdkDisplay *disp = gdk_display_get_default();
		GdkScreen* screen = gdk_display_get_default_screen(disp);
		oldXscreen = gdk_screen_get_width(screen) / 2;
		oldYscreen = gdk_screen_get_height(screen) / 2;
		gdk_display_warp_pointer(disp, screen,
				oldXscreen, oldYscreen);
		CursorWasMoved = TRUE;
	}
	return FALSE;
}
gboolean release_ogl_mbtn(Window *window  __attribute__((unused)),
						GdkEventButton *event  __attribute__((unused))){
	return FALSE;
}
/*
 * Move mouse in 3D mode
 * !!! In game mode angles measures in radians, not degrees !!!
 */
gboolean ogl_motion(Window *window, GdkEventMotion * event){
	if(window->context->drawingMode == GAME_MODE
			&& event->state & GDK_BUTTON1_MASK){
		if(CursorWasMoved){
			CursorWasMoved = FALSE;
			oldXbtn = event->x;
			oldYbtn = event->y;
			return FALSE;
		}
		window->Xangle += (event->x - oldXbtn)/1000.;
		window->Zangle += (event->y - oldYbtn)/1000.;
		GdkDisplay *disp = gdk_display_get_default();
		gdk_display_warp_pointer(disp, gdk_display_get_default_screen(disp),
						oldXscreen, oldYscreen);
		CursorWasMoved = TRUE;
		force_redraw(window->drawingArea);
	}
	return FALSE;
}
gboolean reset_accel(){
	key_acceleration = 1.;
	return FALSE;
}
gboolean press_ogl_key(Window *window, GdkEventKey * event){
	static gint tag = 0;
	if(tag) gtk_timeout_remove(tag);
	tag = gtk_timeout_add(100, reset_accel, NULL);
	key_acceleration += 1.;
	switch(event->keyval){
		case GDK_Left:
			moveXp(window); // Move left - picture moves to the right
		break;
		case GDK_Right:
			moveXm(window);
		break;
		case GDK_Up:
			moveZp(window);
		break;
		case GDK_Down:
			moveZm(window);
		break;
		default:
			return FALSE;
	}
	return TRUE;
}
gboolean release_ogl_key(Window *window  __attribute__((unused)),
						GdkEventKey *event  __attribute__((unused))){
	return FALSE;
}
gboolean ogl_scroll(Window *window  __attribute__((unused)),
					GdkEventScroll *event  __attribute__((unused))){
	return FALSE;
}

/*
 * Redefine rulers limits, x0 & xm are the limits of horizontal ruler
 * y0 & ym - vertical
 * hrule, vrule - appropriate rulers
 */
void set_rules(double x0, double y0, double xm, double ym,
				GtkWidget *hrule, GtkWidget *vrule){
	double yy = (y0 < 0.) ? 0 : y0;
	gtk_ruler_set_range(GTK_RULER (hrule), x0, xm, 0., xm);
	gtk_ruler_set_range(GTK_RULER (vrule), y0, ym, yy, ym);
}
void set_Drulers(double x0, double y0, double xm, double ym, Window *window){
	set_rules(-x0, ym, xm, -y0, window->hRule, window->vRule);
}
void set_Grulers(double y0, double xm, double ym, Window *window){
	//DBG("gXm=%g, gY0=%g, gYm=%g\n", xm, y0, ym);
	set_rules(0., ym, xm, y0, window->hRule, window->vRule);
}

void destroy_window(Window* window){
	FNAME();
	if(!window) return;
	if(window->id == GRAPH_WINDOW)
		window->parent->graphWindow = NULL;
	else{
		destroy_window(window->graphWindow);
		_FREE(window->context);
	}
	if(window->id == GL3D_WINDOW)
		freeGLmemory(window);
	if(window->image){
		if(window->id == GL3D_WINDOW){
			window->image->cLevels = NULL;
		}
		destroy_image(window->image);
	}
	if(window->image_transformed)
		destroy_image(window->image_transformed);
	_FREE(window->texture);
	spots_free(&window->spots);
	free_filelist(window);
	gtk_widget_destroy(window->window);
	free(window);
}

/*
 * Signal handler structure definition
 * array Sighandler[] should end by NULL element
 */
typedef struct{
	gchar *elementName;		// name of element
	gchar *sigName;			// signal name
	GCallback sigHandler;	// function - handler, a kind of void f(Window* w)
} Sighandler;
// Signals for graph window
Sighandler graphSigHandler[] = {
	{"graphCloseMenuItem", "activate", G_CALLBACK(close_subwindow)},
	{"graphGaussMenuItem", "activate", G_CALLBACK(fit_gaussian)},
	{"graphLinMenuItem", "activate", G_CALLBACK(lin_scale)},
	{"graphLogMenuItem", "activate", G_CALLBACK(log_scale)},
	{"graphArea", "expose-event", G_CALLBACK(Gexpose)},
	{"graphArea", "configure-event", G_CALLBACK(Gconfigure)},
	{NULL, "button-press-event", G_CALLBACK(press_sub_mbtn)},
	{NULL, "button-release-event", G_CALLBACK(release_sub_mbtn)},
	{NULL, "motion-notify-event", G_CALLBACK(sub_motion)},
//	{NULL, "key-press-event", G_CALLBACK(press_sub_key)},
	{NULL, NULL, NULL}
};
// Signals for window with picture
Sighandler mainSigHandler[] = {
	{"mainOpenMenuItem", "activate", G_CALLBACK(open_file)},
	{"mainSaveAsMenuItem", "activate", G_CALLBACK(save_as)},
	{"mainOpenInNewMenuItem", "activate", G_CALLBACK(open_in_new)},
	{"mainOpenIn3DMenuItem", "activate", G_CALLBACK(open_in_3D)},
	{"mainHistMenuItem", "activate", G_CALLBACK(show_hist)},
	{"mainHeadersMenuItem", "activate", G_CALLBACK(edit_headers)},
	{"mainSpotsidentMenuItem", "activate", G_CALLBACK(spot_ident)},
	{"mainSaveSpotsMenuItem", "activate", G_CALLBACK(save_spots)},
	{"mainTrackMenuItem", "activate", G_CALLBACK(draw_tracks)},
	{"main3DviewFrameMenuItem", "activate", G_CALLBACK(select_terrain)},
	{"main3DviewFullMenuItem", "activate", G_CALLBACK(full_3D_terrain)},
	{"mainZoomframeMenuItem", "activate", G_CALLBACK(zooom_qframe)},
	{"mainZoomrestoreMenuItem", "activate", G_CALLBACK(zoom_restore)},
	{"mainZoomnearMenuItem", "activate", G_CALLBACK(zoom_near)},
	{"mainZoomfarMenuItem", "activate", G_CALLBACK(zoom_far)},
	{"mainSpotsParamMenuItem", "activate", G_CALLBACK(set_spots_parameters)},
	{"mainFindSpotsMenuItem", "activate", G_CALLBACK(spots_coords)},
	{"mainHartmannMenuItem", "activate", G_CALLBACK(hartmann_spots)},
	{"mainCirclesMenuItem", "activate", G_CALLBACK(circle_params)},
	{"mainHoughMenuItem", "activate", G_CALLBACK(hough_lines)},
	{"mainFilterMenuItem", "activate", G_CALLBACK(filter)},
//	{"", "activate", G_CALLBACK()},
	{NULL, "button-press-event", G_CALLBACK(press_main_mbtn)},
	{NULL, "button-release-event", G_CALLBACK(release_main_mbtn)},
	{NULL, "motion-notify-event", G_CALLBACK(main_motion)},
	{NULL, "key-press-event", G_CALLBACK(press_main_key)},
	{NULL, "scroll-event", G_CALLBACK(main_scroll)},
	{NULL, NULL, NULL}
};
// Signals for 3D window
Sighandler oglSigHandler[] = {
	{"oglOpenMenuItem", "activate", G_CALLBACK(open3D)},
	{"oglCloseMenuItem", "activate", G_CALLBACK(destroy_window)},
	{"oglZoomrestoreMenuItem", "activate", G_CALLBACK(zoom_restore)},
	{"oglXAplusMenuItem", "activate", G_CALLBACK(rotXp)},
	{"oglXAminusMenuItem", "activate", G_CALLBACK(rotXm)},
	{"oglZAplusMenuItem", "activate", G_CALLBACK(rotZp)},
	{"oglZAminusMenuItem", "activate", G_CALLBACK(rotZm)},
	{"oglXplusMenuItem", "activate", G_CALLBACK(moveXp)},
	{"oglXminusMenuItem", "activate", G_CALLBACK(moveXm)},
	{"oglYplusMenuItem", "activate", G_CALLBACK(moveYp)},
	{"oglYminusMenuItem", "activate", G_CALLBACK(moveYm)},
	{"oglZplusMenuItem", "activate", G_CALLBACK(moveZp)},
	{"oglZminusMenuItem", "activate", G_CALLBACK(moveZm)},
	{"oglGameModeMenuItem", "toggled", G_CALLBACK(game_mode)},
	{NULL, "button-press-event", G_CALLBACK(press_ogl_mbtn)},
	{NULL, "button-release-event", G_CALLBACK(release_ogl_mbtn)},
	{NULL, "motion-notify-event", G_CALLBACK(ogl_motion)},
	{NULL, "key-press-event", G_CALLBACK(press_ogl_key)},
	{NULL, "key-release-event", G_CALLBACK(release_ogl_key)},
	{NULL, "scroll-event", G_CALLBACK(ogl_scroll)},
	{NULL, NULL, NULL}
};

inline GtkWidget *initA(gchar *prefix, GtkBuilder *builder){
	gchar text[32];
	g_sprintf(text, "%sArea", prefix);
	GtkWidget *A = GTK_WIDGET(gtk_builder_get_object(builder, text));
	gtk_widget_set_events(A,	GDK_EXPOSURE_MASK		|
								GDK_BUTTON_PRESS_MASK	|
								GDK_POINTER_MOTION_MASK	|
								GDK_POINTER_MOTION_HINT_MASK|
								GDK_BUTTON_RELEASE_MASK);
	return A;
}
inline GtkWidget *initR(GtkWidget *area, gchar *name, GtkBuilder *builder){
	GtkWidget *R = GTK_WIDGET(gtk_builder_get_object(builder, name));
	g_signal_connect_swapped(area, "motion-notify-event",
			G_CALLBACK(EVENT_METHOD(R, motion_notify_event)), R);
	return R;
}

/*
 * Window initialisation
 * parent - parent window
 * handlers - array of signal handlers
 * prefix - window name prefix (main, graph, ...)
 * nBlocks - amount of blocks in status bar
 */
Window *init_window(Window *parent, int winId){
	FNAME();
	guint i;
	GtkBuilder *builder;
	Window *window;
	Sighandler *handlers = NULL;
	gchar *prefix = NULL, *title = NULL;
	guint nBlocks = 0;
	GError *err = NULL;
	gchar text[32];
	builder = gtk_builder_new();
	if(!gtk_builder_add_from_string(builder, ui, -1, &err)){
		g_warning("%d, %s\n", err->code, err->message);
		exit(-1);
	}
	INIT(window, Window);
	window->id = winId;
	if(parent){
		window->parent = parent;
		if(winId == GRAPH_WINDOW) // if window is graph context equal to parent's
			window->context = parent->context;
		else{
			window->context = COPY(parent->context, Context);
			window->context->visualMode = SHOW_ONLY_IMAGE;
		}
	}else
		INIT(window->context, Context);
	switch(winId){
		case MAIN_WINDOW:
		case OPENGL_WINDOW:
			prefix = "main";
			handlers = mainSigHandler;
			nBlocks = 4;
			break;
		case GRAPH_WINDOW:
			prefix = "graph";
			handlers = graphSigHandler;
			nBlocks = 3;
			break;
		case GL3D_WINDOW:
			prefix = "ogl";
			handlers = oglSigHandler;
			break;
	}
	g_sprintf(text, "%sWindow", prefix);
	window->window = GTK_WIDGET(gtk_builder_get_object(builder, text));
	if(winId == MAIN_WINDOW || winId == OPENGL_WINDOW){
		GtkAction *quit = GTK_ACTION(gtk_builder_get_object(builder, "mainQuitMenuItem"));
		gtk_action_disconnect_accelerator(quit);
		if(winId == MAIN_WINDOW){
			g_signal_connect(quit, "activate", chk_quit, NULL);
			gtk_action_set_icon_name(quit, "application-exit");
			gtk_action_set_accel_path(quit, "<Main>/quit");
			gtk_action_set_label(quit, _("Quit"));
		}else{
			gtk_action_set_accel_path(quit, "<Main>/close");
			gtk_action_set_label(quit, _("Close"));
			g_signal_connect_swapped(quit, "activate",
				G_CALLBACK(destroy_window), window);
			gtk_action_set_icon_name (quit, "window-close");
		}
		gtk_action_connect_accelerator(quit);
	}
	window->drawingArea = initA(prefix, builder);
	if(winId != GL3D_WINDOW){
		g_sprintf(text, "%sHruler", prefix);
		window->hRule = initR(window->drawingArea, text, builder);
		g_sprintf(text, "%sVruler", prefix);
		window->vRule = initR(window->drawingArea, text, builder);
	}
	for(i = 0; i < nBlocks; i++){
		g_sprintf(text, "%sStatusBar%d", prefix, i);
		window->SBars[i] = GTK_ENTRY(gtk_builder_get_object(builder, text));
		#ifdef EBUG
		set_status_text(i, text, window);
		#endif
	}
	window->statusBlocks = nBlocks;
	if(winId != GRAPH_WINDOW){
		INIT(window->image, IMAGE);
		INIT(window->image_transformed, IMAGE);
		initGl(window->drawingArea);
	}else{
		window->LinLogMenu[0] = GTK_RADIO_ACTION(gtk_builder_get_object(builder, "graphLinMenuItem"));
		window->LinLogMenu[1] = GTK_RADIO_ACTION(gtk_builder_get_object(builder, "graphLogMenuItem"));
		gtk_toggle_action_set_active(GTK_TOGGLE_ACTION(window->LinLogMenu[0]), TRUE);
	}
	while(handlers->sigHandler){
		void *element;
		if(handlers->elementName)
			element = (void*)gtk_builder_get_object(builder, handlers->elementName);
		else
			element = (void*)window->window;
		g_signal_connect_swapped(element, handlers->sigName,
				handlers->sigHandler, window);
		handlers++;
	}
	gtk_builder_connect_signals(builder, NULL);
	g_object_unref(G_OBJECT(builder));
	if(winId != GRAPH_WINDOW){
		if(winId == MAIN_WINDOW || winId == OPENGL_WINDOW)
			window->graphWindow = init_window(window, GRAPH_WINDOW);
		gtk_widget_show_all(window->window);
	}
	window->Daspect = 1.; window->Zoom = 1.;
	if(winId != OPENGL_WINDOW)
		title = g_strdup_printf("%s (%s)", PRGNAME, prefix);
	else
		title = g_strdup_printf("%s", PRGNAME);
	gtk_window_set_title(GTK_WINDOW(window->window), title);
	g_free(title);
	return window;
}

void init_main_window(int *argc, char ***argv){
	gtk_init(argc, argv);
	gtk_gl_init(argc, argv);
	gtk_accel_map_add_entry("<Main>/close", GDK_w, GDK_CONTROL_MASK);
	gtk_accel_map_add_entry("<Main>/quit", GDK_q, GDK_CONTROL_MASK);
	mainWindow = init_window(NULL, MAIN_WINDOW);
}

/*
 * *prefocal = TRUE if image is prefocal, FALSE otherwisw
 * (Ask user, this function calls only if there's a problems to find VAL_F in FITS keys)
 */
void get_prefocal(Window *window, gboolean *prefocal){
	double foc;
	GtkWidget *dialog, *btn, *focus, *radio;
	void ch_val_f(GtkEntry *e){
		char *eptr = NULL;
		errno = 0;
		foc = strtod(gtk_entry_get_text(e), &eptr);
		if((eptr && *eptr) || errno == ERANGE){
			gchar *txt = g_strdup_printf("%g", foc);
			gtk_entry_set_text(e, txt);
			g_free(txt);
		}
		DBG("F: %g", foc);
	}
	void ch_vals(GtkWidget *w){
		*prefocal = !gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(radio));
		DBG("prefocal: %d\n", *prefocal);
		window->image->val_f = foc;
		gtk_widget_destroy(w);
	}
	GtkBuilder *builder;
	GError *err = NULL;
	gchar *text;
	if(!window || !window->image) return;
	if(window->image->val_f > 0.) text = g_strdup_printf("%.3f", window->image->val_f);
	else text = g_strdup_printf("0.000");
	builder = gtk_builder_new();
	if(!gtk_builder_add_from_string(builder, ui, -1, &err)){
		g_warning("%d, %s\n", err->code, err->message);
		exit(-1);
	}
	dialog = GTK_WIDGET(gtk_builder_get_object(builder, "FocalWindow"));
	radio = GTK_WIDGET(gtk_builder_get_object(builder, "FocalPostRadioBtn"));
	focus = GTK_WIDGET(gtk_builder_get_object(builder, "FocFocusEntry"));
	g_signal_connect(focus, "changed", G_CALLBACK(ch_val_f), NULL);
	gtk_entry_set_text(GTK_ENTRY(focus), text);
	btn = GTK_WIDGET(gtk_builder_get_object(builder, "FocalOK"));
	g_signal_connect_swapped(btn, "clicked",
						G_CALLBACK(ch_vals), dialog);
	run_modal_dialog(GTK_DIALOG(dialog), window);
	g_free(text);
}
