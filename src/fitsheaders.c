//      fitsheaders.c - functions to edit FITS head
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
#include "fitsheaders.h"
#include "fits.h"

GtkListStore *combo_ = NULL;
GtkTreeSelection *sel;
GtkTreeView* Tree;
GtkWidget *del_button;
GtkWidget *window;
// Data types to write FITS keys into menu
gchar* key_types[] = {
	"TSTRING", "TLOGICAL", "TBYTE", "TSHORT", "TUSHORT", "TINT", "TUINT",
	"TLONG", "TULONG", "TLONGLONG", "TFLOAT", "TDOUBLE", "TCOMPLEX",
	"TDBLCOMPLEX"
};
// Values itself
int type_vals[] = {
	TSTRING,
	TLOGICAL,
	TBYTE,
	TSHORT,
	TUSHORT,
	TINT,
	TUINT,
	TLONG,
	TULONG,
	TLONGLONG,
	TFLOAT,
	TDOUBLE,
	TCOMPLEX,
	TDBLCOMPLEX
};

int VALS_N; // Number of items

gchar *colnames[] = {N_("Type"), N_("Name"), N_("Value"), N_("Comment")};
#define VISIBLE_COLS 4 // The amount of visible columns

int find_list_type(char *type){
	int i;
	for(i = 0; i < VALS_N; i++)
		if(strcmp(key_types[i], type) == 0) break;
	if(i == VALS_N) // Error: there's no such value in list
		i = 0;
	return i;
}

gboolean ll_conv(gchar **value, int type){
	long long ll_val, oll_val;
	gboolean ret = TRUE;
	ll_val = oll_val = strtoll(*value, NULL, 0);
	//g_free(*value);
	switch(type){
		case TBYTE:
			if(ll_val > CHAR_MAX || ll_val < CHAR_MIN) ll_val = 0L;
		break;
		case TSHORT:
			if(ll_val > SHRT_MAX || ll_val < SHRT_MIN) ll_val = 0L;
		break;
		case TINT:
			if(ll_val > INT_MAX || ll_val < INT_MIN) ll_val = 0L;
		break;
		case TLONG:
			if(ll_val > LONG_MAX || ll_val < LONG_MIN) ll_val = 0L;
		break;
		case TLONGLONG:
			if(ll_val > LLONG_MAX || ll_val < LLONG_MIN) ll_val = 0L;
		break;
		case TUSHORT:
			if(ll_val > USHRT_MAX || ll_val < 0L) ll_val = 0L;
		break;
		case TUINT:
			if(ll_val > UINT_MAX || ll_val < 0L) ll_val = 0L;
		break;
		case TULONG:
			if((unsigned long)ll_val > ULONG_MAX || ll_val < 0L) ll_val = 0L;
		break;
		default:
			ll_val = 0L;
		break;
	}
	DBG("long long: %lld", ll_val);
	if(ll_val != oll_val){
		g_err(_("Error in value range!"));
		ret = FALSE;
	}
	g_free(*value);
	*value = g_strdup_printf("%lld", ll_val);
	return ret;
}

gboolean mystrtod(gchar **val){
	double d_val;
	gboolean ret = TRUE;
	char *eptr = NULL, *ptr = *val;
	errno = 0;
	d_val = strtod(ptr, &eptr);
	DBG("double: %.16e (was: %s)", d_val, *val);
	if((eptr && *eptr) || errno == ERANGE){
		g_free(*val);
		DBG("try convert %g into val", d_val);
		*val = g_strdup_printf("%.16e", d_val);
		g_err(_("Invalid double number!"));
		ret = FALSE;
	}
	DBG("val: %s", *val);
	return ret;
}

gboolean check_complex_number(gchar **val){
	double re = 0., im = 0.;
	gboolean ret = TRUE;
	if(!cmplx_conv(*val, &re, &im)){
		gchar *msg = g_strdup_printf("%s\n%s",
			_("Format error: complex number must be in format"),
			"RE + iIM, RE + jIM, RE iIM, RE jIM, (RE, IM) or (RE IM)");
		g_err(msg);
		g_free(msg);
		ret = FALSE;
	}
	DBG("complex: (%g, %g)", re, im);
	g_free(*val);
	*val = g_strdup_printf("(%g, %g)", re, im);
	return ret;
}

gboolean check_val(gchar **newval, int type){
	gboolean b_val, ret = TRUE;
	switch(type){
		case TBYTE:
		case TSHORT:
		case TINT:
		case TLONG:
		case TLONGLONG:
		case TUSHORT:
		case TUINT:
		case TULONG:
		// convert *char -> long long vith checking value
			ret = ll_conv(newval, type);
		break;
		case TFLOAT:
		case TDOUBLE:
		// convert *char -> double
			ret = mystrtod(newval);
		break;
		case TCOMPLEX:
		case TDBLCOMPLEX:
		// check whether a number is complex
			ret = check_complex_number(newval);
		break;
		case TLOGICAL:
		// convert value into boolean
			if(g_ascii_strncasecmp("true", *newval, 4)==0) b_val = TRUE;
			else
				b_val = atoi(*newval);
			g_free(*newval);
			*newval = g_strdup_printf("%s", b_val ? "TRUE" : "FALSE");
			DBG("boolean: %s", *newval);
		break;
		default:
			DBG("string: %s", *newval); // all other values a strings
	}
	return ret;
}

gboolean try_change_val(IMAGE *ima, gchar **newval){
	int type;
	GtkTreeIter cur_items;
	gboolean ret;
	GtkTreeModel *model;
	FNAME();
	if(!gtk_tree_selection_get_selected(sel,&model, &cur_items)) return FALSE;
	gtk_tree_model_get(GTK_TREE_MODEL(ima->store), &cur_items, L_TYPE, &type, -1);
	ret = check_val(newval, type);
	gtk_list_store_set(ima->store, &cur_items, L_VAL, *newval, -1);
	DBG("ret");
	return ret;
}
guint cur_column;
void edit_cell(	GtkCellRendererText *cell,
				gchar *path,
				gchar *new_text,
				IMAGE *ima){
	GtkTreeModel *model;
	GtkTreeIter cur_items;
	cur_column =  GPOINTER_TO_UINT(g_object_get_data(G_OBJECT(cell), "column_num"));
	gboolean nxtcol = TRUE, editable = TRUE;
	if(!path) return;
	DBG("path: %s, new: %s, colnum: %d", path, new_text, cur_column);
	gchar *textptr;
	if(!gtk_tree_selection_get_selected(sel,&model, &cur_items)) return;
	switch(cur_column){
		case L_TYPENM: // key type was changed
			gtk_list_store_set(ima->store, &cur_items, L_TYPE,
						type_vals[find_list_type(new_text)], -1);
			gtk_list_store_set(ima->store, &cur_items, cur_column, new_text, -1);
		break;
		case L_VAL: // key value was changed
			textptr = g_strdup_printf("%s", new_text);
			nxtcol = try_change_val(ima, &textptr); // try to change parameter according to its type
			g_free(textptr);
		break;
		default:
			gtk_list_store_set(ima->store, &cur_items, cur_column, new_text, -1);
		break;
	}
	if(nxtcol){
		if(cur_column < L_COMM) // it's not a last column - select next
			cur_column++;
		else{
			cur_column = 0;
			editable = FALSE;
		}
	}
	DBG("set cursor col:%d, e=%d", cur_column, editable);
	gtk_tree_view_set_cursor(Tree,
			gtk_tree_model_get_path(GTK_TREE_MODEL(ima->store), &cur_items),
			gtk_tree_view_get_column(Tree, cur_column), editable);
	DBG("ret");
}
/*
void cancel_edit(GtkCellRenderer *renderer){
	cur_column =  GPOINTER_TO_UINT(g_object_get_data(G_OBJECT(renderer), "column_num"));
	gtk_tree_view_set_cursor(Tree,
				gtk_tree_model_get_path(GTK_TREE_MODEL(store), &cur_items),
				gtk_tree_view_get_column(Tree, cur_column), FALSE);
}
*/
void add_key(	IMAGE *ima,
				gchar *keytype, gchar *keyname, gchar *keyval,
				gchar *keycomment, int keyfitstype){
	GtkTreeIter iter;
	//~ FNAME();
	gtk_list_store_append(ima->store, &iter);
	gtk_list_store_set(ima->store, &iter,
			L_TYPENM,	keytype,
			L_KEY,		keyname,
			L_VAL,		keyval,
			L_COMM,		keycomment,
			L_TYPE,		keyfitstype,
			-1);// the last element must be  -1
}

void init_keylist(IMAGE *ima){
	/* Create a list model.
	 * We will have 4 columns: name of key type, key name, key value, comment, type
	 */
	if(ima->store) g_object_unref(G_OBJECT(ima->store)); // clear list, if it already exists
	ima->store = gtk_list_store_new(L_MAX, G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING, G_TYPE_INT);
}

GtkWidget *create_treeview(IMAGE *ima){
	int i, i_max = sizeof(key_types)/sizeof(char*);;
	GtkTreeIter iter;
	GtkCellRenderer *renderer;
	GtkTreeViewColumn *column;
	if(combo_)
		gtk_list_store_clear(combo_);
	combo_ = gtk_list_store_new(1, G_TYPE_STRING);
	for(i = 0; i < i_max; i++){
		gtk_list_store_append(combo_, &iter);
		gtk_list_store_set(combo_, &iter, 0, key_types[i], -1);
	}
	// create list (treeview) according model store
	Tree = (GtkTreeView*)gtk_tree_view_new_with_model(GTK_TREE_MODEL(ima->store));
	// Create columns in list
	for(i = 0; i < VISIBLE_COLS; i++){
		if(i == 0){
			renderer = gtk_cell_renderer_combo_new();
			g_object_set(renderer, "text-column", 0,
						"editable", TRUE, "has-entry", FALSE, "model",
						GTK_TREE_MODEL(combo_), NULL);
		}
		else{
			renderer = gtk_cell_renderer_text_new ();
			g_object_set(renderer, "editable", TRUE, NULL);
		}
		g_object_set_data(G_OBJECT(renderer), "column_num", GINT_TO_POINTER(i));
		column = gtk_tree_view_column_new_with_attributes(_(colnames[i]),
								renderer, "text", i, NULL);
		gtk_tree_view_append_column(Tree, column); // insert column
		gtk_tree_view_column_set_sort_column_id(column, i);// it may be sorted
		g_signal_connect(renderer, "edited", G_CALLBACK(edit_cell), ima);
	//	g_signal_connect(renderer, "editing-canceled", G_CALLBACK(cancel_edit), NULL);
	}
/*
	renderer = gtk_cell_renderer_text_new ();
	column = gtk_tree_view_column_new_with_attributes ("hidden",renderer, NULL);
	g_object_set(column, "visible", FALSE, NULL);
	gtk_tree_view_append_column(Tree, column);
*/
	return GTK_WIDGET(Tree);
}

static void list_changed(GtkTreeSelection *sel){
	GtkTreeModel *model;
	GtkTreeIter cur_items;
	if(gtk_tree_selection_get_selected(sel,&model, &cur_items))
		gtk_widget_set_sensitive(del_button, TRUE);
	else
		gtk_widget_set_sensitive(del_button, FALSE);
	/*
	GtkTreeModel *model;
	FNAME();
	gboolean showdel = FALSE;
	if(!GTK_IS_TREE_SELECTION(sel)){
		showdel = TRUE;
		if(GTK_IS_TREE_SELECTION(obj)) DBG("!!!");
		cur_column = GPOINTER_TO_UINT(g_object_get_data(G_OBJECT(sel), "column_num"));
		DBG("curcol: %d", cur_column);
	}
	else if(gtk_tree_selection_get_selected(sel,&model, &cur_items)){
		gint column = L_VAL;
		gchar *value;
		gtk_tree_model_get(model, &cur_items, column, &value, -1);
		DBG("selected: %s", value);
		gtk_tree_view_set_cursor(Tree,
			gtk_tree_model_get_path(model, &cur_items),
			gtk_tree_view_get_column(Tree, cur_column), FALSE);
		g_free(value);
		showdel = TRUE;
	}
	gtk_widget_set_sensitive(del_button, showdel);
	DBG("ret");*/
}

void newline(IMAGE *ima){
	FNAME();
	GtkTreeIter cur_items;
	gtk_list_store_append(ima->store, &cur_items);
	gtk_list_store_set(ima->store, &cur_items,
			L_TYPENM,	key_types[0],
			L_KEY,		"NEW_KEY",
			L_VAL,		"",
			L_COMM,		"",
			L_TYPE,		type_vals[0],
			-1);
	// Set cursor onto a first field
	gtk_tree_view_set_cursor(Tree,
			gtk_tree_model_get_path(GTK_TREE_MODEL(ima->store), &cur_items),
			gtk_tree_view_get_column(Tree, 0), TRUE);
	DBG("ret");
}

gboolean del_line(IMAGE *ima){
	GtkTreeIter iter;
	GtkTreeModel *model = GTK_TREE_MODEL(ima->store);
	gboolean ret = FALSE;
	FNAME();
	GtkTreeSelection *selection = gtk_tree_view_get_selection(Tree);
	GtkTreePath *path;
	if(gtk_tree_selection_get_selected(selection, NULL, &iter)){
		path = gtk_tree_model_get_path(model, &iter);
		ret = gtk_list_store_remove(ima->store, &iter);
		// check wherher  &cur_items was last
		if(gtk_list_store_iter_is_valid(ima->store, &iter)){
			gtk_tree_path_free(path);
			path = gtk_tree_model_get_path(model, &iter);
		}
		else gtk_tree_path_prev(path); // goto previous
		if(path && gtk_tree_model_get_iter_first(model, &iter))
				gtk_tree_view_set_cursor(Tree, path, NULL, FALSE);
		else DBG("No entries");
		gtk_tree_path_free(path);
	}
	DBG("ret");
	return ret;
}
/*
void show_menu(){
	GtkWidget *menu, *menuitem;
	gchar *text;
	int i;
	menu = gtk_menu_new();
	for(i = 0; i < 10; i++){
		text = g_strdup_printf("%d", i);
		menuitem = gtk_menu_item_new_with_label(text);
		gtk_menu_shell_append(GTK_MENU_SHELL(menu), menuitem);
		gtk_widget_show(menuitem);
		g_free(text);
	}
	gtk_menu_popup(GTK_MENU(menu), NULL, NULL, NULL, NULL, 3, 0);
}

gboolean treeview_button_press_cb(GtkWidget *wi, GdkEventButton *ev, gpointer user_data)
{
    GtkTreeSelection *sel;
   GtkTreeModel *model;
    GtkTreePath *path, *oldpath;
    GtkTreeViewColumn *column;
    GtkTreeIter iter;
    gchar *colname;

    // if there's no path where the click occurred...
    if(!gtk_tree_view_get_path_at_pos(GTK_TREE_VIEW(wi), ev->x, ev->y, &path, &column, NULL, NULL))
        return FALSE;    // ...return FALSE to allow the event to continue
	colname = strdup(gtk_tree_view_column_get_title(column));
    // get the tree view's selection
    sel = gtk_tree_view_get_selection(GTK_TREE_VIEW(wi));

    switch(ev->button)
    {
        // LMB
        case 1:
            // check to see if this is already the selected path
            if(gtk_tree_selection_get_selected(sel, &model, &iter))
            {
                oldpath = gtk_tree_model_get_path(model, &iter);
                if(gtk_tree_path_compare(oldpath, path) == 0)
                {
                    gtk_tree_path_free(oldpath);
                    gtk_tree_path_free(path);
                    return FALSE;
                }
                gtk_tree_path_free(oldpath);
            }

            // select the clicked path
            gtk_tree_selection_select_path(sel, path);
            return TRUE;
        // RMB
        case 3:
			if(strncmp(colname, "No.", 3)){
				gtk_tree_path_free(path);
				return TRUE;
			}
			show_menu();
//            gtk_menu_popup(GTK_MENU(menu), NULL, NULL, NULL, NULL, 3, ev->time);
            break;
    }


    // free our path
    gtk_tree_path_free(path);

    return FALSE;    // allow it to continue
}
*/

void close_dlg(){
	gtk_widget_destroy(window);
}

static gboolean press_key(IMAGE *ima, GdkEventKey * event){
	gboolean ret = FALSE;
	if(event->state & (GDK_CONTROL_MASK)){
		switch(event->keyval){
			// ctrl+n inserts a line
			case 'n': case 'N':
				newline(ima);
				ret = TRUE;
			break;
			// ctrl+d deletes current line
			case 'd': case 'D':
				ret = del_line(ima);
			break;
			// ctrl+w closes the window
			case 'w': case 'W':
				close_dlg();
		}
	}
	return ret;
}

void edit_headers(Window *parent){
	VALS_N = sizeof(key_types) / sizeof(gchar *);
	GtkWidget *treeview, *vbox, *hbox, *button;
	if(!parent->image || !parent->image->store){
		g_err(_("Open fits file first"));
		return;
	}
	window = gtk_window_new(GTK_WINDOW_TOPLEVEL);
	gtk_window_set_title(GTK_WINDOW(window), "List");
	gtk_signal_connect(GTK_OBJECT(window),"delete_event",GTK_SIGNAL_FUNC(gtk_false),NULL);
	gtk_signal_connect(GTK_OBJECT(window),"destroy",GTK_SIGNAL_FUNC(close_dlg), NULL);
	treeview = create_treeview(parent->image);
	sel = gtk_tree_view_get_selection(GTK_TREE_VIEW(treeview));
	g_signal_connect(G_OBJECT(sel), "changed", G_CALLBACK(list_changed), NULL);
	//g_signal_connect(G_OBJECT(treeview), "cursor-changed", G_CALLBACK(list_changed), NULL);

	vbox = gtk_vbox_new(FALSE, 5);
	gtk_container_add (GTK_CONTAINER (window), vbox);

	GtkWidget *sw = gtk_scrolled_window_new (NULL, NULL);
	gtk_scrolled_window_set_shadow_type(GTK_SCROLLED_WINDOW(sw), GTK_SHADOW_ETCHED_IN);
	gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(sw), GTK_POLICY_AUTOMATIC,
									GTK_POLICY_AUTOMATIC);
	gtk_box_pack_start(GTK_BOX(vbox), sw, TRUE, TRUE, 0);
	gtk_container_add(GTK_CONTAINER(sw), treeview);
	gtk_window_set_default_size(GTK_WINDOW (window), 640, 480);

	//~ gtk_box_pack_start(GTK_BOX(vbox), treeview, TRUE, TRUE, 0);

	hbox = gtk_hbox_new(FALSE, 5);
	gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);
	button = gtk_button_new_with_label(_("New entry (ctrl+n)"));
	g_signal_connect_swapped(G_OBJECT(button), "clicked",
			G_CALLBACK(newline), parent->image);
	gtk_box_pack_end(GTK_BOX(hbox), button, FALSE, FALSE, 0);

	del_button = gtk_button_new_with_label(_("Delete entry (ctrl+d)"));
	g_signal_connect_swapped(G_OBJECT(del_button), "clicked",
			G_CALLBACK(del_line), parent->image);
	gtk_box_pack_end(GTK_BOX(hbox), del_button, FALSE, FALSE, 0);
	gtk_widget_set_sensitive(del_button, FALSE);

	button = gtk_button_new_with_label(_("Close (ctrl+w)"));
	g_signal_connect(G_OBJECT(button), "clicked", G_CALLBACK(close_dlg), NULL);
	gtk_box_pack_end(GTK_BOX(hbox), button, FALSE, FALSE, 0);
	g_signal_connect_swapped(G_OBJECT(window), "key-press-event",
			G_CALLBACK(press_key), parent->image);
	run_modal_window(GTK_WINDOW(window), parent);
}
