#include <stdio.h>
#include <string.h>
#include <caml/mlvalues.h>
#include <caml/memory.h>
#include <caml/alloc.h>

/*
value cleanseq ( value ocamls )
{
	CAMLparam1 (ocamls);
	CAMLlocal1 (caml_strc);
// The original string
	char *s;
	s = String_val (ocamls);
	int l = strlen (s);
// The cleaned string

	char *strc;
	caml_strc = caml_alloc_string (l+1);
	strc = String_val (caml_strc);
	char *ss = strc;
	while (*s) {
		if (*s >= 'a' && *s <= 'z'){
			*ss++ = *s;
		}
		s++;
	}
	*ss = '\0';

// The returned value
	caml_strc = caml_copy_string (strc);
	CAMLreturn (caml_strc);
}

value ltrim ( value ocamls )
{
	CAMLparam1 (ocamls);
	CAMLlocal1 (caml_strc);
// The original string
	char *s;
	s = String_val (ocamls);
	int l = strlen (s);

// The trimed string
	char *strc;
	caml_strc = caml_alloc_string (l+1);
	strc = String_val (caml_strc);
	char *ss = strc;
	while (*s) {
		if (*s != ' ' && *s != '\t' && *s != '\n' && *s != '\r'){
			strcpy (ss, s);
			break;
		}
		s++;
	}

// The returned value
	caml_strc = caml_copy_string (strc);
	CAMLreturn (caml_strc);
}
*/

char *c_revcomp ( char *orig, char *dest )
{
	int l = strlen (orig);
	orig = orig + l - 1;
	while ( l > 0 ){
		switch (*orig) {
		case 'a':
			*dest++ = 't';
			break;
		case 'A' :
			*dest++ = 'T';
			break;
		case 'c' :
			*dest++ = 'g';
			break;
		case 'C' :
			*dest++ = 'G';
			break;
		case 'g' :
			*dest++ = 'c';
			break;
		case 'G' :
			*dest++ = 'C';
			break;
		case 't' :
			*dest++ = 'a';
			break;
		case 'T' :
			*dest++ = 'A';
			break;
		default :
			*dest++ = *orig;
			break;
		}
		orig--;
		l--;
	}
	*dest = '\0';
	return 0;
}

char *c_ltrim ( char *orig, char *dest )
{
	while (*orig) {
		if (*orig != ' ' && *orig != '\t' && *orig != '\n' && *orig != '\r'){
			strcpy (dest, orig);
			break;
		}
		orig++;
	}
	return 0;
}

char *c_cleanseq_light (char *orig, char *dest)
{
	while (*orig) {
          if ((*orig >= 'a' && *orig <= 'z') || (*orig >= 'A' && *orig <= 'Z') || (*orig == '-') || (*orig == '*')){
			*dest++ = *orig;
		}
		orig++;
	}
	*dest = '\0';
	return 0;
}

char *c_cleanseq (char *orig, char *dest)
{
	while (*orig) {
          if ((*orig >= 'a' && *orig <= 'z') || (*orig >= 'A' && *orig <= 'Z')){
			*dest++ = *orig;
		}
		orig++;
	}
	*dest = '\0';
	return 0;
}

value cleanseq ( value ocamls )
{
	CAMLparam1 (ocamls);
	CAMLlocal1 (caml_strc);
// The original string
	char *o;
	o = String_val (ocamls);
	int l = strlen (o);
// The cleaned string
	char *d = (char *) malloc ((l+1) * sizeof (char));
	c_cleanseq (o,d);
	caml_strc = caml_copy_string (d);
	free (d);
	CAMLreturn (caml_strc);
}

value cleanseq_light (value ocamls)
{
	CAMLparam1 (ocamls);
	CAMLlocal1 (caml_strc);
// The original string
	char *o;
	o = String_val (ocamls);
	int l = strlen (o);
// The cleaned string
	char *d = (char *) malloc ((l+1) * sizeof (char));
	c_cleanseq_light (o,d);
	caml_strc = caml_copy_string (d);
	free (d);
	CAMLreturn (caml_strc);
}

value ltrim ( value ocamls )
{
	CAMLparam1 (ocamls);
	CAMLlocal1 (caml_strc);
// The original string
	char *o;
	o = String_val (ocamls);
	int l = strlen (o);

// The trimed string
	char *d = (char *) malloc ((l+1) * sizeof (char));
	c_ltrim (o,d);
	caml_strc = caml_copy_string (d);
	free (d);
	CAMLreturn (caml_strc);
}

value revcomp ( value ocamls )
{
	CAMLparam1 (ocamls);
	CAMLlocal1 (caml_strc);
	char *o;
	o = String_val (ocamls);
	int l = strlen (o);

	char *d = (char *) malloc ((l+1) * sizeof (char));
	c_revcomp (o,d);
	caml_strc = caml_copy_string (d);
	free (d);
	CAMLreturn (caml_strc);
}

