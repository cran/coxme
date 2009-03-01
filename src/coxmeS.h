/* $Id$ */
/*
**   A few things to make code work on both Splus and R
*/
#include "R.h"
#include "Rinternals.h"

#ifdef USING_R
#define S_EVALUATOR   /* turn this into a blank line in R */
  /*typedef int Sint   now a part of std R */
#define ALLOC(a,b)  R_alloc(a,b)

#else
typedef long Sint
/*
** Memory defined with S_alloc is removed automatically by S.
**  That with "CALLOC" I have to remove myself.  Use the
**  latter for objects that need to to persist between calls
*/
#if( defined(SPLUS_VERSION) && SPLUS_VERSION >= 5000)
#define ALLOC(a,b)  S_alloc(a,b,S_evaluator)
#define CALLOC(a,b) S_ok_calloc((size_t)(a), b, S_evaluator)
#endif
#endif
