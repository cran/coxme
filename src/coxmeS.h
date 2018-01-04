/*
**   A common header
*/
#include "R.h"
#include "Rinternals.h"

#ifdef USING_R
#define ALLOC(a,b)  R_alloc(a,b)
#endif
