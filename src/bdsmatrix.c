/*
** The source code for establishing links to the bdsmatrix C routines is
**  found in the .h file below, which is included automatically from the 
**  inst/include directory of the bdsmatrix package by the  "LinkingTo"
**  line in ../DESCRIPTION.
** All we need is a hook (this file) to cause the lines to be compiled.
*/
/* 
** Later note: I've added a copy of bdsmatrix/inst/include/bdsmatrix_stub.h
**  to this directory.  It makes this more self contained, but has the danger
**  that if the bdsmatrix routines change I may have an error.  (But then, if
**  they change we have more updates to do than just this.)
*/
#include "bdsmatrix_stub.h"
