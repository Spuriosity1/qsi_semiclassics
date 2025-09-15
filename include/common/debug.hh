#ifndef DEBUG_HH
#define DEBUG_HH

#include <stdio.h>

#ifdef VERBOSE
  #define db_printerr(...) fprintf(stderr, __VA_ARGS__)
#else
  #define db_printerr(...)
#endif

#endif