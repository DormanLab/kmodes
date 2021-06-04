#ifndef __KMODES_R_H__
#define __KMODES_R_H__

#ifdef MATHLIB_STANDALONE

#define kmodes_eprintf(...) fprintf(stderr, __VA_ARGS__)
#define kmodes_fprintf(fp, ...) fprintf((fp), __VA_ARGS__)
#define kmodes_printf(...) printf(__VA_ARGS__)
#define kmodes_vfprintf(fp, ...) vfprintf((fp), __VA_ARGS__)
#define kmodes_vprintf(...) vprintf(__VA_ARGS__)

#else

#include <R.h>
#include <R_ext/Memory.h>
#define kmodes_fprintf(fp, ...) REprintf(__VA_ARGS__)
#define kmodes_eprintf(...) REprintf(__VA_ARGS__)
#define kmodes_printf(...) Rprintf(__VA_ARGS__)
#define kmodes_vfprintf(fp, ...) REvprintf(__VA_ARGS__)
#define kmodes_vprintf(...) Rvprintf(__VA_ARGS__)

#endif


#endif
