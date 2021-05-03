#ifndef __SAMPLE_H__
#define __SAMPLE_H__

#include "constants.h"

void sample(unsigned int N, unsigned int n, unsigned int *idx);
int heap_sample(unsigned int n, unsigned int k, double *w, unsigned int *idx, unsigned int indic);
int random_sample(size_t N, size_t n,  size_t *D_idx, size_t *s_idx);

#endif
