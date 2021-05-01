//
//  entropy.h
//  kmodes
//
//  Created by Yudi Zhang on 4/19/21.
//

#ifndef entropy_h
#define entropy_h

#include <float.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "kmodes.h"
#include "array.h"
#include "order.h"
#include "error.h"
#include "math.h"
#include "io.h"
#include "io_kmodes.h"
#include "sample.h"
#include "run_kmodes.h"

int perturbs(data *dat, options *opt);
int initialize_by_abundance(data *dat, unsigned int K, unsigned int k1, data_t **seeds, unsigned int *sd_idx);
int mask_nhash(data *dat, unsigned int K, int abunk);

#endif /* entropy_h */
