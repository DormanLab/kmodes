//
//  entropy.h
//  kmodes
//
//  Created by Yudi Zhang on 4/19/21.
//

#ifndef __ENTROPY_H__
#define __ENTROPY_H__

#include "run_kmodes.h"

int perturb(data *dat, options *opt);
int initialize_by_abundance(data *dat, unsigned int K, unsigned int k1, data_t **seeds, unsigned int *sd_idx);
int mask_nhash(data *dat, options *opt);
int perturb_general(data *dat, options *opt);

#endif /* __ENTROPY_H__ */
