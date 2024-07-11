#ifndef __RMATH_H__
#define __RMATH_H__

void set_seed(unsigned int i1, unsigned int i2);

double unif_rand(void);
double exp_rand(void);
double norm_rand(void);

double rgamma(double a, double scale);

double fmin2(double x, double y);
double fmax2(double x, double y);

double qnorm5(double p, double mu, double sigma, int lower_tail, int log_p);


#endif
