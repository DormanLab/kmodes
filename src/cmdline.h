/**
 * @file cmdline.h
 * @author Karin Dorman, kdorman@iastate.edu
 * @date Mon Dec 10 16:20:59 CST 2012
 *
 * Header file for cmdline.c.
 */

#ifndef __H_CMDLINE__
#define __H_CMDLINE__

#include <errno.h>
#include <stdio.h>
#include <string.h>

/**
 * Check if option argument.
 *
 * @param argc	number of command line arguments
 * @param argv	command line arguments
 * @param aidx	command line argument currently checking for argument
 * @return	error status
 */
inline int has_argument(int argc, const char **argv, int aidx)
{
	if (aidx + 1 == argc || argv[aidx + 1][0] == '-')
		return 0;
	return 1;
} /* has_argument */

int usage_error(const char **argv, int i, void *obj);

/* check single entries */
int is_numeric(const char *argv);

/* read single entries */
int read_int(int argc, const char **argv, int i, void *obj);
unsigned int read_uint(int argc, const char **argv, int i, void *obj);
long read_long(int argc, const char **argv, int i, void *obj);
unsigned long read_ulong(int argc, const char **argv, int i, void *obj);
unsigned long read_ulonglong(int argc, const char **argv, int i, void *obj);
double read_cmdline_double(int argc, const char **argv, int i, void *obj);
char read_cmdline_char(int argc, const char **argv, int i, void *obj);

/* read multiple enries */
unsigned int read_cmdline_doubles(unsigned int argc, const char **argv, unsigned int i, double **dret, void *obj);
unsigned int read_ulongs(unsigned int argc, const char **argv, unsigned int i, unsigned long **sret, void *obj);
unsigned int read_cmdline_uints(unsigned int argc, const char **argv, unsigned int i, unsigned int **uiret, void *obj);
unsigned int read_cmdline_strings(unsigned int argc, const char ** argv, unsigned int i, char const ***sret, void *obj);

void print_usage(const char *);
void fprint_usage(FILE *fp, const char *program_name, void *obj);

#endif
