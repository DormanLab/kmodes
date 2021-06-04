#ifndef TIMING_MACH_H
#define TIMING_MACH_H

/* C99 check: C99 or greater required for inline functions */
#if defined(__STDC__)
# if defined(__STDC_VERSION__)
#  if (__STDC_VERSION__ >= 199901L)
#   define TIMING_C99
#  endif
# endif
#endif

#include <time.h>

#define TIMING_GIGA (1000000000)
#define TIMING_NANO (1e-9)

/* MacOSX: old MacOSX (<10.12) not POSIX compliant, but can time. */
#if defined(__MACH__) && !defined(CLOCK_REALTIME)

#define __OLD_MACH__

/* only CLOCK_REALTIME and CLOCK_MONOTONIC are emulated */
#ifndef CLOCK_REALTIME
# define CLOCK_REALTIME 0
#endif
#ifndef CLOCK_MONOTONIC
# define CLOCK_MONOTONIC 1
#endif

/* typdef POSIX clockid_t */
typedef int clockid_t;

/* clock_gettime - emulate POSIX */
int clock_gettime(const clockid_t id, struct timespec *tspec);

/* POSIX-compliant */
#elif defined(__unix__) || defined(__MACH__)

#include <unistd.h>

#if _POSIX_VERSION < 199309L
# define __NON_POSIX__
#endif

/* Not POSIX-compliant: resort to C standard */
#else			

#define __NON_POSIX__

struct timespec {
	clock_t time;
};

#endif

/* helper function prototypes */
extern void timespec_mark(struct timespec *ts);
extern double timespec_elapsed(const struct timespec *ts_start);

#define MARK_TIME(a)	timespec_mark((a))
#define ELAP_TIME(a)	timespec_elapsed((a))
#define TIME_STRUCT	struct timespec

#endif
