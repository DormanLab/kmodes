/**
 * @file timing_mach.c
 *
 * Attempt at universal timing code, using real time and higher precision 
 * when possible, and resort to clock() and CLOCKS_PER_SEC in the worst
 * case.
 */

#include "timing_mach.h"

/* Old MacOSX: see timing_mach.h */
#ifdef __OLD_MACH__

#include <mach/mach_time.h>
#include <mach/mach.h>
#include <mach/clock.h>

/* timing struct for osx */
static struct TimingMach {
    mach_timebase_info_data_t timebase;
    clock_serv_t cclock;
} timing_mach_g;

/* mach clock port */
extern mach_port_t clock_port;

int clock_gettime(clockid_t id, struct timespec *tspec) {
    mach_timespec_t mts;
    if (id == CLOCK_REALTIME) {
        if (clock_get_time(timing_mach_g.cclock, &mts) != 0)
            return -1;
        tspec->tv_sec = mts.tv_sec;
        tspec->tv_nsec = mts.tv_nsec;
    } else if (id == CLOCK_MONOTONIC) {
        if (clock_get_time(clock_port, &mts) != 0)
            return -1;
        tspec->tv_sec = mts.tv_sec;
        tspec->tv_nsec = mts.tv_nsec;
    } else {
        /* only CLOCK_MONOTONIC and CLOCK_REALTIME clocks supported */
        return -1;
    }
    return 0;
}

/* newer POSIX machines */
#elif !defined(__NON_POSIX__)

/* no inline functions if not at least C99 */
#ifndef TIMING_C99
# define inline
#endif

inline void timespec_mark(struct timespec *ts)
{
	clock_gettime(CLOCK_MONOTONIC, ts);
} /* timespec_mark */

inline double timespec_elapsed(const struct timespec *ts_start)
{
	struct timespec ts_stop, ts_temp;

	timespec_mark(&ts_stop);

	if (ts_stop.tv_nsec - ts_start->tv_nsec < 0) {
		ts_temp.tv_sec = ts_stop.tv_sec - ts_start->tv_sec - 1;
		ts_temp.tv_nsec = TIMING_GIGA + ts_stop.tv_nsec - ts_start->tv_nsec;
	} else {
		ts_temp.tv_sec = ts_stop.tv_sec - ts_start->tv_sec;
		ts_temp.tv_nsec = ts_stop.tv_nsec - ts_start->tv_nsec;
	}

	return ((double) ts_temp.tv_sec) + ((double) ts_temp.tv_nsec) * TIMING_NANO;
} /* timespec_elapsed */

/* clean up define 'inline' */
#ifndef TIMING_C99
# undef inline
#endif

/* older POSIX or non-POSIX machines: resort to C standard */
#else

#ifndef TIMING_C99
# define inline
#endif

inline void timespec_mark(struct timespec *ts)
{
	ts->time = clock();
} /* timespec_mark */

inline double timespec_elapsed(const struct timespec *ts_start)
{
	struct timespec ts_stop;

	timespec_mark(&ts_stop);

	return (double) (ts_stop.time - ts_start.time) / CLOCKS_PER_SECOND;
} /* timespec_elapsed */

/* clean up define 'inline' */
#ifndef TIMING_C99
# undef inline
#endif

#endif	/* OS choice */
