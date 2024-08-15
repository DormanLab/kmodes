#ifndef __ERROR_H__
#define __ERROR_H__

#include <stdio.h>
#include <stdarg.h>
#include <errno.h>
#ifdef USE_CURSES
#include <curses.h>
#endif

#include "kmodes_r.h"

/** Types of messages.
 */
enum {	NO_MSG,		/*!< no message */
	INFO_MSG,	/*!< informative message */
	DEBUG_MSG,	/*!< debugging message */
	WARNING_MSG,	/*!< warning message */
	ERROR_MSG	/*!< error message */
};

/** Types of errors.
 */
enum {	NO_ERROR,		/*!< no error */
	CUSTOM_ERROR,		/*!< customer error */
	NO_DATA,		/*!< necessary data not provided */
	MEMORY_ALLOCATION,	/*!< memory allocation error */
	FILE_NOT_FOUND,		/*!< file not found */
	FILE_OPEN_ERROR,	/*!< file open error */
	END_OF_FILE,		/*!< premature end of file error */
	FILE_FORMAT_ERROR,	/*!< invalid file format error */
	INVALID_CMDLINE,	/*!< invalid command line, specific below */
	INVALID_SUBCMD,		/*!< invalid subcommand */
	INVALID_CMD_OPTION,	/*!< invalid command-line option, e.g. -v */
	INVALID_CMD_PARAMETER,	/*!< invalid parameter to command-line argument */
	INVALID_USER_INPUT,	/*!< invalid user setup */
	INTERNAL_MISMATCH,	/*!< data does not match in some way it must */
	INTERNAL_ERROR,		/*!< internal error */
	CLUSTER_SIZE_OVERFLOW,	/*!< overflow in no. of clusters? */
	STATE_SPACE_OVERFLOW,	/*!< memory overflow due to state space size */
	OUT_OF_TIME,		/*!< ran out of time */
	MEMORY_USAGE_LIMIT,	/*!< ran out of memory */
	MEMCPY_ERROR,		/*!< error in call to memcpy() */
	EXCEED_ITERATIONS,	/*!< exceed some iteration limit */
	PIPE_OPEN_ERROR,	/*!< pipe open error */
	PIPE_READ_ERROR,	/*!< pipe read error */
	NUM_ERRORS		/*!< total number of errors */
};

/** Level of verbosity.
 */
enum {	ABSOLUTE_SILENCE,	/*!< only output through files */
	QUIET,			/*!< try to be quiet */
	MINIMAL,		/*!< minimal verbosity */
	VERBOSE,		/*!< verbose output */
	DEBUG_I,		/*!< debugging output */
	DEBUG_II,		/*!< debugging output */
	DEBUG_III,		/*!< debugging output */
	DEBUG_OVERRIDE		/*!< ignores global level */
};

extern int global_debug_level;

int message(FILE *, const char *, const char *, int, int, int, const char *, ...);
int r_message(const char *, const char *, int, int, int, const char *, ...);
#ifdef USE_CURSES
int wmessage(WINDOW *, const char *, const char *, int, int, int, const char *, ...);
#endif

/**
 * Print a formatted message to stderr. Have to avoid stderr for R.
 */
#ifdef MATHLIB_STANDALONE
#define mmessage(type, err, ...) message(stderr, __FILE__, __func__,  __LINE__, (type), (err),  __VA_ARGS__)
#else
#define mmessage(type, err, ...) r_message(__FILE__, __func__,  __LINE__, (type), (err),  __VA_ARGS__)
#endif

/**
 * Conditionally print a formatted message to stderr. Kludge a solution for R.
 */
#ifdef MATHLIB_STANDALONE
#define debug_msg(condition, level, ...) do {                                  \
	if ((condition) || ((level) && (level) <= global_debug_level))         \
		message(stderr, __FILE__, __func__, __LINE__, level >= DEBUG_I \
			? DEBUG_MSG : INFO_MSG, NO_ERROR, __VA_ARGS__);        \
} while (0)
#else
#define debug_msg(condition, level, ...) do {                                  \
	if ((condition) || ((level) && (level) <= global_debug_level)) {       \
		kmodes_printf("%s [%s::%s(%4d)]: ", level >= DEBUG_I ?         \
			"DEBUG" : "INFO", __FILE__, __func__, __LINE__);       \
		kmodes_printf(__VA_ARGS__);                                    \
	}                                                                      \
} while (0)
#endif

/**
 * Conditionally print a continuing message to stderr.
 */
#define debug_msg_cont(condition, level, ...) do {                             \
	if ((condition) || ((level) && (level) <= global_debug_level))         \
		kmodes_eprintf(__VA_ARGS__);                                   \
} while(0)

/**
 * Conditionally call a function.
 */
#define debug_call(condition, level, fxn_call) do {                            \
	if ((condition) || ((level) && (level) <= global_debug_level))         \
		(fxn_call);                                                    \
} while (0)

#endif
