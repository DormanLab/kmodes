#ifndef __MATH_NOT_H__
#define __MATH_NOT_H__

#include <stddef.h>

/**
 * Type of vectorization, by row or column.
 */
enum {
	ROW_ORDER,
	COLUMN_ORDER
};

/**
 * Simple mathematical functions.
 */
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define SQ(x) ((x) * (x))
#define SIGN(x) (((x) > 0) - ((x) < 0))

#endif
