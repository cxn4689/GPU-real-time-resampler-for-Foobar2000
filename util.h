/* General purpose, i.e. non SoX specific, utility functions and macros.
 *
 * (c) 2006-8 Chris Bagwell and SoX contributors
 *
 * This library is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation; either version 2.1 of the License, or (at
 * your option) any later version.
 *
 * This library is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser
 * General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this library; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

#define HAVE_SYS_TYPES_H 1

#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h> /* For off_t not found in stdio.h */
#endif

#ifdef HAVE_SYS_STAT_H
#include <sys/stat.h> /* Needs to be included before we redefine off_t. */
#endif

#include "xmalloc.h"

/*---------------------------- Portability stuff -----------------------------*/

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

#ifdef _MSC_VER
#define __STDC__ 1
#define inline __inline
#elif defined(__MINGW32__)
#endif

#ifdef WORDS_BIGENDIAN
  #define MACHINE_IS_BIGENDIAN 1
  #define MACHINE_IS_LITTLEENDIAN 0
#else
  #define MACHINE_IS_BIGENDIAN 0
  #define MACHINE_IS_LITTLEENDIAN 1
#endif

/*--------------------------- Language extensions ----------------------------*/

/* Compile-time ("static") assertion */
/*   e.g. assert_static(sizeof(int) >= 4, int_type_too_small)    */
#define assert_static(e,f) enum {assert_static__##f = 1/(e)}
#define array_length(a) (sizeof(a)/sizeof(a[0]))

/*------------------------------- Maths stuff --------------------------------*/

#include <math.h>

#ifdef min
#undef min
#endif
#define min(a, b) ((a) <= (b) ? (a) : (b))

#ifdef max
#undef max
#endif
#define max(a, b) ((a) >= (b) ? (a) : (b))

#define range_limit(x, lower, upper) (min(max((x), (lower)), (upper)))

#ifndef M_PI
#define M_PI    3.14159265358979323846
#endif
#ifndef M_PI_2
#define M_PI_2  1.57079632679489661923  /* pi/2 */
#endif
#ifndef M_LN10
#define M_LN10  2.30258509299404568402  /* natural log of 10 */
#endif
#ifndef M_LN2
#define M_LN2  0.693147180559945309417  /* natural log of 2 */
#endif
#ifndef M_SQRT2
#define M_SQRT2  1.41421356237309504880
#endif

#define sqr(a) ((a) * (a))
#define sign(x) ((x) < 0? -1 : 1)
