/* libSoX Library Public Interface
 *
 * Copyright 1999-2009 Chris Bagwell and SoX Contributors.
 *
 * This source code is freely redistributable and may be used for
 * any purpose.  This copyright notice must be maintained.
 * Chris Bagwell And SoX Contributors are not responsible for
 * the consequences of using this software.
 */

#ifndef SOX_H
#define SOX_H

#include <stddef.h> /* Ensure NULL etc. are available throughout SoX */
#include <stdlib.h>

#ifdef HAVE_STDINT_H
  #include <stdint.h>
#else
  #ifdef HAVE_INTTYPES_H
    #include <inttypes.h>
  #else
    #ifdef _MSC_VER
      typedef __int64 int64_t;
      typedef unsigned __int64 uint64_t;
    #else
      typedef long long int64_t;
      typedef unsigned long long uint64_t;
    #endif
    typedef long int32_t;
    typedef unsigned long uint32_t;
  #endif
#endif

/* Avoid warnings about unused parameters. */
#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

#define SOX_SUCCESS 0
#define SOX_EOF (-1)             /* End Of File or other error */

/* Boolean type, assignment (but not necessarily binary) compatible with C++ bool */
typedef enum {sox_false = 0, sox_true = 1} sox_bool;

#endif
