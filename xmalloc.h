/* libSoX Memory allocation functions
 *
 * Copyright (c) 2005-2006 Reuben Thomas.  All rights reserved.
 * Copyright (c) 2010 lvqcl
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

#ifndef LSX_MALLOC_H
#define LSX_MALLOC_H

//#include <stddef.h>
#include <string.h>
#include <malloc.h>


#ifdef __cplusplus
extern "C" {
#endif

void set_malloc_error_handler(void (*h)(void));

#if defined(__GNUC__)
void __cdecl _aligned_free(void * memblock);
#endif

void * lsx_realloc(void * ptr, size_t newsize);
static __inline void* lsx_malloc(size_t size) { return lsx_realloc(NULL, size); }
static __inline void* lsx_calloc(size_t n ,size_t s) { if (n*s) return memset(lsx_realloc(NULL, n*s), 0, n*s); else return NULL; }
void lsx_free(void* p);
/*static __inline void lsx_free(void* p) { free(p); }*/


void * lsx_aligned_malloc(size_t size);
void * lsx_aligned_realloc(void * ptr, size_t newsize, size_t oldsize);
static __inline void* lsx_aligned_calloc(size_t n, size_t s) { if (n*s) return memset(lsx_aligned_malloc(n*s), 0, n*s); else return NULL; }
static __inline void  lsx_aligned_free(void* p) { _aligned_free(p); }

/*#define DONT_USE(x) { enum {do_not_use_##x = 1/ }; }*/
/*#define malloc { enum {do_not_use_malloc = 1/ }; }
#define calloc { enum {do_not_use_calloc = 1/ }; }
#define realloc { enum {do_not_use_realloc = 1/ }; }
#define free { enum {do_not_use_free = 1/ }; }*/

#ifdef __cplusplus
}
#endif

#endif
