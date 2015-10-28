/* SoX Memory allocation functions
 *
 * Copyright (c) 2005-2006 Reuben Thomas.  All rights reserved.
 * Copyright (C) 2008-10 lvqcl
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

#include "StdAfx.h"
#include "xmalloc.h"
/*#undef malloc
#undef calloc
#undef realloc
#undef free*/

#include <stdlib.h>

void (*malloc_error_handler)(void) = NULL;
void set_malloc_error_handler(void (*h)(void))
{
    malloc_error_handler = h;
}


/* Resize an allocated memory area; call throw_new_exception if not possible.
 *
 * For malloc, `If the size of the space requested is zero, the behavior is
 * implementation defined: either a null pointer is returned, or the
 * behavior is as if the size were some nonzero value, except that the
 * returned pointer shall not be used to access an object'
 */

void *lsx_realloc(void *ptr, size_t newsize)
{
    void * p;

    if (ptr && newsize == 0) { free(ptr); return NULL; }

    if ((p = realloc(ptr, newsize)) == NULL)
    {
        free(ptr);
        malloc_error_handler();
    }

    return p;
}

void lsx_free(void *ptr)
{
    free(ptr);
}

void* lsx_aligned_malloc(size_t size)
{
    void * p = _aligned_malloc(size, 16);
    if (p == NULL) malloc_error_handler();
    return p;
}

void * lsx_aligned_realloc(void * ptr, size_t newsize, size_t oldsize)
{
    void * p;
    if (ptr == NULL) return lsx_aligned_malloc(newsize); // & assert(oldsize==0)

    p = _aligned_realloc(ptr, newsize, 16);
    if (p) return p;

    p = _aligned_malloc(newsize, 16);
    if (p) memcpy(p, ptr, oldsize);
    _aligned_free(ptr);
    if (p == NULL) malloc_error_handler();
    return p;
}
