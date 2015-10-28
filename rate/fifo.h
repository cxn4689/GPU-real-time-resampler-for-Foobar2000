/* Addressable FIFO buffer    Copyright (c) 2007 robs@users.sourceforge.net
 * Copyright (c) 2010-11 lvqcl
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

#ifndef fifo_included
#define fifo_included

#include <string.h>
#include "xmalloc.h"
#include "util.h"

#ifndef FIFO_SIZE_T
#define FIFO_SIZE_T size_t
#endif

/* _mm_lddqu_si128: "A 32-byte block is loaded"
This intrinsic can read extra 8 bytes (for arrays of doubles) so the size of buffers was increased by RESERV == 8
I think that this is not necessary but not sure... Should read Intel/AMD docs or do some tests */

/* _mm_loadu_ps:
this intrinsic can read extra 8 bytes (for arrays of floats) so the size of buffers was increased by RESERV == 8*/

#define RESERV 8
#define MASK (((size_t)-1)^0xf)
#define CEIL16(x) (((x)+0xf)&MASK)
#define SHIFT(x) (CEIL16(x) - (x))

typedef struct {
  char* data;
  size_t allocation;   /* Number of bytes allocated for data. */
  size_t item_size;    /* Size of each item in data */
  size_t begin;        /* Offset of the first byte to read. */
  size_t end;          /* 1 + Offset of the last byte byte to read. */
} fifo_t;

#define FIFO_MIN 0x08000

UNUSED static __inline void fifo_clear(fifo_t * f)
{
  f->end = f->begin = 0;
}

UNUSED static void * fifo_reserve(fifo_t * f, FIFO_SIZE_T n)
{
    size_t old_allocation = f->allocation;

    n *= f->item_size;

    if (f->begin == f->end)
        fifo_clear(f);

    for (;;)
    {
        if (f->end + n <= f->allocation)
        {
            void *p = f->data + f->end;

            f->end += n;
            return p;
        }
        if (f->begin > FIFO_MIN)
        {
            memmove(f->data, f->data + f->begin, f->end - f->begin);
            f->end -= f->begin;
            f->begin = 0;
            continue;
        }
        f->allocation += n;
        f->data = (char*)lsx_aligned_realloc(f->data, f->allocation + RESERV, old_allocation/* + RESERV*/);
        if (f->data == NULL) {f->allocation = 0; fifo_clear(f); return NULL; }
    }
}

UNUSED static void * fifo_reserve_aligned(fifo_t * f, FIFO_SIZE_T n)
{
    size_t old_allocation = f->allocation;
    n *= f->item_size;

    if (f->begin == f->end)
        fifo_clear(f);

    for (;;)
    {
        if (CEIL16(f->end) + n <= f->allocation)
        {
            void * p;
            if (SHIFT(f->end))
            {
                size_t shift = SHIFT(f->end - f->begin);
                memmove(f->data + shift, f->data + f->begin, f->end - f->begin);
                f->end -= f->begin;
                f->end += shift;   // f->end = f->end - f->begin + shift;
                f->begin = shift;
            }
            p = f->data + f->end;

            f->end += n;
            return p;
        }

        if (f->begin > FIFO_MIN)
        {
            size_t shift = SHIFT(f->end - f->begin);
            memmove(f->data + shift, f->data + f->begin, f->end - f->begin);
            f->end -= f->begin;
            f->end += shift;   // f->end = f->end - f->begin + shift;
            f->begin = shift;
            continue;
        }

        f->allocation += n + SHIFT(f->end);
        f->data = (char*)lsx_aligned_realloc(f->data, f->allocation + RESERV, old_allocation/* + RESERV*/);
        if (f->data == NULL) {f->allocation = 0; fifo_clear(f); return NULL; }
    }
}

UNUSED static __inline void * fifo_write(fifo_t * f, FIFO_SIZE_T n, void const * data)
{
  void * s = fifo_reserve(f, n);
  if (s == NULL) return NULL;
  if (data)
    memcpy(s, data, n * f->item_size);
  else
    memset(s, 0, n * f->item_size);
  return s;
}

UNUSED static __inline void fifo_trim_to(fifo_t * f, FIFO_SIZE_T n)
{
  n *= f->item_size;
  f->end = f->begin + n;
}

UNUSED static __inline void fifo_trim_by(fifo_t * f, FIFO_SIZE_T n)
{
  n *= f->item_size;
  f->end -= n;
}

UNUSED static __inline FIFO_SIZE_T fifo_occupancy(fifo_t * f)
{
  return (f->end - f->begin) / f->item_size;
}

UNUSED static __inline void * fifo_read(fifo_t * f, FIFO_SIZE_T n, void * data)
{
  char * ret = f->data + f->begin;
  n *= f->item_size;
  if (n > (FIFO_SIZE_T)(f->end - f->begin))
    return NULL;
  if (data)
    memcpy(data, ret, (size_t)n);
  f->begin += n;
  return ret;
}

/*#define fifo_read_ptr(f) fifo_read(f, (FIFO_SIZE_T)0, NULL)*/

UNUSED static __inline void * fifo_read_ptr(fifo_t * f)
{
  return f->data + f->begin;
}

UNUSED static __inline void fifo_delete(fifo_t * f)
{
  lsx_aligned_free(f->data);
}

UNUSED static int fifo_create(fifo_t * f, FIFO_SIZE_T item_size)
{
  f->item_size = item_size;
  f->allocation = FIFO_MIN;
  fifo_clear(f);
  f->data = (char*)lsx_aligned_malloc(f->allocation + RESERV);
  if (f->data == NULL) {
    f->allocation = 0;
    return 0;
  }
  return 1;
}

#endif
