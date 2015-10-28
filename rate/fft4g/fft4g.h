/* This library is free software; you can redistribute it and/or modify it
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

#ifndef FFT4G_H
#define FFT4G_H

#define _USE_MATH_DEFINES
#include <math.h>

typedef struct FFTcontext_tag
{
    int n;
    //int bits;
    int *br;
    double *sc;
}FFTcontext;

#ifdef __cplusplus
extern "C" {
#endif

int  lsx_rdft_init(int n, FFTcontext* z);
void lsx_rdft_close(FFTcontext* z);

void lsx_rdft_SSE3(int isgn, double *a, FFTcontext z);
void lsx_rdft_generic(int isgn, double *a, FFTcontext z);

/* Over-allocate h by 2 to use these macros */
#define LSX_PACK(h, n)   h[1] = h[n]
#define LSX_UNPACK(h, n) h[n] = h[1], h[n + 1] = h[1] = 0

#ifdef __cplusplus
}
#endif

#ifdef __cplusplus
template<typename T> static inline void lsx_pack(T* a, unsigned n) { a[1] = a[n]; }
template<typename T> static inline void lsx_unpack(T* a, unsigned n) { a[n] = a[1]; a[n+1] = a[1] = (T)0; }
#endif

#endif
