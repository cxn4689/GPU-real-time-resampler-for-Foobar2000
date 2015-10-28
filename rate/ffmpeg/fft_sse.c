/*
 * FFT/MDCT transform with SSE optimizations
 * Copyright (c) 2008 Loren Merritt
 *
 * This file is part of FFmpeg.
 *
 * FFmpeg is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * FFmpeg is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with FFmpeg; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
 */

#include "ffti.h"

void ff_fft_dispatch_interleave_sse(FFTComplex *z, int nbits);

#if !defined (__GNUC__)
#include <xmmintrin.h>
#endif

#if defined (__GNUC__)
#if defined (__GNUC__)
void ff_fft_calc_sse(FFTContext *s, FFTComplex *z)
{
    int n = 1 << s->nbits;

    ff_fft_dispatch_interleave_sse(z, s->nbits);

    if(n <= 16) {
        x86_reg i = -8*n;
        __asm__ volatile(
            "1: \n"
            "movaps     (%0,%1), %%xmm0 \n"
            "movaps      %%xmm0, %%xmm1 \n"
            "unpcklps 16(%0,%1), %%xmm0 \n"
            "unpckhps 16(%0,%1), %%xmm1 \n"
            "movaps      %%xmm0,   (%0,%1) \n"
            "movaps      %%xmm1, 16(%0,%1) \n"
            "add $32, %0 \n"
            "jl 1b \n"
            :"+r"(i)
            :"r"(z+n)
            :"memory"
        );
    }
}
#else
void ff_fft_calc_sse(FFTContext *s, FFTComplex *z)
{
    int n = 1 << s->nbits;
    ff_fft_dispatch_interleave_sse(z, s->nbits);

    if(n <= 16)
    {
        __m128 xmm0, xmm1;
        int i = -2*n;
        float* p = (float*)&z[n];
        do
        {
            xmm0 = _mm_load_ps(&p[i]);
            xmm1 = xmm0;
            xmm0 = _mm_unpacklo_ps(xmm0, *(__m128*)&p[i+4]);
            xmm1 = _mm_unpackhi_ps(xmm1, *(__m128*)&p[i+4]);
            _mm_store_ps(&p[i], xmm0);
            _mm_store_ps(&p[i+4], xmm1);
            i+=8;
        }
        while(i<0);
    }
}
#endif
#endif

#if defined (__GNUC__)
#if defined (__GNUC__)
void ff_fft_permute_sse(FFTContext *s, FFTComplex *z, FFTComplex *tmp_buf)
{
    int n = 1 << s->nbits;
    int i;
    for(i=0; i<n; i+=2) {
        __asm__ volatile(
            "movaps %2, %%xmm0 \n"
            "movlps %%xmm0, %0 \n"
            "movhps %%xmm0, %1 \n"
            :"=m"(tmp_buf[s->revtab[i]]),
             "=m"(tmp_buf[s->revtab[i+1]])
            :"m"(z[i])
        );
    }
    memcpy(z, tmp_buf, n*sizeof(FFTComplex));
}
#else
void ff_fft_permute_sse(FFTContext *s, FFTComplex *z, FFTComplex *tmp_buf)
{
    int n = 1 << s->nbits;
    int i;
    __m128 xmm0;
    for(i=0; i<n; i+=2)
    {
        xmm0 = _mm_load_ps((float*)&z[i]);
        _mm_storel_pi((__m64*)(&tmp_buf[s->revtab[i]]), xmm0);
        _mm_storeh_pi((__m64*)(&tmp_buf[s->revtab[i+1]]), xmm0);
    }
    memcpy(z, tmp_buf, n*sizeof(FFTComplex));
}
#endif
#endif
