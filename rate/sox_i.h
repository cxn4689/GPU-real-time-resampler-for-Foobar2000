/* libSoX Internal header
 *
 *   This file is meant for libSoX internal use only
 *
 * Copyright 2001-2008 Chris Bagwell and SoX Contributors
 *
 * This source code is freely redistributable and may be used for
 * any purpose.  This copyright notice must be maintained.
 * Chris Bagwell And SoX Contributors are not responsible for
 * the consequences of using this software.
 * (C) 2008-11 lvqcl - some (hmm...) changes
 */

#ifndef SOX_I_H
#define SOX_I_H

#if defined(__GNUC__)
extern __inline int __attribute__((__gnu_inline__, __always_inline__, __artificial__)) _bit_scan_reverse (int __X)
{
    return __builtin_ia32_bsrsi (__X);
}
static __inline int dlog2(int n)
{
    return _bit_scan_reverse(n);
}
#else
#include <intrin.h>
#pragma intrinsic(_BitScanReverse)
static __inline int dlog2(int n)
{
    unsigned long l;
    _BitScanReverse(&l, n);
    return (int)l;
}
#endif

#include "sox.h"
#include "util.h"
#include "fft4g/fft4g.h"
#include "ffmpeg/fft_ffmpeg.h"

int lsx_set_dft_length(int num_taps);
int x_lsx_set_dft_length(int num_taps);

extern FFTcontext fftx[17];

extern void (*lsx_rdft)(int isgn, double *a, FFTcontext z);

static __inline void lsx_safe_rdft_SSE3(int n, int type, double * d)
{
    lsx_rdft_SSE3(type, d, fftx[dlog2(n)]);
}

static __inline void lsx_safe_rdft_generic(int n, int type, double * d)
{
    lsx_rdft_generic(type, d, fftx[dlog2(n)]);
}

static __inline void x_lsx_safe_rdft_SSE3(int bits, int type, double * d)
{
    lsx_rdft_SSE3(type, d, fftx[bits]);
}

static __inline void x_lsx_safe_rdft_generic(int bits, int type, double * d)
{
    lsx_rdft_generic(type, d, fftx[bits]);
}

/*static __inline void lsx_safe_rdft_A(int n, int type, double * d)
{
    if(n<=65536)
    {
        lsx_rdft(type, d, fftx[dlog2(n)]);
        return;
    }
    else
    {
        FFTcontext z;
        lsx_rdft_init(n, &z);
        lsx_rdft(type, d, z);
        lsx_rdft_close(&z);
    }
}*/

extern RDFTContext ff_fwd[17];
extern RDFTContext ff_bkd[17];

static __inline void ff_rdft_SSE(int n, int type, float * d, FFTComplex * tmp_buf)
{
    ff_rdft_calc_sse( &( type==1 ? ff_fwd : ff_bkd )[dlog2(n)], d, tmp_buf );
}

static __inline void ff_rdft_generic(int n, int type, float * d, FFTComplex * tmp_buf)
{
    ff_rdft_calc_c( &( type==1 ? ff_fwd : ff_bkd )[dlog2(n)], d, tmp_buf );
}

static __inline void x_ff_rdft_SSE(int bits, int type, float * d, FFTComplex * tmp_buf)
{
    ff_rdft_calc_sse( &( type==1 ? ff_fwd : ff_bkd )[bits], d, tmp_buf );
}

static __inline void x_ff_rdft_generic(int bits, int type, float * d, FFTComplex * tmp_buf)
{
    ff_rdft_calc_c( &( type==1 ? ff_fwd : ff_bkd )[bits], d, tmp_buf );
}

double * lsx_design_lpf(
    double Fp,      /* End of pass-band; ~= 0.01dB point */
    double Fc,      /* Start of stop-band */
    double Fn,      /* Nyquist freq; e.g. 0.5, 1, PI */
    sox_bool allow_aliasing,
    double att,     /* Stop-band attenuation in dB */
    int * num_taps, /* (Single phase.)  0: value will be estimated */
    int k,          /* Number of phases; 0 for single-phase */
    double beta);

void lsx_fir_to_phase(double * * h, int * len, int * post_len, double phase0);

#define LSX_TO_6dB .5869
#define LSX_TO_3dB ((2/3.) * (.5 + LSX_TO_6dB))           /* 0.7246 */
#define LSX_MAX_TBW0 36.
#define LSX_MAX_TBW0A (LSX_MAX_TBW0 / (1 + LSX_TO_3dB))   /* 20.8744 */
//#define LSX_MAX_TBW3 floor(LSX_MAX_TBW0 * LSX_TO_3dB)   /* 26 => 100-26 = 74% */
//#define LSX_MAX_TBW3A floor(LSX_MAX_TBW0A * LSX_TO_3dB) /* 15 => 100-15 = 85% */


#endif
