/*
 * (I)RDFT transforms
 * Copyright (c) 2009 Alex Converse <alex dot converse at gmail dot com>
 * Copyright (c) 2010-11 lvqcl
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

/**
 * @file
 * (Inverse) Real Discrete Fourier Transforms.
 */

/* sin(2*pi*x/n) for 0<=x<n/4, followed by n/2<=x<3n/4 */
#if !CONFIG_HARDCODED_TABLES
SINTABLE(16);
SINTABLE(32);
SINTABLE(64);
SINTABLE(128);
SINTABLE(256);
SINTABLE(512);
SINTABLE(1024);
SINTABLE(2048);
SINTABLE(4096);
SINTABLE(8192);
SINTABLE(16384);
SINTABLE(32768);
SINTABLE(65536);
#endif
static SINTABLE_CONST FFTSample * const ff_sin_tabs[] = {
    NULL, NULL, NULL, NULL,
    ff_sin_16, ff_sin_32, ff_sin_64, ff_sin_128, ff_sin_256, ff_sin_512, ff_sin_1024,
    ff_sin_2048, ff_sin_4096, ff_sin_8192, ff_sin_16384, ff_sin_32768, ff_sin_65536,
};

/** Map one real FFT into two parallel real even and odd FFTs. Then interleave
 * the two real FFTs into one complex FFT. Unmangle the results.
 * ref: http://www.engineeringproductivitytools.com/stuff/T0001/PT10.HTM
 */
void ff_rdft_calc_c(RDFTContext* s, FFTSample* data, FFTComplex *tmp_buf)
{
    int i, i1, i2;
    FFTComplex ev, od;
    const int n = 1 << s->nbits;
    const float k1 = 0.5f;
    const float k2 = (float)(0.5f - s->inverse);
    const FFTSample *tcos = s->tcos;
    const FFTSample *tsin = s->tsin;

    if (!s->inverse) {
        ff_fft_permute_c(&s->fft, (FFTComplex*)data, tmp_buf);
        ff_fft_calc_c(&s->fft, (FFTComplex*)data);
    }
    /* i=0 is a special case because of packing, the DC term is real, so we
       are going to throw the N/2 term (also real) in with it. */
    ev.re = data[0];
    data[0] = ev.re+data[1];
    data[1] = ev.re-data[1];
    for (i = 1; i < (n>>2); i++) {
        i1 = 2*i;
        i2 = n-i1;
        /* Separate even and odd FFTs */
        ev.re =  k1*(data[i1  ]+data[i2  ]);
        od.im = -k2*(data[i1  ]-data[i2  ]);
        ev.im =  k1*(data[i1+1]-data[i2+1]);
        od.re =  k2*(data[i1+1]+data[i2+1]);
        /* Apply twiddle factors to the odd FFT and add to the even FFT */
        data[i1  ] =  ev.re + od.re*tcos[i] - od.im*tsin[i];
        data[i1+1] =  ev.im + od.im*tcos[i] + od.re*tsin[i];
        data[i2  ] =  ev.re - od.re*tcos[i] + od.im*tsin[i];
        data[i2+1] = -ev.im + od.im*tcos[i] + od.re*tsin[i];
    }
    data[2*i+1]=s->sign_convention*data[2*i+1];
    if (s->inverse) {
        data[0] *= k1;
        data[1] *= k1;
        ff_fft_permute_c(&s->fft, (FFTComplex*)data, tmp_buf);
        ff_fft_calc_c(&s->fft, (FFTComplex*)data);
    }
}

#include <xmmintrin.h>
#define PM128(x) (*(__m128*)(x))
static _MM_ALIGN16 const unsigned long PCS_RNRN[4] = {0x00000000, 0x80000000, 0x00000000, 0x80000000};
static _MM_ALIGN16 const unsigned long PCS_NRNR[4] = {0x80000000, 0x00000000, 0x80000000, 0x00000000};
static _MM_ALIGN16 const float PPPP[4] = {+0.5f, +0.5f, +0.5f, +0.5f};
static _MM_ALIGN16 const float PMPM[4] = {+0.5f, -0.5f, +0.5f, -0.5f};
static _MM_ALIGN16 const float MPMP[4] = {-0.5f, +0.5f, -0.5f, +0.5f};

#pragma warning(disable:4700)
void ff_rdft_calc_sse(RDFTContext* s, FFTSample* data, FFTComplex *tmp_buf)
{
    int i, i1, i2;
    FFTComplex ev, od;
    //FFTComplex e2, o2;
    const int n = 1 << s->nbits;
    const float k1 = 0.5f;
    //const float k2 = (float)(0.5 - s->inverse);
    const float k2 = s->inverse? -0.5f : 0.5f;
    const FFTSample *tcos = s->tcos;
    const FFTSample *tsin = s->tsin;
    int n4 = n>>2;
    __m128 XMM0, XMM1, XMM2, XMM3, XMM4;
    const __m128 kk2 = s->inverse ? PM128(PMPM): PM128(MPMP);
    const __m128 rnrn = PM128(PCS_RNRN);

#ifdef _DEBUG
    XMM3 = XMM4 = _mm_setzero_ps();
#endif

    if (!s->inverse) {
        ff_fft_permute_sse(&s->fft, (FFTComplex*)data, tmp_buf);
        ff_fft_calc_sse(&s->fft, (FFTComplex*)data);
    }
    /* i=0 is a special case because of packing, the DC term is real, so we
       are going to throw the N/2 term (also real) in with it. */
    ev.re = data[0];
    data[0] = ev.re+data[1];
    data[1] = ev.re-data[1];
    /* Separate even and odd FFTs */
    ev.re =  k1*(data[2   ]+data[n-2 ]);
    od.im = -k2*(data[2   ]-data[n-2 ]);
    ev.im =  k1*(data[3   ]-data[n-1 ]);
    od.re =  k2*(data[3   ]+data[n-1 ]);
    /* Apply twiddle factors to the odd FFT and add to the even FFT */
    data[2   ] =  ev.re + od.re*tcos[1] - od.im*tsin[1];
    data[3   ] =  ev.im + od.im*tcos[1] + od.re*tsin[1];
    data[n-2 ] =  ev.re - od.re*tcos[1] + od.im*tsin[1];
    data[n-1 ] = -ev.im + od.im*tcos[1] + od.re*tsin[1];
#if 0
    for (i = 2; i < (n>>2); i++) {
        i1 = 2*i;
        i2 = n-i1;
        /* Separate even and odd FFTs */
        ev.re =  k1*(data[i1  ]+data[i2  ]);
        od.im = -k2*(data[i1  ]-data[i2  ]);
        ev.im =  k1*(data[i1+1]-data[i2+1]);
        od.re =  k2*(data[i1+1]+data[i2+1]);
        /* Apply twiddle factors to the odd FFT and add to the even FFT */
        data[i1  ] =  ev.re + od.re*tcos[i] - od.im*tsin[i];
        data[i1+1] =  ev.im + od.im*tcos[i] + od.re*tsin[i];
        data[i2  ] =  ev.re - od.re*tcos[i] + od.im*tsin[i];
        data[i2+1] = -ev.im + od.im*tcos[i] + od.re*tsin[i];
    }
#elif 0
    for (i = 2; i < (n>>2); i+=2) {
        i1 = 2*i;
        i2 = n-i1;
        /* Separate even and odd FFTs */
        ev.re =  k1*(data[i1  ]+data[i2  ]);
        ev.im =  k1*(data[i1+1]-data[i2+1]);
        e2.re =  k1*(data[i1+2]+data[i2-2]);
        e2.im =  k1*(data[i1+3]-data[i2-1]);

        od.im = -k2*(data[i1  ]-data[i2  ]);
        od.re =  k2*(data[i1+1]+data[i2+1]);
        o2.im = -k2*(data[i1+2]-data[i2-2]);
        o2.re =  k2*(data[i1+3]+data[i2-1]);

        /* Apply twiddle factors to the odd FFT and add to the even FFT */
        data[i1  ] =  ev.re + od.re*tcos[i  ] - od.im*tsin[i  ];
        data[i1+1] =  ev.im + od.im*tcos[i  ] + od.re*tsin[i  ];
        data[i1+2] =  e2.re + o2.re*tcos[i+1] - o2.im*tsin[i+1];
        data[i1+3] =  e2.im + o2.im*tcos[i+1] + o2.re*tsin[i+1];

        data[i2-2] =  e2.re - o2.re*tcos[i+1] + o2.im*tsin[i+1];
        data[i2-1] = -e2.im + o2.im*tcos[i+1] + o2.re*tsin[i+1];
        data[i2  ] =  ev.re - od.re*tcos[i  ] + od.im*tsin[i  ];
        data[i2+1] = -ev.im + od.im*tcos[i  ] + od.re*tsin[i  ];

    }
#else
    // Intrinsics by lvqcl
    for (i=2,i1=4,i2=n-4; i < n4; i+=2,i1+=4,i2-=4)
    {
        /* Separate even and odd FFTs */
        XMM0 = _mm_load_ps(data+i1); // data[i1]
        XMM1 = _mm_loadu_ps((data+i2-2));
        XMM1 = _mm_shuffle_ps(XMM1, XMM1, _MM_SHUFFLE(1, 0, 3, 2)); // data[i2]
        XMM1 = _mm_xor_ps(XMM1, rnrn); // +-data[i2]
        XMM2 = XMM0; //data[i1]

        XMM0 = _mm_add_ps(XMM0, XMM1); // data[i1] +- data[i2]
        XMM2 = _mm_sub_ps(XMM2, XMM1); // data[i1] -+ data[i2]

        XMM0 = _mm_mul_ps(XMM0, PM128(PPPP));                     //  ev.re,  ev.im,  e2.re,  e2.im
        XMM2 = _mm_mul_ps(XMM2, kk2);// PMPM or MPMP              //  od.im,  od.re,  o2.im,  o2.re
        XMM1 = XMM2;
        XMM1 = _mm_shuffle_ps(XMM1, XMM1, _MM_SHUFFLE(2,3,0,1));  //  od.re, od.im, o2.re, o2.im

        /* Apply twiddle factors to the odd FFT and add to the even FFT */
        XMM3 = _mm_loadl_pi(XMM3, (__m64*)(tcos+i));
        XMM4 = _mm_loadl_pi(XMM4, (__m64*)(tsin+i));
#pragma warning(default:4700)
        XMM3 = _mm_unpacklo_ps(XMM3, XMM3); // tcos[i], tcos[i], tcos[i+1], tcos[i+1]
        XMM4 = _mm_unpacklo_ps(XMM4, XMM4); // tsin[i], tsin[i], tsin[i+1], tsin[i+1]

        XMM1 = _mm_mul_ps(XMM1, XMM3); // od*cos
        XMM2 = _mm_mul_ps(XMM2, XMM4); // od*sin

        XMM4 = XMM0; // ev

        XMM0 = _mm_add_ps(XMM0, XMM1); // ev + od*cos
        XMM4 = _mm_sub_ps(XMM4, XMM1); // ev - od*cos

        XMM4 = _mm_xor_ps(XMM4, rnrn); // +-ev -+ od*cos
        XMM4 = _mm_add_ps(XMM4, XMM2); // +-ev -+ od*cos + od*sin
        XMM4 = _mm_shuffle_ps(XMM4, XMM4, _MM_SHUFFLE(1, 0, 3, 2));
        _mm_storeu_ps((data+i2-2), XMM4);

        XMM2 = _mm_xor_ps(XMM2, rnrn); // +-od*sin
        XMM0 = _mm_sub_ps(XMM0, XMM2); // ev + od*cos +- od*sin

        _mm_store_ps(data+i1, XMM0);
    }
#endif
    data[2*i+1] *= s->sign_convention;
    if (s->inverse) {
        data[0] *= k1;
        data[1] *= k1;
        ff_fft_permute_sse(&s->fft, (FFTComplex*)data, tmp_buf);
        ff_fft_calc_sse(&s->fft, (FFTComplex*)data);
    }
}


int ff_rdft_init(RDFTContext *s, int nbits, enum RDFTransformType trans, int sse)
{
    int n = 1 << nbits;
    int i;
    const double theta = (trans == DFT_R2C || trans == DFT_C2R ? -1 : 1)*2*M_PI/n;

    if (nbits < 4 || nbits > 16)
        return -1;

    if (ff_fft_init(&s->fft, nbits-1, trans == IDFT_C2R || trans == IDFT_R2C, sse) < 0)
        return -1;

    s->nbits           = nbits;
    s->inverse         = trans == IDFT_C2R || trans == DFT_C2R;
    s->sign_convention = trans == IDFT_R2C || trans == DFT_C2R ? 1 : -1;

    ff_init_ff_cos_tabs(nbits);
    s->tcos = ff_cos_tabs[nbits];
    s->tsin = ff_sin_tabs[nbits]+(trans == DFT_R2C || trans == DFT_C2R)*(n>>2);
#if !CONFIG_HARDCODED_TABLES
    if(s->tsin[1]==0.0) {
        for (i = 0; i < (n>>2); i++) {
            s->tsin[i] = (FFTSample)sin(i*theta);
        }
    }
#endif

    return 0;
}

void ff_rdft_end(RDFTContext *s)
{
    ff_fft_end(&s->fft);
}
