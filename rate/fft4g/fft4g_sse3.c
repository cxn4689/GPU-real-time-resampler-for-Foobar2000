/* Copyright(C) Takuya OOURA
 * email: ooura@mmm.t.u-tokyo.ac.jp
 * download: http://momonga.t.u-tokyo.ac.jp/~ooura/fft.html
 * You may use, copy, modify this code for any purpose and 
 * without fee. You may distribute this ORIGINAL package. 
 *
 * + some changes from SoX project
 * + SSE accelerated functions: Copyright (C) lvqcl
*/

#include "fft4g.h"

#include <pmmintrin.h>
#define _mm_lddqu_pd(x) (_mm_castsi128_pd(_mm_lddqu_si128((__m128i*)(x))))
#define _mm_shuffle_lh(x, y) (_mm_shuffle_pd((x), (y), _MM_SHUFFLE2(0, 1)))

#define _mm_neg_lo(x) (_mm_xor_pd((x), _mm_set_sd(-0.0)))
#define _mm_neg_hi(x) (_mm_xor_pd((x), _mm_set_pd(-0.0, 0.0)))
#define _mm_neg_lo2(x) (_mm_xor_pd((x), (neg_lo)))
#define _mm_neg_hi2(x) (_mm_xor_pd((x), (neg_hi)))

static __inline __m128d ZMUL(__m128d ai, __m128d b)
{
    __m128d ar;

    ar = _mm_movedup_pd(ai);      /* ar = [a.r a.r] */
    ai = _mm_unpackhi_pd(ai, ai); /* ai = [a.i a.i] */
    ar = _mm_mul_pd(ar, b);       /* ar = [a.r*b.r a.r*b.i] */
    b  = _mm_shuffle_lh(b, b);    /* b = [b.i b.r] */
    ai = _mm_mul_pd(ai, b);       /* ai = [a.i*b.i a.i*b.r] */
    return _mm_addsub_pd(ar, ai); /* [a.r*b.r-a.i*b.i a.r*b.i+a.i*b.r] */
}

static __inline __m128d Z111_hi(__m128d ai, __m128d b, __m128d neg_hi)
{
    __m128d ar;

    ar = _mm_movedup_pd(ai);      /* ar = [a.r a.r] */
    ai = _mm_unpackhi_pd(ai, ai); /* ai = [a.i a.i] */
    ar = _mm_mul_pd(ar, b);       /* ar = [a.r*b.r a.r*b.i] */
    b  = _mm_shuffle_lh(b, b);    /* b = [b.i b.r] */
    ai = _mm_mul_pd(ai, b);       /* ai = [a.i*b.i a.i*b.r] */
    ai = _mm_neg_hi2(ai);          /* ai = [a.i*b.i -a.i*b.r] */
    return _mm_add_pd(ar, ai); /* [a.r*b.r+a.i*b.i a.r*b.i-a.i*b.r] */
    //ai = _mm_neg_lo2(ai);          /* ai = [-a.i*b.i a.i*b.r] */
    //return _mm_sub_pd(ar, ai); /* [a.r*b.r+a.i*b.i a.r*b.i-a.i*b.r] */
    //return _mm_addsub_pd(ar, ai); /* [a.r*b.r-a.i*b.i a.r*b.i+a.i*b.r] */
}

static __inline void makebr(int n, int *ip);

static __inline void bitrv2(int n, int const *ip, double *a);
static __inline void cftbsub(int n, double *a, double const *w);
static __inline void cftfsub(int n, double *a, double const *w);
static __inline void rftbsub(int n, double *a, int nc, double const *c);
static __inline void rftfsub(int n, double *a, int nc, double const *c);

static __inline void cft1st(int n, double *a, double const *w);
static __inline void cftmdl(int n, int l, double *a, double const *w);


void lsx_rdft_SSE3(int isgn, double *a, FFTcontext z)
{
    double xi;
    int n = z.n;
    int nw = n >> 2;

    if (isgn >= 0)
    {
        if (n > 4)
        {
            bitrv2(n, z.br, a);
            cftfsub(n, a, z.sc);
            rftfsub(n, a, nw, z.sc + nw);
        }
        else if (n == 4)
        {
            cftfsub(n, a, z.sc);
        }
        xi = a[0] - a[1];
        a[0] += a[1];
        a[1] = xi;
    }
    else
    {
        a[1] = 0.5 * (a[0] - a[1]);
        a[0] -= a[1];
        if (n > 4)
        {
            rftbsub(n, a, nw, z.sc + nw);
            bitrv2(n, z.br, a);
            cftbsub(n, a, z.sc);
        }
        else if (n == 4)
        {
            cftfsub(n, a, z.sc);
        }
    }
}


/* -------- child routines -------- */


static __inline void bitrv2(int n, int const *ip, double *a)
{
    int j, j1, k, k1, l, m, m2;
    //double xr, xi, yr, yi;
    __m128d x1, y1, x2, y2;

    l = n; m = 1;
    while ((m << 3) < l) { l >>= 1; m <<= 1; }

    m2 = 2 * m;
    if ((m << 3) == l) {
        for (k = 0; k < m; k++) {
            for (j = 0; j < k; j++) {
                j1 = 2 * j + ip[k];
                k1 = 2 * k + ip[j];

                //xr = a[j1];
                //xi = a[j1 + 1];
                //yr = a[k1];
                //yi = a[k1 + 1];
                x1 = _mm_load_pd(a+j1);
                y1 = _mm_load_pd(a+k1);
                //a[j1] = yr;
                //a[j1 + 1] = yi;
                //a[k1] = xr;
                //a[k1 + 1] = xi;
                _mm_store_pd(a+j1, y1);
                _mm_store_pd(a+k1, x1);

                j1 += m2;
                k1 += 2 * m2;

                //xr = a[j1];
                //xi = a[j1 + 1];
                //yr = a[k1];
                //yi = a[k1 + 1];
                x2 = _mm_load_pd(a+j1);
                y2 = _mm_load_pd(a+k1);
                //a[j1] = yr;
                //a[j1 + 1] = yi;
                //a[k1] = xr;
                //a[k1 + 1] = xi;
                _mm_store_pd(a+j1, y2);
                _mm_store_pd(a+k1, x2);

                j1 += m2;
                k1 -= m2;

                //xr = a[j1];
                //xi = a[j1 + 1];
                //yr = a[k1];
                //yi = a[k1 + 1];
                x1 = _mm_load_pd(a+j1);
                y1 = _mm_load_pd(a+k1);
                //a[j1] = yr;
                //a[j1 + 1] = yi;
                //a[k1] = xr;
                //a[k1 + 1] = xi;
                _mm_store_pd(a+j1, y1);
                _mm_store_pd(a+k1, x1);

                j1 += m2;
                k1 += 2 * m2;

                //xr = a[j1];
                //xi = a[j1 + 1];
                //yr = a[k1];
                //yi = a[k1 + 1];
                x2 = _mm_load_pd(a+j1);
                y2 = _mm_load_pd(a+k1);
                //a[j1] = yr;
                //a[j1 + 1] = yi;
                //a[k1] = xr;
                //a[k1 + 1] = xi;
                _mm_store_pd(a+j1, y2);
                _mm_store_pd(a+k1, x2);
            }
            j1 = 2 * k + m2 + ip[k];
            k1 = j1 + m2;

            //xr = a[j1];
            //xi = a[j1 + 1];
            //yr = a[k1];
            //yi = a[k1 + 1];
            x1 = _mm_load_pd(a+j1);
            y1 = _mm_load_pd(a+k1);
            //a[j1] = yr;
            //a[j1 + 1] = yi;
            //a[k1] = xr;
            //a[k1 + 1] = xi;
            _mm_store_pd(a+j1, y1);
            _mm_store_pd(a+k1, x1);
        }
    } else {
        for (k = 1; k < m; k++) {
            for (j = 0; j < k; j++) {
                j1 = 2 * j + ip[k];
                k1 = 2 * k + ip[j];

                //xr = a[j1];
                //xi = a[j1 + 1];
                //yr = a[k1];
                //yi = a[k1 + 1];
                x1 = _mm_load_pd(a+j1);
                y1 = _mm_load_pd(a+k1);
                //a[j1] = yr;
                //a[j1 + 1] = yi;
                //a[k1] = xr;
                //a[k1 + 1] = xi;
                _mm_store_pd(a+j1, y1);
                _mm_store_pd(a+k1, x1);

                j1 += m2;
                k1 += m2;

                //xr = a[j1];
                //xi = a[j1 + 1];
                //yr = a[k1];
                //yi = a[k1 + 1];
                x2 = _mm_load_pd(a+j1);
                y2 = _mm_load_pd(a+k1);
                //a[j1] = yr;
                //a[j1 + 1] = yi;
                //a[k1] = xr;
                //a[k1 + 1] = xi;
                _mm_store_pd(a+j1, y2);
                _mm_store_pd(a+k1, x2);
            }
        }
    }
}


static __inline void cftfsub(int n, double *a, double const *w)
{
    int j, j1, j2, j3, l;
    //double x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i;
    __m128d x0, x1, x2, x3, t1, t2;
    __m128d neg_lo = _mm_set_sd(-0.0);

    l = 2;
    if (n > 8) {
        cft1st(n, a, w);
        l = 8;
        while ((l << 2) < n) {
            cftmdl(n, l, a, w);
            l <<= 2;
        }
    }
    if ((l << 2) == n) {
        for (j = 0; j < l; j += 2) {
            j1 = j + l;
            j2 = j1 + l;
            j3 = j2 + l;

            //x0r = a[j    ] + a[j1    ];
            //x0i = a[j + 1] + a[j1 + 1];
            //x1r = a[j    ] - a[j1    ];
            //x1i = a[j + 1] - a[j1 + 1];
            x0 = _mm_load_pd(a+j );
            t1  = _mm_load_pd(a+j1);
            x1 = _mm_sub_pd(x0, t1);
            x0 = _mm_add_pd(x0, t1);
            //x2r = a[j2    ] + a[j3    ];
            //x2i = a[j2 + 1] + a[j3 + 1];
            //x3r = a[j2    ] - a[j3    ];
            //x3i = a[j2 + 1] - a[j3 + 1];
            x2 = _mm_load_pd(a+j2);
            t2 = _mm_load_pd(a+j3);
            x3 = _mm_sub_pd(x2, t2);
            x2 = _mm_add_pd(x2, t2);
            //a[j    ] = x0r + x2r;
            //a[j + 1] = x0i + x2i;
            _mm_store_pd(a+j , _mm_add_pd(x0, x2));
            //a[j2    ] = x0r - x2r;
            //a[j2 + 1] = x0i - x2i;
            _mm_store_pd(a+j2, _mm_sub_pd(x0, x2));
            //a[j1    ] = x1r - x3i;
            //a[j1 + 1] = x1i + x3r;
            x3 = _mm_neg_lo2(_mm_shuffle_lh(x3, x3));
            _mm_store_pd(a+j1, _mm_add_pd(x1, x3));
            //a[j3    ] = x1r + x3i;
            //a[j3 + 1] = x1i - x3r;
            _mm_store_pd(a+j3, _mm_sub_pd(x1, x3));
        }
    } else {
        for (j = 0; j < l; j += 2) {
            j1 = j + l;

            //x0r = a[j    ] - a[j1    ];
            //x0i = a[j + 1] - a[j1 + 1];
            x0 = _mm_load_pd(a+j );
            x1 = _mm_load_pd(a+j1);
            //a[j    ] += a[j1    ];
            //a[j + 1] += a[j1 + 1];
            _mm_store_pd(a+j, _mm_add_pd(x0, x1));
            //a[j1    ] = x0r;
            //a[j1 + 1] = x0i;
            _mm_store_pd(a+j1, _mm_sub_pd(x0, x1));
            //    vv
            //a[j     ] = a[j    ] + a[j1    ];
            //a[j  + 1] = a[j + 1] + a[j1 + 1];
            //a[j1    ] = a[j    ] - a[j1    ];
            //a[j1 + 1] = a[j + 1] - a[j1 + 1];
        }
    }
}


static __inline void cftbsub(int n, double *a, double const *w)
{
    int j, j1, j2, j3, l;
    //double x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i;
    __m128d x0, x1, x2, x3, t1, t2;
    __m128d neg_hi = _mm_set_pd(-0.0, 0.0);

    l = 2;
    if (n > 8) {
        cft1st(n, a, w);
        l = 8;
        while ((l << 2) < n) {
            cftmdl(n, l, a, w);
            l <<= 2;
        }
    }
    if ((l << 2) == n) {
        for (j = 0; j < l; j += 2) {
            j1 = j + l;
            j2 = j1 + l;
            j3 = j2 + l;

            //x0r =  (a[j    ] + a[j1    ]);
            //x0i = -(a[j + 1] + a[j1 + 1]);
            //x1r =  (a[j    ] - a[j1    ]);
            //x1i = -(a[j + 1] - a[j1 + 1]);
            x0 = _mm_load_pd(a+j );
            t1 = _mm_load_pd(a+j1);
            x1 = _mm_sub_pd(x0, t1);
            x0 = _mm_add_pd(x0, t1);
            x1 = _mm_neg_hi2(x1);
            x0 = _mm_neg_hi2(x0);


            //x2r = a[j2    ] + a[j3    ];
            //x2i = a[j2 + 1] + a[j3 + 1];
            //x3r = a[j2    ] - a[j3    ];
            //x3i = a[j2 + 1] - a[j3 + 1];
            x2 = _mm_load_pd(a+j2);
            t2 = _mm_load_pd(a+j3);
            x3 = _mm_sub_pd(x2, t2);
            x2 = _mm_add_pd(x2, t2);

            //a[j     ] = x0r + x2r;
            //a[j  + 1] = x0i - x2i;
            t1 = _mm_add_pd(x0, _mm_neg_hi2(x2));
            _mm_store_pd(a+j, t1);

            //a[j2    ] = x0r - x2r;
            //a[j2 + 1] = x0i + x2i;
            x0 = _mm_addsub_pd(x0, x2);
            _mm_store_pd(a+j2, x0);

            //a[j1    ] = x1r - x3i;
            //a[j1 + 1] = x1i - x3r;
            x3 = _mm_shuffle_lh(x3, x3);
            _mm_store_pd(a+j1, _mm_sub_pd(x1, x3));

            //a[j3    ] = x1r + x3i;
            //a[j3 + 1] = x1i + x3r;
            _mm_store_pd(a+j3, _mm_add_pd(x1, x3));
        }
    } else {
        for (j = 0; j < l; j += 2) {
            j1 = j + l;

            //x0r =  (a[j    ] - a[j1    ]);
            //x0i = -(a[j + 1] - a[j1 + 1]);
            x0 = _mm_load_pd(a+j );
            t1 = _mm_load_pd(a+j1);
            x1 = _mm_sub_pd(x0, t1);
            x0 = _mm_add_pd(x0, t1);

            //a[j    ] =  (a[j    ] + a[j1    ]);
            //a[j + 1] = -(a[j + 1] + a[j1 + 1]);
            _mm_store_pd(a+j, _mm_neg_hi2(x0));

            //a[j1    ] =  (a[j    ] - a[j1    ]);
            //a[j1 + 1] = -(a[j + 1] - a[j1 + 1]);
            _mm_store_pd(a+j1, _mm_neg_hi2(x1));
        }
    }
}


static __inline void rftfsub(int n, double *a, int nc, double const *c)
{
    int j, k, kk, ks, m;
    //double wkr, wki, xr, xi, yr, yi;
    __m128d wk, aj, ak, x, y;
    __m128d neg_hi = _mm_set_pd(-0.0, 0.0);

    m = n >> 1;
    ks = 2 * nc / m; //if nc == n/4  =>  ks == 1;
    kk = 0;
    for (j = 2; j < m; j += 2) {
        k = n - j;
        kk += ks;

        //wkr = 0.5 - c[nc - kk];
        //wki = c[kk];
        wk = _mm_set_sd(0.5 - c[nc - kk]);
        wk = _mm_loadh_pd(wk, c+kk);
        //xr = a[j    ] - a[k    ];
        //xi = a[j + 1] + a[k + 1];
        aj = _mm_load_pd(a+j);
        ak = _mm_load_pd(a+k);
        x = _mm_addsub_pd(aj, ak);

        //yr = wkr * xr - wki * xi;
        //yi = wkr * xi + wki * xr;
        y = ZMUL(wk, x);

        //a[j    ] = a[j    ] - yr;
        //a[j + 1] = a[j + 1] - yi;
        _mm_store_pd(a+j, _mm_sub_pd(aj, y));
        //a[k    ] = a[k    ] + yr;
        //a[k + 1] = a[k + 1] - yi;
        ak = _mm_add_pd(ak, _mm_neg_hi2(y));
        _mm_store_pd(a+k, ak);
    }
}


static __inline void rftbsub(int n, double *a, int nc, double const *c)
{
    int j, k, kk, ks, m;
    //double wkr, wki, xr, xi, yr, yi;
    __m128d wk, aj, ak, x, y;
    __m128d neg_hi = _mm_set_pd(-0.0, 0.0);

    a[1] = -a[1];
    m = n >> 1;
    ks = 2 * nc / m; //if nc == n/4  =>  ks == 1;
    kk = 0;
    for (j = 2; j < m; j += 2) {
        k = n - j;
        kk += ks;

        //wkr = 0.5 - c[nc - kk];
        //wki = c[kk];
        wk = _mm_set_sd(0.5 - c[nc - kk]);
        wk = _mm_loadh_pd(wk, c+kk);
        //xr = a[j    ] - a[k    ];
        //xi = a[j + 1] + a[k + 1];
        aj = _mm_load_pd(a+j);
        ak = _mm_load_pd(a+k);
        ak = _mm_neg_hi2(ak);
        x = _mm_sub_pd(aj, ak);

        //yr = wkr * xr + wki * xi;
        //yi = wkr * xi - wki * xr;
        y = Z111_hi(wk, x, neg_hi);

        //a[j    ] =  (a[j    ] - yr);
        //a[j + 1] = -(a[j + 1] - yi);
        aj = _mm_sub_pd(aj, y);
        _mm_store_pd(a+j, _mm_neg_hi2(aj));
        //a[k    ] =  a[k    ] + yr;
        //a[k + 1] = -a[k + 1] + yi;
        _mm_store_pd(a+k, _mm_add_pd(ak, y));
    }
    a[m + 1] = -a[m + 1];
}


static __inline void cft1st(int n, double *a, double const *w)
{
    int j, k1, k2;
    //double wk1r, wk1i, wk2r, wk2i, wk3r, wk3i;
    //double x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i;
    __m128d x0, x1, x2, x3, t1, t2;

    //x0r = a[0] + a[2];
    //x0i = a[1] + a[3];
    //x1r = a[0] - a[2];
    //x1i = a[1] - a[3];
    x0 = _mm_load_pd(a+0);
    t1  = _mm_load_pd(a+2);
    x1 = _mm_sub_pd(x0, t1);
    x0 = _mm_add_pd(x0, t1);

    //x2r = a[4] + a[6];
    //x2i = a[5] + a[7];
    //x3r = a[4] - a[6];
    //x3i = a[5] - a[7];
    x2 = _mm_load_pd(a+4);
    t2 = _mm_load_pd(a+6);
    x3 = _mm_sub_pd(x2, t2);
    x2 = _mm_add_pd(x2, t2);

    //a[0] = x0r + x2r;
    //a[1] = x0i + x2i;
    _mm_store_pd(a+0, _mm_add_pd(x0, x2));
    //a[4] = x0r - x2r;
    //a[5] = x0i - x2i;
    _mm_store_pd(a+4, _mm_sub_pd(x0, x2));

    t1 = _mm_neg_lo(_mm_shuffle_lh(x3, x3));
    //a[2] = x1r - x3i;
    //a[3] = x1i + x3r;
    _mm_store_pd(a+2, _mm_add_pd(x1, t1));
    //a[6] = x1r + x3i;
    //a[7] = x1i - x3r;
    _mm_store_pd(a+6, _mm_sub_pd(x1, t1));

    //x0r = a[8] + a[10];
    //x0i = a[9] + a[11];
    //x1r = a[8] - a[10];
    //x1i = a[9] - a[11];
    x0 = _mm_load_pd(a+8);
    t1  = _mm_load_pd(a+10);
    x1 = _mm_sub_pd(x0, t1);
    x0 = _mm_add_pd(x0, t1);

    //x2r = a[12] + a[14];
    //x2i = a[13] + a[15];
    //x3r = a[12] - a[14];
    //x3i = a[13] - a[15];
    x2 = _mm_load_pd(a+12);
    t2 = _mm_load_pd(a+14);
    x3 = _mm_sub_pd(x2, t2);
    x2 = _mm_add_pd(x2, t2);

    //a[8] = x0r + x2r;
    //a[9] = x0i + x2i;
    _mm_store_pd(a+8, _mm_add_pd(x0, x2));
    //a[12] = x2i - x0i; = - (x0i - x2i)
    //a[13] = x0r - x2r; = + (x0r - x2r)
    x0 = _mm_sub_pd(x0, x2);
    _mm_store_pd(a+12, _mm_neg_lo(_mm_shuffle_lh(x0, x0)));

    t1 = _mm_shuffle_lh(x3, x3);
    //x0r = x1r - x3i;
    //x0i = x1i + x3r;
    x0 = _mm_addsub_pd(x1, t1);
    //x2r = x3i + x1r;
    //x2i = x3r - x1i;
    x2 = _mm_add_pd(t1, _mm_neg_hi(x1));

    //wk1r = w[2];
    t2 = _mm_load1_pd(w+2);

    //a[10] = wk1r * (x0r - x0i);
    //a[11] = wk1r * (x0r + x0i);
    x0 = _mm_hadd_pd(_mm_neg_hi(x0), x0);
    _mm_store_pd(a+10, _mm_mul_pd(t2, x0));

    //a[14] = wk1r * (x2i - x2r); = wk1r * (-x2r + x2i)
    //a[15] = wk1r * (x2i + x2r); = wk1r * (+x2r + x2i)
    x2 = _mm_hadd_pd(_mm_neg_lo(x2), x2);
    _mm_store_pd(a+14, _mm_mul_pd(t2, x2));

    k1 = 0;
    for (j = 16; j < n; j += 16) {
        __m128d wk1, wk2, wk3;
        k1 += 2;
        k2 = 2 * k1;
        //wk2r = w[k1];
        //wk2i = w[k1 + 1];
        wk2 = _mm_load_pd(w+k1);
        //wk1r = w[k2];
        //wk1i = w[k2 + 1];
        wk1 = _mm_load_pd(w+k2);

        //wk3r =   wk1r - 2 * wk2i * wk1i;
        //wk3i = - wk1i + 2 * wk2i * wk1r;
        t2 = _mm_unpackhi_pd(wk2, wk2);
        t2 = _mm_mul_pd(t2, wk1);
        t2 = _mm_add_pd(t2, t2);
        //wk3r = + wk1r - t_i; = + (wk1r - t_i)
        //wk3i = - wk1i + t_r; = - (wk1i - t_r)
        wk3 = _mm_sub_pd(wk1, _mm_shuffle_lh(t2, t2));
        wk3 = _mm_neg_hi(wk3);

        //x0r = a[j    ] + a[j + 2];
        //x0i = a[j + 1] + a[j + 3];
        //x1r = a[j    ] - a[j + 2];
        //x1i = a[j + 1] - a[j + 3];
        x0 = _mm_load_pd(a+j+0);
        t1  = _mm_load_pd(a+j+2);
        x1 = _mm_sub_pd(x0, t1);
        x0 = _mm_add_pd(x0, t1);

        //x2r = a[j + 4] + a[j + 6];
        //x2i = a[j + 5] + a[j + 7];
        //x3r = a[j + 4] - a[j + 6];
        //x3i = a[j + 5] - a[j + 7];
        x2 = _mm_load_pd(a+j+4);
        t2 = _mm_load_pd(a+j+6);
        x3 = _mm_sub_pd(x2, t2);
        x2 = _mm_add_pd(x2, t2);

        //a[j]     = x0r + x2r;
        //a[j + 1] = x0i + x2i;
        _mm_store_pd(a+j, _mm_add_pd(x0, x2));

        //x0r -= x2r;
        //x0i -= x2i;
        x0 = _mm_sub_pd(x0, x2);
        //a[j + 4] = wk2r * x0r - wk2i * x0i;
        //a[j + 5] = wk2r * x0i + wk2i * x0r;
        _mm_store_pd(a+j+4, ZMUL(wk2, x0));

        //x0r = x1r + x3i;
        //x0i = x1i - x3r;
        //x1r = x1r - x3i;
        //x1i = x1i + x3r;
        t1 = _mm_neg_lo(x3);
        t1 = _mm_shuffle_lh(t1, t1);
        x0 = _mm_add_pd(x1, t1);
        x1 = _mm_sub_pd(x1, t1);
        //a[j + 2] = wk1r * x1r - wk1i * x1i;
        //a[j + 3] = wk1r * x1i + wk1i * x1r;
        _mm_store_pd(a+j+2, ZMUL(wk1, x1));
        //a[j + 6] = wk3r * x0r - wk3i * x0i;
        //a[j + 7] = wk3r * x0i + wk3i * x0r;
        _mm_store_pd(a+j+6, ZMUL(wk3, x0));

        //wk1r = w[k2 + 2];
        //wk1i = w[k2 + 3];
        wk1 = _mm_load_pd(w+k2+2);

        //wk3r = + wk1r - 2 * wk2r * wk1i;
        //wk3i = - wk1i + 2 * wk2r * wk1r;
        t2 = _mm_unpacklo_pd(wk2, wk2);
        t2 = _mm_mul_pd(t2, wk1);
        t2 = _mm_add_pd(t2, t2);
        //wk3r = + wk1r - t_i; = + (wk1r - t_i)
        //wk3i = - wk1i + t_r; = - (wk1i - t_r)
        wk3 = _mm_sub_pd(wk1, _mm_shuffle_lh(t2, t2));
        wk3 = _mm_neg_hi(wk3);

        //x0r = a[j + 8] + a[j + 10];
        //x0i = a[j + 9] + a[j + 11];
        //x1r = a[j + 8] - a[j + 10];
        //x1i = a[j + 9] - a[j + 11];
        x0 = _mm_load_pd(a+j+8);
        t1  = _mm_load_pd(a+j+10);
        x1 = _mm_sub_pd(x0, t1);
        x0 = _mm_add_pd(x0, t1);

        //x2r = a[j + 12] + a[j + 14];
        //x2i = a[j + 13] + a[j + 15];
        //x3r = a[j + 12] - a[j + 14];
        //x3i = a[j + 13] - a[j + 15];
        x2 = _mm_load_pd(a+j+12);
        t2 = _mm_load_pd(a+j+14);
        x3 = _mm_sub_pd(x2, t2);
        x2 = _mm_add_pd(x2, t2);

        //a[j + 8] = x0r + x2r;
        //a[j + 9] = x0i + x2i;
        _mm_store_pd(a+j+8, _mm_add_pd(x0, x2));

        //x0r -= x2r;
        //x0i -= x2i;
        x0 = _mm_sub_pd(x0, x2);
        //a[j + 12] = -wk2i * x0r - wk2r * x0i;
        //a[j + 13] = -wk2i * x0i + wk2r * x0r;
        t1 = ZMUL(wk2, x0);
        t1 = _mm_shuffle_lh(t1, t1);
        t1 = _mm_neg_lo(t1);
        _mm_store_pd(a+j+12, t1);

        //x0r = x1r + x3i;
        //x0i = x1i - x3r;
        //x1r = x1r - x3i;
        //x1i = x1i + x3r;
        t2 = _mm_neg_lo(x3);
        t2 = _mm_shuffle_lh(t2, t2);
        x0 = _mm_add_pd(x1, t2);
        x1 = _mm_sub_pd(x1, t2);

        //a[j + 10] = wk1r * x1r - wk1i * x1i;
        //a[j + 11] = wk1r * x1i + wk1i * x1r;
        _mm_store_pd(a+j+10, ZMUL(wk1, x1));

        //a[j + 14] = wk3r * x0r - wk3i * x0i;
        //a[j + 15] = wk3r * x0i + wk3i * x0r;
        _mm_store_pd(a+j+14, ZMUL(wk3, x0));
    }
}


static __inline void cftmdl(int n, int l, double *a, double const *w)
{
    int j, j1, j2, j3, k, k1, k2, m, m2;
    //double wk1r, wk1i, wk2r, wk2i, wk3r, wk3i;
    //double x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i;
    __m128d x0, x1, x2, x3, wk1;

    m = l << 2;
    for (j = 0; j < l; j += 2)
    {
        __m128d t1, t2;
        j1 = j  + l;
        j2 = j1 + l;
        j3 = j2 + l;

        //x0r = a[j    ] + a[j1    ];
        //x0i = a[j + 1] + a[j1 + 1];
        //x1r = a[j    ] - a[j1    ];
        //x1i = a[j + 1] - a[j1 + 1];
        x0 = _mm_load_pd(a+j );
        t1  = _mm_load_pd(a+j1);
        x1 = _mm_sub_pd(x0, t1);
        x0 = _mm_add_pd(x0, t1);

        //x2r = a[j2    ] + a[j3    ];
        //x2i = a[j2 + 1] + a[j3 + 1];
        //x3r = a[j2    ] - a[j3    ];
        //x3i = a[j2 + 1] - a[j3 + 1];
        x2 = _mm_load_pd(a+j2);
        t2 = _mm_load_pd(a+j3);
        x3 = _mm_sub_pd(x2, t2);
        x2 = _mm_add_pd(x2, t2);

        //a[j    ] = x0r + x2r;
        //a[j + 1] = x0i + x2i;
        _mm_store_pd(a+j , _mm_add_pd(x0, x2));
        //a[j2    ] = x0r - x2r;
        //a[j2 + 1] = x0i - x2i;
        _mm_store_pd(a+j2, _mm_sub_pd(x0, x2));
        //a[j1    ] = x1r - x3i;
        //a[j1 + 1] = x1i + x3r;
        t1 = _mm_neg_lo(_mm_shuffle_lh(x3, x3));
        _mm_store_pd(a+j1, _mm_add_pd(x1, t1));
        //a[j3    ] = x1r + x3i;
        //a[j3 + 1] = x1i - x3r;
        _mm_store_pd(a+j3, _mm_sub_pd(x1, t1));
    }
    //wk1r = w[2];
    wk1 = _mm_load1_pd(w+2);

    for (j = m; j < l + m; j += 2)
    {
        __m128d t1, t2;
        j1 = j  + l;
        j2 = j1 + l;
        j3 = j2 + l;

        //x0r = a[j] + a[j1];
        //x0i = a[j + 1] + a[j1 + 1];
        //x1r = a[j] - a[j1];
        //x1i = a[j + 1] - a[j1 + 1];
        x0 = _mm_load_pd(a+j );
        t1  = _mm_load_pd(a+j1);
        x1 = _mm_sub_pd(x0, t1);
        x0 = _mm_add_pd(x0, t1);

        //x2r = a[j2] + a[j3];
        //x2i = a[j2 + 1] + a[j3 + 1];
        //x3r = a[j2] - a[j3];
        //x3i = a[j2 + 1] - a[j3 + 1];
        x2 = _mm_load_pd(a+j2);
        t2 = _mm_load_pd(a+j3);
        x3 = _mm_sub_pd(x2, t2);
        x2 = _mm_add_pd(x2, t2);

        //a[j    ] = x0r + x2r;
        //a[j + 1] = x0i + x2i;
        _mm_store_pd(a+j, _mm_add_pd(x0, x2));
        //a[j2    ] = x2i - x0i; = - (x0i - x2i)
        //a[j2 + 1] = x0r - x2r; = + (x0r - x2r)
        t1 = _mm_sub_pd(x0, x2);
        t1 = _mm_neg_lo(_mm_shuffle_lh(t1, t1));
        _mm_store_pd(a+j2, t1);

        //x0r = x1r - x3i;
        //x0i = x1i + x3r;
        x0 = _mm_addsub_pd(x1, _mm_shuffle_lh(x3, x3));
        //a[j1    ] = wk1r * (x0r - x0i);
        //a[j1 + 1] = wk1r * (x0r + x0i);
        x0 = _mm_hadd_pd(_mm_neg_hi(x0), x0);
        _mm_store_pd(a+j1, _mm_mul_pd(wk1, x0));

        //x0r = x3i + x1r;
        //x0i = x3r - x1i;
        x0 = _mm_addsub_pd(x3, _mm_shuffle_lh(x1, x1));
        x0 = _mm_shuffle_lh(x0, x0);
        //a[j3    ] = wk1r * (x0i - x0r); = wk1r * (- x0r + x0i)
        //a[j3 + 1] = wk1r * (x0i + x0r); = wk1r * (+ x0r + x0i)
        x0 = _mm_hadd_pd(_mm_neg_lo(x0), x0);
        _mm_store_pd(a+j3, _mm_mul_pd(wk1, x0));
    }

    k1 = 0;
    m2 = 2 * m;

    for (k = m2; k < n; k += m2)
    {
        __m128d wk2, wk3, t;
        k1 += 2;
        k2 = 2 * k1;

        //wk2r = w[k1];
        //wk2i = w[k1 + 1];
        wk2 = _mm_load_pd(w+k1);
        //wk1r = w[k2];
        //wk1i = w[k2 + 1];
        wk1 = _mm_load_pd(w+k2);

        //wk3r = wk1r - 2 * wk2i * wk1i; = + (wk1r - 2 * wk2i * wk1i)
        //wk3i = 2 * wk2i * wk1r - wk1i; = - (wk1i - 2 * wk2i * wk1r)
        t = _mm_unpackhi_pd(wk2, wk2);
        t = _mm_mul_pd(t, wk1);
        t = _mm_add_pd(t, t);
        //wk3r = + (wk1r - t_i)
        //wk3i = - (wk1i - t_r)
        wk3 = _mm_sub_pd(wk1, _mm_shuffle_lh(t, t));
        wk3 = _mm_neg_hi(wk3);

        for (j = k; j < l + k; j += 2)
        {
            j1 = j  + l;
            j2 = j1 + l;
            j3 = j2 + l;

            //x0r = a[j    ] + a[j1    ];
            //x0i = a[j + 1] + a[j1 + 1];
            //x1r = a[j    ] - a[j1    ];
            //x1i = a[j + 1] - a[j1 + 1];
            x0 = _mm_load_pd(a+j );
            t  = _mm_load_pd(a+j1);
            x1 = _mm_sub_pd(x0, t);
            x0 = _mm_add_pd(x0, t);
            //x2r = a[j2] + a[j3];
            //x2i = a[j2 + 1] + a[j3 + 1];
            //x3r = a[j2] - a[j3];
            //x3i = a[j2 + 1] - a[j3 + 1];
            x2 = _mm_load_pd(a+j2);
            t  = _mm_load_pd(a+j3);
            x3 = _mm_sub_pd(x2, t);
            x2 = _mm_add_pd(x2, t);

            //a[j    ] = x0r + x2r;
            //a[j + 1] = x0i + x2i;
            _mm_store_pd(a+j, _mm_add_pd(x0, x2));

            //x0r -= x2r;
            //x0i -= x2i;
            x0 = _mm_sub_pd(x0, x2);
            //a[j2    ] = wk2r * x0r - wk2i * x0i;
            //a[j2 + 1] = wk2r * x0i + wk2i * x0r;
            _mm_store_pd(a+j2, ZMUL(wk2, x0));

            //x0r = x1r + x3i;
            //x0i = x1i - x3r;
            //x1r = x1r - x3i;
            //x1i = x1i + x3r;
            t = _mm_neg_lo(x3);
            t = _mm_shuffle_lh(t, t);
            x0 = _mm_add_pd(x1, t);
            x1 = _mm_sub_pd(x1, t);

            //a[j1    ] = wk1r * x1r - wk1i * x1i;
            //a[j1 + 1] = wk1r * x1i + wk1i * x1r;
            _mm_store_pd(a+j1, ZMUL(wk1, x1));
            //a[j3    ] = wk3r * x0r - wk3i * x0i;
            //a[j3 + 1] = wk3r * x0i + wk3i * x0r;
            _mm_store_pd(a+j3, ZMUL(wk3, x0));
        }
        //wk1r = w[k2 + 2];
        //wk1i = w[k2 + 3];
        wk1 = _mm_load_pd(w+k2+2);

        //wk3r = + wk1r - 2 * wk2r * wk1i;
        //wk3i = - wk1i + 2 * wk2r * wk1r;
        t = _mm_unpacklo_pd(wk2, wk2);
        t = _mm_mul_pd(t, wk1);
        t = _mm_add_pd(t, t);
        //wk3r = + wk1r - t_i; = + (wk1r - t_i)
        //wk3i = - wk1i + t_r; = - (wk1i - t_r)
        wk3 = _mm_sub_pd(wk1, _mm_shuffle_lh(t, t));
        wk3 = _mm_neg_hi(wk3);

        for (j = k + m; j < l + (k + m); j += 2)
        {
            j1 = j  + l;
            j2 = j1 + l;
            j3 = j2 + l;

            //x0r = a[j    ] + a[j1    ];
            //x0i = a[j + 1] + a[j1 + 1];
            //x1r = a[j    ] - a[j1    ];
            //x1i = a[j + 1] - a[j1 + 1];
            x0 = _mm_load_pd(a+j );
            t  = _mm_load_pd(a+j1);
            x1 = _mm_sub_pd(x0, t);
            x0 = _mm_add_pd(x0, t);

            //x2r = a[j2    ] + a[j3    ];
            //x2i = a[j2 + 1] + a[j3 + 1];
            //x3r = a[j2    ] - a[j3    ];
            //x3i = a[j2 + 1] - a[j3 + 1];
            x2 = _mm_load_pd(a+j2);
            t  = _mm_load_pd(a+j3);
            x3 = _mm_sub_pd(x2, t);
            x2 = _mm_add_pd(x2, t);

            //a[j    ] = x0r + x2r;
            //a[j + 1] = x0i + x2i;
            _mm_store_pd(a+j , _mm_add_pd(x0, x2));

            //x0r -= x2r;
            //x0i -= x2i;
            x0 = _mm_sub_pd(x0, x2);
            //a[j2    ] = -wk2i * x0r - wk2r * x0i;
            //a[j2 + 1] = -wk2i * x0i + wk2r * x0r;
            t = ZMUL(wk2, x0);
            t = _mm_neg_lo(_mm_shuffle_lh(t, t));
            _mm_store_pd(a+j2, t);

            //x0r = x1r + x3i;
            //x0i = x1i - x3r;
            //x1r = x1r - x3i;
            //x1i = x1i + x3r;
            t = _mm_neg_lo(x3);
            t = _mm_shuffle_lh(t, t);
            x0 = _mm_add_pd(x1, t);
            x1 = _mm_sub_pd(x1, t);
            //a[j1    ] = wk1r * x1r - wk1i * x1i;
            //a[j1 + 1] = wk1r * x1i + wk1i * x1r;
            _mm_store_pd(a+j1, ZMUL(wk1, x1));
            //a[j3    ] = wk3r * x0r - wk3i * x0i;
            //a[j3 + 1] = wk3r * x0i + wk3i * x0r;
            _mm_store_pd(a+j3, ZMUL(wk3, x0));
        }
    }
}
