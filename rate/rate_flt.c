/* Effect: change sample rate     Copyright (c) 2008 robs@users.sourceforge.net
 * Copyright (C) 2008-11 lvqcl - some changes for thread safety & foobar2000
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

/* Inspired by, and builds upon some of the ideas presented in:
 * `The Quest For The Perfect Resampler' by Laurent De Soras;
 * http://ldesoras.free.fr/doc/articles/resampler-en.pdf */

#include <assert.h>
#include <string.h>

#include "sox_i.h"

#define  FIFO_SIZE_T int
#include "fifo.h"

typedef double raw_coef_t;
typedef float sox_sample_t;

typedef struct {
  int dft_length, num_taps, post_peak;
  sox_sample_t * coefs;

  FFTComplex * tmp_buf;
} dft_filter_t;

#if 0
static int lsx_lpf_num_taps(double att, double tr_bw, int k)
{                    /* TODO this could be cleaner, esp. for k != 0 */
  int n;
  if (att <= 80)
    n = (int)(.25 / M_PI * (att - 7.95) / (2.285 * tr_bw) + .5);
  else {
    double n160 = (.0425* att - 1.4) / tr_bw;   /* Half order for att = 160 */
    n = (int)(n160 * (16.556 / (att - 39.6) + .8625) + .5);  /* For att [80,160) */
  }
  if (k && (n&1)) ++n; // for SSE; we can also use "#ifdef SSE_"
  return k? 2 * n : 2 * (n + (n & 1)) + 1; /* =1 %4 (0 phase 1/2 band) */
}
#endif

#define  coef(coef_p, interp_order, fir_len, phase_num, coef_interp_num, fir_coef_num) coef_p[(fir_len) * ((interp_order) + 1) * (phase_num) + ((interp_order) + 1) * (fir_coef_num) + (interp_order - coef_interp_num)]

static sox_sample_t * prepare_coefs(raw_coef_t const * coefs, int num_coefs,
    int num_phases, int interp_order, int multiplier, int q)
{
  int i, j, length = num_coefs * num_phases;
  sox_sample_t * result = (sox_sample_t*)lsx_aligned_malloc(length * (interp_order + 1) * sizeof(*result) +RESERV);
  double fm1 = coefs[0], f1 = 0, f2 = 0;
  if (result == NULL) return NULL;

  for (i = num_coefs - 1; i >= 0; --i)
    for (j = num_phases - 1; j >= 0; --j) {
      double f0 = fm1, b = 0, c = 0, d = 0; /* = 0 to kill compiler warning */
      int pos = i * num_phases + j - 1;
      fm1 = (pos > 0 ? coefs[pos - 1] : 0) * multiplier;
      switch (interp_order) {
        case 1: b = f1 - f0; break;
        case 2: b = f1 - (.5 * (f2+f0) - f1) - f0; c = .5 * (f2+f0) - f1; break;
        case 3: c=.5*(f1+fm1)-f0;d=(1/6.)*(f2-f1+fm1-f0-4*c);b=f1-f0-d-c; break;
        default: if (interp_order) assert(0);
      }
      #define coef_coef(x) \
        coef(result, interp_order, num_coefs, j, x, num_coefs - 1 - i)
      coef_coef(0) = f0;
      if (interp_order > 0) coef_coef(1) = b;
      if (interp_order > 1) coef_coef(2) = c;
      if (interp_order > 2) coef_coef(3) = d;
      #undef coef_coef
      f2 = f1, f1 = f0;
    }
#ifdef SSE_
  if (q==1) {
    if (interp_order==2) {  //d150_2, u150_2
      sox_sample_t * const c = result;
      sox_sample_t t;
      for(i=0; i < num_coefs*num_phases/4; i++) {
        t = c[12*i + 1]; c[12*i + 1] = c[12*i + 3]; c[12*i + 3] = c[12*i + 9]; c[12*i + 9] = c[12*i + 5];  c[12*i + 5]  = c[12*i + 4]; c[12*i + 4] = t;
        t = c[12*i + 2]; c[12*i + 2] = c[12*i + 6]; c[12*i + 6] = c[12*i + 7]; c[12*i + 7] = c[12*i + 10]; c[12*i + 10] = c[12*i + 8]; c[12*i + 8] = t;
      }
    }
  }
  else {
    if (interp_order==1) {  //d120_1, u120_1
      sox_sample_t * const c = result;
      sox_sample_t t;
      for(i=0; i < num_coefs*num_phases/4; i++) {
        t = c[8*i + 1]; c[8*i + 1] = c[8*i + 2]; c[8*i + 2] = c[8*i + 4]; c[8*i + 4] = t;
        t = c[8*i + 3]; c[8*i + 3] = c[8*i + 6]; c[8*i + 6] = c[8*i + 5]; c[8*i + 5] = t;
      }
    }
  }
#endif
  return result;
}

typedef struct {    /* Data that are shared between channels and stages */
  sox_sample_t   * poly_fir_coefs;
  dft_filter_t half_band[2];   /* Div or mul by n */
} rate_shared_t;

struct stage;
typedef void (* stage_fn_t)(struct stage * input, fifo_t * output);
typedef struct stage {
  rate_shared_t * shared;
  fifo_t     fifo;
  int        pre;              /* Number of past samples to store */
  int        pre_post;         /* pre + number of future samples to store */
  int        preload;          /* Number of zero samples to pre-load the fifo */
  int        which;            /* Which, if any, of the 2 dft filters to use */
  stage_fn_t fn;
                               /* For poly_fir & spline: */
  union {                      /* 32bit.32bit fixed point arithmetic */
    #if defined(WORDS_BIGENDIAN)
    struct {int32_t integer; uint32_t fraction;} parts;
    #else
    struct {uint32_t fraction; int32_t integer;} parts;
    #endif
    int64_t all;
    #define MULT32 (65536.f * 65536.f)
  } at, step;
  int        divisor;          /* For step: > 1 for rational; 1 otherwise */
  double     out_in_ratio;
  int        rem, tuple;
} stage_t;

#define stage_occupancy(s) max(0, fifo_occupancy(&(s)->fifo) - (s)->pre_post)
#define stage_read_p(s) ((sox_sample_t *)fifo_read_ptr(&(s)->fifo) + (s)->pre)

#ifdef SSE_
#if defined(__GNUC__)
#define _MM_MALLOC_H_INCLUDED /*hack to avoid malloc clash*/
#ifndef _MM_ALIGN16
#define _MM_ALIGN16 __attribute__((aligned(16)))
#endif
#endif

#include <xmmintrin.h>
#ifdef SSE3_
#include <pmmintrin.h>
#endif

static _MM_ALIGN16 const unsigned long PCS_RNRN[4] = {0x00000000, 0x80000000, 0x00000000, 0x80000000};
static _MM_ALIGN16 const unsigned long PCS_NRNR[4] = {0x80000000, 0x00000000, 0x80000000, 0x00000000};
#define SIGN _mm_load_ps((float*)PCS_RNRN)

static __inline __m128 ZMUL2(__m128 a, __m128 b, __m128 sign)
{
#ifdef SSE3_
    // a = a1.r  a1.i  a2.r  a2.i
    // b = b1.r  b1.i  b2.r  b2.i
    __m128 ar;

    ar = _mm_moveldup_ps(a);        // ar = a1.r  a1.r  a2.r  a2.r
    a = _mm_movehdup_ps(a);         // a  = a1.i  a1.i  a2.i  a2.i
    ar = _mm_mul_ps(ar, b);         // ar = a1.r*b1.r  a1.r*b1.i  a2.r*b2.r  a2.r*b2.i
    
    b  = _mm_shuffle_ps(b, b, _MM_SHUFFLE(2, 3, 0, 1)); // b  = b1.i  b1.r  b2.i  b2.r
    a = _mm_mul_ps(a, b);           // ai = a1.i*b1.i  a1.i*b1.r  a2.i*b2.i  a2.i*b2.r

    return _mm_addsub_ps(ar, a);    // a1.r*b1.r-a1.i*b1.i  a1.r*b1.i+a1.i*b1.r  a2.r*b2.r-a2.i*b2.i  a2.r*b2.i+a2.i*b2.r
#else
    // a = a1.r  a1.i  a2.r  a2.i
    // b = b1.r  b1.i  b2.r  b2.i
    __m128 ar;

    ar = _mm_shuffle_ps(a, a, _MM_SHUFFLE(2, 2, 0, 0));     // ar = a1.r  a1.r  a2.r  a2.r
    a  = _mm_shuffle_ps(a, a, _MM_SHUFFLE(3, 3, 1, 1));     // ai = a1.i  a1.i  a2.i  a2.i
    ar = _mm_mul_ps(ar, b);                                 // ar = +a1.r*b1.r  +a1.r*b1.i  +a2.r*b2.r  +a2.r*b2.i
    
    a  = _mm_xor_ps(a, sign);                             // ai = a1.i  -a1.i  a2.i  -a2.i
    a  = _mm_mul_ps(a, b);                                // ai = a1.i*b1.r  -a1.i*b1.i  a2.i*b2.r  -a2.i*b2.i
    a  = _mm_shuffle_ps(a, a, _MM_SHUFFLE(2, 3, 0, 1));  // ai = -a1.i*b1.i  +a1.i*b1.r  -a2.i*b2.i  +a2.i*b2.r

    return _mm_add_ps(ar, a);   // a1.r*b1.r-a1.i*b1.i  a1.r*b1.i+a1.i*b1.r  a2.r*b2.r-a2.i*b2.i  a2.r*b2.i+a2.i*b2.r
#endif
}
#endif

#define CAT2(a,b) a##b
#if defined (SSE3_)
#define RATE_POSTFIX(a) CAT2(a, _SSE3)
#define FFT_POSTFIX(a)  CAT2(a, _SSE)
#elif defined (SSE_)
#define RATE_POSTFIX(a) CAT2(a, _SSE)
#define FFT_POSTFIX(a)  CAT2(a, _SSE)
#else
#define RATE_POSTFIX(a) CAT2(a, _generic)
#define FFT_POSTFIX(a)  CAT2(a, _generic)
#endif

#define ff_rdft_x FFT_POSTFIX(ff_rdft)
#define x_ff_rdft_x FFT_POSTFIX(x_ff_rdft)


static void half_sample(stage_t * p, fifo_t * output_fifo)
{
  sox_sample_t * output;
  int i, j, num_in = max(0, fifo_occupancy(&p->fifo));
  rate_shared_t const * s = p->shared;
  dft_filter_t const * f = &s->half_band[p->which];
  int const overlap = f->num_taps - 1;
#ifdef SSE_
  const float * const coeff = f->coefs;
  sox_sample_t tmp;
  __m128 coef, outp, sign;
  sign = SIGN;
#endif

  while (num_in >= f->dft_length) {
    sox_sample_t const * input = fifo_read_ptr(&p->fifo);
    fifo_read(&p->fifo, f->dft_length - overlap, NULL);
    num_in -= f->dft_length - overlap;

    output = fifo_reserve_aligned(output_fifo, f->dft_length);
    memcpy(output, input, f->dft_length * sizeof(*output));

    ff_rdft_x(f->dft_length, 1, output, f->tmp_buf);
#ifdef SSE_
    output[0] *= coeff[0];
    output[1] *= coeff[1];
    tmp = output[2];
    output[2] = coeff[2] * tmp - coeff[3] * output[3];
    output[3] = coeff[3] * tmp + coeff[2] * output[3];

    for (i = 4; i < f->dft_length; i += 4)
    {
      outp = _mm_load_ps(output+i);
      coef = _mm_load_ps(coeff+i);
      _mm_store_ps(output+i, ZMUL2(outp, coef, sign));
    }
#else
    output[0] *= f->coefs[0];
    output[1] *= f->coefs[1];
    for (i = 2; i < f->dft_length; i += 2) {
      sox_sample_t tmp = output[i];
      output[i  ] = f->coefs[i  ] * tmp - f->coefs[i+1] * output[i+1];
      output[i+1] = f->coefs[i+1] * tmp + f->coefs[i  ] * output[i+1];
    }
#endif
    ff_rdft_x(f->dft_length, -1, output, f->tmp_buf);

    for (j = 0, i = p->rem; i < f->dft_length - overlap; ++j, i += p->tuple)
      output[j] = output[i];
    p->rem = i - (f->dft_length - overlap);
    fifo_trim_by(output_fifo, f->dft_length - j);
  }
}

static void double_sample(stage_t * p, fifo_t * output_fifo)
{
  sox_sample_t * output;
  int i, j, num_in = max(0, fifo_occupancy(&p->fifo));
  rate_shared_t const * s = p->shared;
  dft_filter_t const * f = &s->half_band[p->which];
  int const overlap = f->num_taps - 1;
#ifdef SSE_
  const float * const coeff = f->coefs;
  sox_sample_t tmp;
  __m128 coef, outp, sign;
  sign = SIGN;
#endif

  while (p->rem + p->tuple * num_in >= f->dft_length) {
    div_t divd = div(f->dft_length - overlap - p->rem + p->tuple - 1, p->tuple);
    sox_sample_t const * input = fifo_read_ptr(&p->fifo);
    fifo_read(&p->fifo, divd.quot, NULL);
    num_in -= divd.quot;

    output = fifo_reserve_aligned(output_fifo, f->dft_length);
    fifo_trim_by(output_fifo, overlap);
    memset(output, 0, f->dft_length * sizeof(*output));
    for (j = 0, i = p->rem; i < f->dft_length; ++j, i += p->tuple)
      output[i] = input[j];
    p->rem = p->tuple - 1 - divd.rem;

    ff_rdft_x(f->dft_length, 1, output, f->tmp_buf);
#ifdef SSE_
    output[0] *= coeff[0];
    output[1] *= coeff[1];
    tmp = output[2];
    output[2] = coeff[2] * tmp - coeff[3] * output[3];
    output[3] = coeff[3] * tmp + coeff[2] * output[3];

    for (i = 4; i < f->dft_length; i += 4)
    {
      outp = _mm_load_ps(output+i);
      coef = _mm_load_ps(coeff+i);
      _mm_store_ps(output+i, ZMUL2(outp, coef, sign));
    }
#else
    output[0] *= f->coefs[0];
    output[1] *= f->coefs[1];
    for (i = 2; i < f->dft_length; i += 2) {
      sox_sample_t tmp = output[i];
      output[i  ] = f->coefs[i  ] * tmp - f->coefs[i+1] * output[i+1];
      output[i+1] = f->coefs[i+1] * tmp + f->coefs[i  ] * output[i+1];
    }
#endif
    ff_rdft_x(f->dft_length, -1, output, f->tmp_buf);
  }
}

static void init_dft_filter(rate_shared_t * p, unsigned which, int num_taps,
    sox_sample_t const h[], double Fp, double Fc, double Fn, double att,
    int multiplier, double phase, sox_bool allow_aliasing)
{
  dft_filter_t * f = &p->half_band[which];
  int dft_length, i;

  if (f->num_taps)
    return;
  if (h) { // for half_fir_coefs_low
    dft_length = lsx_set_dft_length(num_taps);
    f->coefs = lsx_aligned_calloc(dft_length, sizeof(*f->coefs));
    for (i = 0; i < num_taps; ++i)
      f->coefs[(i + dft_length - num_taps + 1) & (dft_length - 1)]
          = h[abs(num_taps / 2 - i)] / dft_length * 2 * multiplier;
    f->post_peak = num_taps / 2;
  }
  else {
    double * h2 = lsx_design_lpf(Fp, Fc, Fn, allow_aliasing, att, &num_taps, 0, -1.);

    if (phase != 50)
      lsx_fir_to_phase(&h2, &num_taps, &f->post_peak, phase);
    else f->post_peak = num_taps / 2;

    dft_length = lsx_set_dft_length(num_taps);
    f->coefs = lsx_aligned_calloc(dft_length, sizeof(*f->coefs));
    for (i = 0; i < num_taps; ++i)
      f->coefs[(i + dft_length - num_taps + 1) & (dft_length - 1)]
          = h2[i] / dft_length * 2 * multiplier;
    lsx_free(h2);
  }
  assert(num_taps & 1);
  f->num_taps = num_taps;
  f->dft_length = dft_length;
  f->tmp_buf = lsx_aligned_malloc(dft_length*sizeof(FFTComplex)/2);
  ff_rdft_x(dft_length, 1, f->coefs, f->tmp_buf);
}

#ifdef SSE_
#include "rate_filters_sse.h"
#else
#include "rate_filters_flt.h"
#endif

typedef struct {
  double     factor;
  uint64_t   samples_in, samples_out;
  int        level, input_stage_num, output_stage_num;
  sox_bool   upsample;
  stage_t    * stages;

  unsigned int in_samplerate, out_samplerate;
} rate_t;

#define pre_stage p->stages[-1]
#define last_stage p->stages[p->level]
#define post_stage p->stages[p->level + 1]

typedef enum {Default = -1, High = 3, Very} quality_t;

static void rate_init(rate_t * p, rate_shared_t * shared, double factor,
    quality_t quality, int interp_order, double phase, double bandwidth,
    sox_bool allow_aliasing, sox_bool old_behaviour)
{
  int i, tuple, mult, divisor = 1, two_or_three = 2;
  sox_bool two_factors = sox_false;
  const int max_divisor = 2048;      /* Keep coef table size ~< 500kb */
  double epsilon;

  assert(factor > 0);
  p->factor = factor;
  if (quality < High || quality > Very)
    quality = High;
  p->upsample = (sox_bool)(factor < 1);
  for (i = (int)(factor), p->level = 0; i >>= 1; ++p->level); /* log base 2 */
  factor /= 1 << (p->level + !p->upsample);
  epsilon = fabs((uint32_t)(factor * MULT32 + .5) / (factor * MULT32) - 1);
  for (i = 2; i <= max_divisor && divisor == 1; ++i) {
    double try_d = factor * i;
    int try_i = (int)(try_d + .5);
    if (fabs(try_i / try_d - 1) <= epsilon) { /* N.B. beware of long doubles */
      if (try_i == i) {
        factor = 1, divisor = 2;
        if (p->upsample)
          p->upsample = sox_false;
        else ++p->level;
      }
    else factor = try_i, divisor = i;
    }
  }
  if (!old_behaviour && factor == 3 && divisor == 4 && p->level == 1)
    two_or_three = 3;
  p->stages = (stage_t *)lsx_calloc((size_t)p->level + 4, sizeof(*p->stages)) + 1;
  for (i = -1; i <= p->level + 1; ++i) p->stages[i].shared = shared;
  last_stage.step.all = (int64_t)(factor * MULT32 + .5);
  last_stage.out_in_ratio = MULT32 * divisor / last_stage.step.all;

  if (divisor != 1)
    assert(!last_stage.step.parts.fraction);
  else
    assert(!last_stage.step.parts.integer);

  tuple = 1 + p->upsample;
  if (!old_behaviour && p->upsample) {
    if (factor == 2 && divisor == 3)
      two_factors = sox_true, tuple = divisor;
    else if (factor == 1) {
      if (divisor < 6)
        two_factors = sox_true, tuple = divisor;
      else for (i = 2; !two_factors && i < 20; ++i)
        if (!(divisor % i))
          two_factors = sox_true, tuple = i;
    }
  }
  mult = tuple;

  p->input_stage_num = -p->upsample;
  p->output_stage_num = p->level;

  if (two_or_three != 3 && last_stage.out_in_ratio != 2 && !two_factors) {
    poly_fir_t const * f;
    poly_fir1_t const * f1;
    int n = 2 * p->upsample + (quality - High);
    //if (interp_order < 0)
    interp_order = quality > High;
    interp_order = divisor == 1? 1 + interp_order : 0;
    last_stage.divisor = divisor;
    p->output_stage_num += 2;
    f = &poly_firs[n];
    f1 = &f->interp[interp_order];
    if (!last_stage.shared->poly_fir_coefs) {
      int num_taps = f->num_coefs, phases = divisor == 1? (1 << f1->phase_bits) : divisor;
      raw_coef_t * coefs = lsx_design_lpf(
          f->pass, f->stop, 1., sox_false, f->att, &num_taps, phases, -1.);
      assert(num_taps == f->num_coefs * phases - 1);
      last_stage.shared->poly_fir_coefs =
          prepare_coefs(coefs, f->num_coefs, phases, interp_order, mult, quality - High);
      lsx_free(coefs);
    }
    last_stage.fn = f1->fn;
    last_stage.pre_post = f->num_coefs - 1;
    last_stage.pre = 0;
    last_stage.preload = last_stage.pre_post >> 1;
    mult = 1;
  }
  {
    typedef struct {int len; sox_sample_t const * h; double bw, a;} filter_t;
    static filter_t const filters[] = {{0, NULL, .931, 125}, {0, NULL, .931, 150}};
    filter_t const * f = &filters[quality - High];
    double att = allow_aliasing? (34./33)* f->a : f->a; /* negate att degrade */
    double bw = bandwidth? 1 - (1 - bandwidth / 100) / LSX_TO_3dB : f->bw;
    double min = 1 - (allow_aliasing? LSX_MAX_TBW0A : LSX_MAX_TBW0) / 100;
    assert((size_t)(quality - High) < array_length(filters));
    init_dft_filter(shared, p->upsample, 0/*f->len*/, NULL/*f->h*/, bw, 1., (double)max(tuple, two_or_three), att, mult, phase, allow_aliasing);
    if (p->upsample) {
      pre_stage.fn = double_sample; /* Finish off setting up pre-stage */
      pre_stage.preload = shared->half_band[1].post_peak / tuple;
      pre_stage.rem     = shared->half_band[1].post_peak % tuple;
      pre_stage.tuple = tuple;
      pre_stage.which = 1;
      if (two_factors && divisor != tuple) {
        int other = divisor / tuple;
        ++p->output_stage_num;
        init_dft_filter(shared, 0, 0, NULL, 1., (double)tuple, (double)divisor, att, other, phase, allow_aliasing);
        last_stage.fn = double_sample;
        last_stage.preload = shared->half_band[0].post_peak / other;
        last_stage.rem     = shared->half_band[0].post_peak % other;
        last_stage.tuple = other;
        last_stage.which = 0;
      }
      else {
        /* Start setting up post-stage; TODO don't use dft for short filters */
        if ((1 - p->factor) / (1 - bw) > 2)
          init_dft_filter(shared, 0, 0, NULL, max(p->factor, min), 1., 2., att, 1, phase, allow_aliasing);
        else shared->half_band[0] = shared->half_band[1];
        if (two_factors && factor == 2) {
          ++p->output_stage_num;
          last_stage.fn = half_sample;
          last_stage.preload = shared->half_band[0].post_peak;
          last_stage.tuple = 2;
        } else {
          post_stage.fn = half_sample;
          post_stage.preload = shared->half_band[0].post_peak;
          post_stage.tuple = 2;
        }
      }
    }
    else {
      if (p->level > 0 && p->output_stage_num > p->level) {
        double pass = bw * divisor / factor / 2;
        if ((1 - pass) / (1 - bw) > 2)
          init_dft_filter(shared, 1, 0, NULL, max(pass, min), 1., 2., att, 1, phase, allow_aliasing);
      }
      post_stage.fn = half_sample;
      post_stage.preload = shared->half_band[0].post_peak;
      post_stage.tuple = two_or_three;
    }
  }
  if (p->level > 0) {
    stage_t * s = & p->stages[p->level - 1];
    if (shared->half_band[1].num_taps) {
      s->fn = half_sample;
      s->preload = shared->half_band[1].post_peak;
      s->tuple = 2;
      s->which = 1;
    }
    else *s = post_stage;
  }
  for (i = p->input_stage_num; i <= p->output_stage_num; ++i) {
    stage_t * s = &p->stages[i];
    if (i >= 0 && i < p->level - 1) {
      s->fn = half_sample_25;
      s->pre_post = 2 * (array_length(half_fir_coefs_25) - 1);
      s->preload = s->pre = s->pre_post >> 1;
    }
    fifo_create(&s->fifo, (int)sizeof(sox_sample_t));
    fifo_write(&s->fifo, s->preload, NULL);
  }
}

static void rate_process(rate_t * p)
{
  stage_t * stage = p->stages + p->input_stage_num;
  int i;

  for (i = p->input_stage_num; i < p->output_stage_num; ++i, ++stage)
    stage->fn(stage, &(stage+1)->fifo);
}

static sox_sample_t * rate_input(rate_t * p, sox_sample_t const * samples, size_t n)
{
  p->samples_in += n;
  while ((p->samples_in > p->in_samplerate) && (p->samples_out > p->out_samplerate))
  {
      p->samples_in -= p->in_samplerate;
      p->samples_out -= p->out_samplerate;
  }
  return fifo_write(&p->stages[p->input_stage_num].fifo, (int)n, samples);
}

static sox_sample_t const * rate_output(rate_t * p, sox_sample_t * samples, size_t * n)
{
  fifo_t * fifo = &p->stages[p->output_stage_num].fifo;
  p->samples_out += *n = min(*n, (size_t)fifo_occupancy(fifo));
  while ((p->samples_in > p->in_samplerate) && (p->samples_out > p->out_samplerate))
  {
      p->samples_in -= p->in_samplerate;
      p->samples_out -= p->out_samplerate;
  }
  return fifo_read(fifo, (int)*n, samples);
}

static void rate_flush(rate_t * p)
{
  uint64_t samples_out = (size_t)(p->samples_in / p->factor + .5);
  size_t remaining = samples_out - p->samples_out;

  if ((int)remaining > 0) {
    fifo_t * fifo = &p->stages[p->output_stage_num].fifo;
    while ((size_t)fifo_occupancy(fifo) < remaining) {
      rate_input(p, NULL, (size_t) 1024);
      rate_process(p);
    }
    fifo_trim_to(fifo, (int)remaining);
    p->samples_in = 0;
  }
}

static void rate_reset(rate_t * p)
{
  int i;
  for (i = p->input_stage_num; i <= p->output_stage_num; ++i) {
    fifo_clear(&p->stages[i].fifo);
    fifo_write(&p->stages[i].fifo, p->stages[i].preload, NULL);
    p->stages[i].at.all = 0;
  }
  p->samples_in = p->samples_out = 0;
}

static void rate_close_shared(rate_shared_t * shared)
{
  lsx_aligned_free(shared->half_band[0].coefs);
  if (shared->half_band[1].coefs != shared->half_band[0].coefs)
    lsx_aligned_free(shared->half_band[1].coefs);
  
  lsx_aligned_free(shared->half_band[0].tmp_buf);
  if (shared->half_band[1].tmp_buf != shared->half_band[0].tmp_buf)
    lsx_aligned_free(shared->half_band[1].tmp_buf);

  lsx_aligned_free(shared->poly_fir_coefs);
  memset(shared, 0, sizeof(*shared));
}

static void rate_close(rate_t * p)
{
  int i;
  for (i = p->input_stage_num; i <= p->output_stage_num; ++i)
    fifo_delete(&p->stages[i].fifo);

  lsx_free(p->stages - 1);
}

/*------------------------------- Wrapper for foobar2000 --------------------------------*/
/* Copyright (c) 2008-09 lvqcl.
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

#include "rate_i.h"

typedef struct RR_internal_tag
{
	RR_vtable2 x;

	int nchannels;
	rate_t* rate;
	rate_shared_t shared;
} RR_internal;

#define RR_init_x RATE_POSTFIX(RR_init)
#define RR_flow_x RATE_POSTFIX(RR_flow)
#define RR_drain_x RATE_POSTFIX(RR_drain)
#define RR_close_x RATE_POSTFIX(RR_close)
#define RR_reset_x RATE_POSTFIX(RR_reset)

static void RR_close_x(RR_handle *h) /*stop*/
{
	RR_internal* p = (RR_internal*)h;
	int i;
	if (p == NULL) return;

	for(i=0; i < p->nchannels; i++)
		rate_close(&p->rate[i]);
	lsx_free(p->rate);
	rate_close_shared(&p->shared);

	lsx_free(p);
}

static int RR_init_x(RR_handle* h, const RR_config* config, int nchannels) /*create+start*/
{
	RR_internal* p = (RR_internal*)h;
	int i;
	quality_t qual = (config->quality == RR_good) ? Very : High;
	if (p == NULL) return RR_NULLHANDLE;

	p->rate = (rate_t*) lsx_calloc(nchannels, sizeof(rate_t));
	p->nchannels = nchannels;

	for(i=0; i < nchannels; i++)
	{
		p->rate[i].in_samplerate = config->in_rate;
		p->rate[i].out_samplerate = config->out_rate;
		rate_init(&p->rate[i], &p->shared, (double)config->in_rate/(double)config->out_rate,
			qual, (int)0 - 1, config->phase, config->bandwidth, (sox_bool)config->allow_aliasing, sox_false);
	}

	return RR_OK;
}

static int RR_flow_x(RR_handle *h, const fb_sample_t* const ibuf, fb_sample_t* const obuf, size_t isamp, size_t osamp,
															size_t* i_used, size_t* o_done) /*flow*/
{
	RR_internal* p = (RR_internal*)h;
	size_t n;
	size_t itmp;
	int i;
	sox_sample_t const * s;
	size_t odone, odone2;
	assert(p);
	if (p == NULL) return RR_NULLHANDLE;
	if (ibuf == NULL) { isamp = 0; if (i_used) *i_used = 0; else i_used = &itmp; }
	if (isamp > 1048576) isamp = 1048576;

	for(i=0; i < p->nchannels; i++)
	{
		*i_used = 0;
		odone = osamp;
		s = rate_output(&p->rate[i], NULL, &odone);
		for (n = 0; n < odone; ++n) obuf[n*p->nchannels+i] = (fb_sample_t) s[n];
		*o_done = odone;

		if (isamp && odone < osamp)
		{
			sox_sample_t * t = rate_input(&p->rate[i], NULL, isamp);
			for (n = 0; n < isamp; ++n) t[n] = (sox_sample_t) ibuf[n*p->nchannels+i]; /*for (n = isamp; n; --n) t[isamp-n] = (sox_sample_t) ibuf[(isamp-n)*p->nchannels+i];*/

			rate_process(&p->rate[i]);
			*i_used = isamp;

			odone2 = osamp - odone;
			s = rate_output(&p->rate[i], NULL, &odone2);
			for (n = 0; n < odone2; ++n) obuf[(n+odone)*p->nchannels+i] = (fb_sample_t) s[n];
			*o_done += odone2;
		}
	}

	return RR_OK;
}

//in rate_flush:
//samples_out = p->samples_in / p->factor + .5;
//size_t remaining = samples_out - p->samples_out;

static int RR_drain_x(RR_handle *h) /*drain*/
{
	RR_internal* p = (RR_internal*)h;
	int i;
	assert(p);
	if (p == NULL) return RR_NULLHANDLE;

	for(i=0; i < p->nchannels; i++)
		rate_flush(&p->rate[i]);
	return RR_OK;
}

static int RR_reset_x(RR_handle *h) /*set all buffers to empty state*/
{
	RR_internal* p = (RR_internal*)h;
	int i;
	assert(p);
	if (p == NULL) return RR_NULLHANDLE;

	for(i=0; i < p->nchannels; i++)
		rate_reset(&p->rate[i]);

	return RR_OK;
}

#ifdef SSE_
RR_handle* RR_it_SSE(const RR_config* config, int nchannels)
#else
RR_handle* RR_it_float(const RR_config* config, int nchannels)
#endif
{
	RR_internal* p;
	p = (RR_internal*)lsx_calloc(1, sizeof (RR_internal) );
	if (p == NULL) return NULL;

	p->x.public_.flow = RR_flow_x;
	p->x.public_.drain = RR_drain_x;
	p->x.public_.reset = RR_reset_x;

	p->x.init = RR_init_x;
	p->x.close = RR_close_x;

	RR_init_x((RR_handle*)p, config, nchannels);

	return (RR_handle*)p;
}
