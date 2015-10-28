/* Effect: change sample rate     Copyright (c) 2008 robs@users.sourceforge.net
 * Copyright (C) 2011 lvqcl - SSE intrinsics
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
/* Note from lvqcl:
   This file is preprocessed version of original rate_filters.h
   with several functions rewritten with SSE intrinsincs
*/

#include <xmmintrin.h>
#if defined(SSE3_)
#include <pmmintrin.h>
#endif

#if defined(__GNUC__)
#ifndef _MM_ALIGN16
#define _MM_ALIGN16 __attribute__((aligned(16)))
#endif
typedef unsigned int intptr_t;
#endif

#pragma warning(disable : 4305)
static _MM_ALIGN16 const float half_fir_coefs_25B[] = {
0, 0, 7.3979325233687461e-008, 1.1833367010222812e-006,
-8.1340436298087893e-007, -1.4284332593063177e-005, 4.6764104835321042e-006,
9.0580351350892191e-005, -1.8501044952475473e-005, -3.9872042837864422e-004,
5.6110366313398705e-005, 1.3634218103234187e-003, -1.3803431143314762e-004,
-3.8562347294894628e-003, 2.8491539998284476e-004, 9.4223774565849357e-003,
-5.0429677622613805e-004, -2.0673365323361139e-002, 7.7661461450703555e-004,
4.2764945027796687e-002, -1.0507348255277846e-003, -9.2035726038137103e-002,
1.2567743716165585e-003, 3.1333582318860204e-001, 4.9866643051942178e-001,
};
#define half_fir_coefs_25A (half_fir_coefs_25B+24)
#define coeff half_fir_coefs_25A
#pragma warning(default : 4305)

#include "rate_helpers_sse.h"

static __inline div_t __cdecl div1(int numer, int denom)
{
    div_t result;

    result.quot = numer / denom;
    result.rem = numer % denom;

    return result;
}


static void half_sample_25(stage_t * p, fifo_t * output_fifo)
{
    int num_out = (fifo_occupancy(&p->fifo) - p->pre_post + 1) / 2;
    if (num_out > 0)
    {
        float const * input = ((float *)fifo_read_ptr(&p->fifo) + p->pre);
        int i;
        float * output = (float*)fifo_reserve(output_fifo, num_out);

        for (i = 0; i < num_out; ++i, input += 2)
        {
            __m128 xmm0, xmm1, xmm2, xmm3, xmm4, xmm5, xmm7;

            //sum = input[0] * coeff[0];
            xmm3 = _mm_lddqu_ps(input); //  0   1   2   3
            xmm7 = _mm_load_ss(coeff); // coeff[0]  '0'  '0'  '0'
            xmm7 = _mm_mul_ss(xmm7, xmm3); // == in[0]*coef[0]

            SUM1( 4);
            SUM2( 8);

            SUM1(12);
            SUM2(16);

            SUM1(20);
            SUM2(24);

            xmm7 = _mm_add_horz_ss(xmm7);
            _mm_store_ss(output+i, xmm7);
        }
        fifo_read(&p->fifo, 2 * num_out, NULL);
    }
}


#define d120_l 32
#define d150_l 40
#define u120_l 16
#define u150_l 20

#define d120_1_b 10
#define d150_2_b 10
#define u120_1_b 10
#define u150_2_b 9

//#define M_N_O (int)(1.0001 + num_in*p->out_in_ratio)

static void d120_0(stage_t * p, fifo_t * output_fifo)
{
    int num_in = fifo_occupancy(&p->fifo) - p->pre_post;
    if (num_in > 0)
    {
        int i;
        //const int max_num_out = M_N_O;
        const int max_num_out = (num_in * p->divisor - p->at.parts.integer + p->step.parts.integer - 1)/p->step.parts.integer;
        float const * input = ((float *)fifo_read_ptr(&p->fifo) + p->pre);
        float * output = (float*)fifo_reserve(output_fifo, max_num_out);
        div_t divided2;

        for (i = 0; p->at.parts.integer < num_in * p->divisor; ++i, p->at.parts.integer += p->step.parts.integer)
        {
            div_t divided = div1(p->at.parts.integer, p->divisor);
            float const * at = input + divided.quot;
            float const * coef = &p->shared->poly_fir_coefs[d120_l*divided.rem];
            __m128 xmm0, xmm1, xmm2, xmm3, xmm4;

            ADD0A();
            ADD0B(8);
            ADD0B(16);
            ADD0B(24);

            xmm0 = _mm_add_horz_ss(xmm0);
            _mm_store_ss(output+i, xmm0);
        }
        assert(max_num_out - i == 0);
        //fifo_trim_by(output_fifo, max_num_out - i);
        divided2 = div1(p->at.parts.integer, p->divisor);
        fifo_read(&p->fifo, divided2.quot, NULL);
        p->at.parts.integer -= divided2.quot * p->divisor;
    }
}


static void d120_1(stage_t * p, fifo_t * output_fifo)
{
    int num_in = fifo_occupancy(&p->fifo) - p->pre_post;
    if (num_in > 0)
    {
        int i;
        //const int max_num_out = M_N_O;
        const int max_num_out = ( (((int64_t)num_in) << 32) - p->at.all + p->step.all - 1)/p->step.all;
        float const * input = ((float *)fifo_read_ptr(&p->fifo) + p->pre);
        float * output = (float*)fifo_reserve(output_fifo, max_num_out);

        for (i = 0; p->at.parts.integer < num_in; ++i, p->at.all += p->step.all)
        {
            float const * at = input + p->at.parts.integer;
            uint32_t fraction = p->at.parts.fraction;
            int phase = fraction >> (32 - d120_1_b);
            float const * coef = &p->shared->poly_fir_coefs[d120_l*2*phase];

            __m128 sum, xmm1, xmm2, xmm3, xmm4, xmm5, xmm6;
            __m128 x = _mm_set1_ps((fraction << d120_1_b) * (1 / MULT32));
            sum = _mm_setzero_ps();

            ADD1p(0);
            ADD1p(8);
            ADD1p(16);
            ADD1p(24);

            sum = _mm_add_horz_ss(sum);
            _mm_store_ss(output+i, sum);
        }
        assert(max_num_out - i == 0);
        //fifo_trim_by(output_fifo, max_num_out - i);
        fifo_read(&p->fifo, p->at.parts.integer, NULL);
        p->at.parts.integer = 0;
    }
}


static void d150_0(stage_t * p, fifo_t * output_fifo)
{
    int num_in = fifo_occupancy(&p->fifo) - p->pre_post;
    if (num_in > 0)
    {
        int i;
        //const int max_num_out = M_N_O;
        const int max_num_out = (num_in * p->divisor - p->at.parts.integer + p->step.parts.integer - 1)/p->step.parts.integer;
        float const * input = ((float *)fifo_read_ptr(&p->fifo) + p->pre);
        float * output = (float*)fifo_reserve(output_fifo, max_num_out);
        div_t divided2;

        for (i = 0; p->at.parts.integer < num_in * p->divisor; ++i, p->at.parts.integer += p->step.parts.integer)
        {
            div_t divided = div1(p->at.parts.integer, p->divisor);
            float const * at = input + divided.quot;
            float const * coef = &p->shared->poly_fir_coefs[d150_l*divided.rem];
            __m128 xmm0, xmm1, xmm2, xmm3, xmm4;

            ADD0A();
            ADD0B(8);
            ADD0B(16);
            ADD0B(24);
            ADD0B(32);

            xmm0 = _mm_add_horz_ss(xmm0);
            _mm_store_ss(output+i, xmm0);
        }
        assert(max_num_out - i == 0);
        //fifo_trim_by(output_fifo, max_num_out - i);
        divided2 = div1(p->at.parts.integer, p->divisor);
        fifo_read(&p->fifo, divided2.quot, NULL);
        p->at.parts.integer -= divided2.quot * p->divisor;
    }
}


static void d150_2(stage_t * p, fifo_t * output_fifo)
{
    int num_in = fifo_occupancy(&p->fifo) - p->pre_post;
    if (num_in > 0)
    {
        int i;
        //const int max_num_out = M_N_O;
        const int max_num_out = ( (((int64_t)num_in) << 32) - p->at.all + p->step.all - 1)/p->step.all;
        float const * input = ((float *)fifo_read_ptr(&p->fifo) + p->pre);
        float * output = (float*)fifo_reserve(output_fifo, max_num_out);

        for (i = 0; p->at.parts.integer < num_in; ++i, p->at.all += p->step.all)
        {
            float const * at = input + p->at.parts.integer;
            uint32_t fraction = p->at.parts.fraction;
            int phase = fraction >> (32 - d150_2_b);
            float const * coef = &p->shared->poly_fir_coefs[d150_l*3*phase];

            __m128 sum, xmm1, xmm2, xmm3, xmm4;
            __m128 x = _mm_set1_ps((fraction << d150_2_b) * (1 / MULT32));
            __m128 x2 = _mm_mul_ps(x, x);
            sum = _mm_setzero_ps();

            ADD2( 0);
            ADD2( 4);
            ADD2( 8);
            ADD2(12);
            ADD2(16);
            ADD2(20);
            ADD2(24);
            ADD2(28);
            ADD2(32);
            ADD2(36);

            sum = _mm_add_horz_ss(sum);
            _mm_store_ss(output+i, sum);
        }
        assert(max_num_out - i == 0);
        //fifo_trim_by(output_fifo, max_num_out - i);
        fifo_read(&p->fifo, p->at.parts.integer, NULL);
        p->at.parts.integer = 0;
    }
}


static void u120_0(stage_t * p, fifo_t * output_fifo)
{
    int num_in = fifo_occupancy(&p->fifo) - p->pre_post;
    if (num_in > 0)
    {
        int i;
        //const int max_num_out = M_N_O;
        const int max_num_out = (num_in * p->divisor - p->at.parts.integer + p->step.parts.integer - 1)/p->step.parts.integer;
        float const * input = ((float *)fifo_read_ptr(&p->fifo) + p->pre);
        float * output = (float*)fifo_reserve(output_fifo, max_num_out);
        div_t divided2;

        for (i = 0; p->at.parts.integer < num_in * p->divisor; ++i, p->at.parts.integer += p->step.parts.integer)
        {
            div_t divided = div1(p->at.parts.integer, p->divisor);
            float const * at = input + divided.quot;
            float const * coef = &p->shared->poly_fir_coefs[u120_l*divided.rem];
            __m128 xmm0, xmm1, xmm2, xmm3, xmm4;

            ADD0A();
            ADD0B(8);

            xmm0 = _mm_add_horz_ss(xmm0);
            _mm_store_ss(output+i, xmm0);
        }
        assert(max_num_out - i == 0);
        //fifo_trim_by(output_fifo, max_num_out - i);
        divided2 = div1(p->at.parts.integer, p->divisor);
        fifo_read(&p->fifo, divided2.quot, NULL);
        p->at.parts.integer -= divided2.quot * p->divisor;
    }
}


static void u120_1(stage_t * p, fifo_t * output_fifo)
{
    int num_in = fifo_occupancy(&p->fifo) - p->pre_post;
    if (num_in > 0)
    {
        int i;
        //const int max_num_out = M_N_O;
        const int max_num_out = ( (((int64_t)num_in) << 32) - p->at.all + p->step.all - 1)/p->step.all;
        float const * input = ((float *)fifo_read_ptr(&p->fifo) + p->pre);
        float * output = (float*)fifo_reserve(output_fifo, max_num_out);

        for (i = 0; p->at.parts.integer < num_in; ++i, p->at.all += p->step.all)
        {
            float const * at = input + p->at.parts.integer;
            uint32_t fraction = p->at.parts.fraction;
            int phase = fraction >> (32 - u120_1_b);
            float const * coef = &p->shared->poly_fir_coefs[u120_l*2*phase];

            __m128 sum, xmm1, xmm2, xmm3, xmm4, xmm5, xmm6;
            __m128 x = _mm_set1_ps((fraction << u120_1_b) * (1 / MULT32));
            sum = _mm_setzero_ps();

            ADD1p(0);
            ADD1p(8);

            sum = _mm_add_horz_ss(sum);
            _mm_store_ss(output+i, sum);
        }
        assert(max_num_out - i == 0);
        //fifo_trim_by(output_fifo, max_num_out - i);
        fifo_read(&p->fifo, p->at.parts.integer, NULL);
        p->at.parts.integer = 0;
    }
}


static void u150_0(stage_t * p, fifo_t * output_fifo)
{
    int num_in = fifo_occupancy(&p->fifo) - p->pre_post;
    if (num_in > 0)
    {
        int i;
        //const int max_num_out = M_N_O;
        const int max_num_out = (num_in * p->divisor - p->at.parts.integer + p->step.parts.integer - 1)/p->step.parts.integer;
        float const * input = ((float *)fifo_read_ptr(&p->fifo) + p->pre);
        float * output = (float*)fifo_reserve(output_fifo, max_num_out);
        div_t divided2;

        for (i = 0; p->at.parts.integer < num_in * p->divisor; ++i, p->at.parts.integer += p->step.parts.integer)
        {
            div_t divided = div1(p->at.parts.integer, p->divisor);
            float const * at = input + divided.quot;
            float const * coef = &p->shared->poly_fir_coefs[u150_l*divided.rem];
            __m128 xmm0, xmm1, xmm2, xmm3, xmm4;

            ADD0A();
            ADD0B(8);

            xmm1 = _mm_load_ps(coef+16);
            xmm2 = _mm_lddqu_ps(at+16);
            xmm1 = _mm_mul_ps(xmm1, xmm2);
            xmm0 = _mm_add_ps(xmm0, xmm1);

            xmm0 = _mm_add_horz_ss(xmm0);
            _mm_store_ss(output+i, xmm0);
        }
        assert(max_num_out - i == 0);
        //fifo_trim_by(output_fifo, max_num_out - i);
        divided2 = div1(p->at.parts.integer, p->divisor);
        fifo_read(&p->fifo, divided2.quot, NULL);
        p->at.parts.integer -= divided2.quot * p->divisor;
    }
}

static void u150_2(stage_t * p, fifo_t * output_fifo)
{
    int num_in = fifo_occupancy(&p->fifo) - p->pre_post;
    if (num_in > 0)
    {
        int i;
        //const int max_num_out = M_N_O;
        const int max_num_out = ( (((int64_t)num_in) << 32) - p->at.all + p->step.all - 1)/p->step.all;
        float const * input = ((float *)fifo_read_ptr(&p->fifo) + p->pre);
        float * output = (float*)fifo_reserve(output_fifo, max_num_out);

        for (i = 0; p->at.parts.integer < num_in; ++i, p->at.all += p->step.all)
        {
            float const * at = input + p->at.parts.integer;
            uint32_t fraction = p->at.parts.fraction;
            int phase = fraction >> (32 - u150_2_b);
            float const * coef = &p->shared->poly_fir_coefs[u150_l*3*phase];

            __m128 sum, xmm1, xmm2, xmm3, xmm4;
            __m128 x = _mm_set1_ps((fraction << u150_2_b) * (1 / MULT32));
            __m128 x2 = _mm_mul_ps(x, x);
            sum = _mm_setzero_ps();

            //sum += (coef[0]*x*x + coef[ 1]*x + coef[ 2])*at[0];
            //sum += (coef[3]*x*x + coef[ 4]*x + coef[ 5])*at[1];
            //sum += (coef[6]*x*x + coef[ 7]*x + coef[ 8])*at[2];
            //sum += (coef[9]*x*x + coef[10]*x + coef[11])*at[3];
            //   vvvv
            //sum += (coef[0]*x*x + coef[4]*x + coef[ 8])*at[0];
            //sum += (coef[1]*x*x + coef[5]*x + coef[ 9])*at[1];
            //sum += (coef[2]*x*x + coef[6]*x + coef[10])*at[2];
            //sum += (coef[3]*x*x + coef[7]*x + coef[11])*at[3];
            /*
            t = c[1]; c[1] = c[3]; c[3] = c[9]; c[9] = c[5]; c[5] = c[4]; c[4] = t;
            t = c[2]; c[2] = c[6]; c[6] = c[7]; c[7] = c[10]; c[10] = c[8]; c[8] = t;
            */

            ADD2( 0);
            ADD2( 4);
            ADD2( 8);
            ADD2(12);
            ADD2(16);

            sum = _mm_add_horz_ss(sum);
            _mm_store_ss(output+i, sum);
        }
        assert(max_num_out - i == 0);
        //fifo_trim_by(output_fifo, max_num_out - i);
        fifo_read(&p->fifo, p->at.parts.integer, NULL);
        p->at.parts.integer = 0;
    }
}


typedef struct {int phase_bits; stage_fn_t fn;} poly_fir1_t;
typedef struct {int num_coefs; double pass, stop, att; poly_fir1_t interp[4];} poly_fir_t;
static poly_fir_t const poly_firs[] = {
  {d120_l,  1, 1.5, 133, {{0, d120_0}, {d120_1_b, d120_1}, {       0,   NULL}, {0, NULL}}},
  {d150_l,  1, 1.5, 150, {{0, d150_0}, {       0,   NULL}, {d150_2_b, d150_2}, {0, NULL}}},
  {u120_l, .5, 1.5, 125, {{0, u120_0}, {u120_1_b, u120_1}, {       0,   NULL}, {0, NULL}}},
  {u150_l, .5, 1.5, 150, {{0, u150_0}, {       0,   NULL}, {u150_2_b, u150_2}, {0, NULL}}},
};
