/* Effect: change sample rate     Copyright (c) 2008 robs@users.sourceforge.net
 * Copyright (C) 2008-11 lvqcl - SSE2/SSE3 intrinsics
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
   with several functions rewritten with SSE2/SSE3 intrinsincs
*/

#if defined(__GNUC__)
#ifndef _MM_ALIGN16
#define _MM_ALIGN16 __attribute__((aligned(16)))
#endif
typedef unsigned int intptr_t;
#endif


static _MM_ALIGN16 const double half_fir_coefs_25B[] = {
7.3979325233687461e-008, 1.1833367010222812e-006,
-8.1340436298087893e-007, -1.4284332593063177e-005, 4.6764104835321042e-006,
9.0580351350892191e-005, -1.8501044952475473e-005, -3.9872042837864422e-004,
5.6110366313398705e-005, 1.3634218103234187e-003, -1.3803431143314762e-004,
-3.8562347294894628e-003, 2.8491539998284476e-004, 9.4223774565849357e-003,
-5.0429677622613805e-004, -2.0673365323361139e-002, 7.7661461450703555e-004,
4.2764945027796687e-002, -1.0507348255277846e-003, -9.2035726038137103e-002,
1.2567743716165585e-003, 3.1333582318860204e-001, 4.9866643051942178e-001,
};
#define half_fir_coefs_25A (&half_fir_coefs_25B[22])

#include "rate_helpers_sse3.h"

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
        const double * input = ((double *)fifo_read_ptr(&p->fifo) + p->pre); //p->pre == 22
        int i;
        double * output = (double*)fifo_reserve(output_fifo, num_out);

        if ((((intptr_t)input)&0xf)==0)
        {                                                //in my tests we're always here...
            for (i = 0; i < num_out; ++i, input += 2)
            {
                __m128d XMM0, XMM1, XMM2, XMM3, XMM4, XMM5, XMM6, XMM7;
                const double * const coeff = half_fir_coefs_25A;

                SUM02(22, _mm_load_pd);

                SUM1(20, _mm_load_pd);
                SUM2(16, _mm_load_pd);
                SUM1(12, _mm_load_pd);
                SUM2( 8, _mm_load_pd);
                SUM1( 4, _mm_load_pd);

                XMM1 = _mm_load_sd(coeff);     //+0, ---
                XMM1 = _mm_mul_sd(XMM1, XMM7); //+0, +1
                XMM0 = _mm_hadd_pd(XMM0, XMM0);
                XMM0 = _mm_add_sd(XMM0, XMM1);

                _mm_store_sd(output+i, XMM0);
            }
        }
        else
        {
            for (i = 0; i < num_out; ++i, input += 2)
            {
                __m128d XMM0, XMM1, XMM2, XMM3, XMM4, XMM5, XMM6, XMM7;
                const double * const coeff = half_fir_coefs_25A;

                SUM02(22, _mm_lddqu_pd);

                SUM1(20, _mm_lddqu_pd);
                SUM2(16, _mm_lddqu_pd);
                SUM1(12, _mm_lddqu_pd);
                SUM2( 8, _mm_lddqu_pd);
                SUM1( 4, _mm_lddqu_pd);

                XMM1 = _mm_load_sd(coeff);     //+0, --
                XMM1 = _mm_mul_sd(XMM1, XMM7); //+0, +1
                XMM0 = _mm_hadd_pd(XMM0, XMM0);
                XMM0 = _mm_add_sd(XMM0, XMM1);

                _mm_store_sd(output+i, XMM0);
            }
        }

        fifo_read(&p->fifo, 2 * num_out, NULL);
    }
}


#define d150_l 38
#define u150_l 24

#define d150_2_b 10
#define u150_2_b 9

//#define M_N_O (int)(1.0001 + num_in*p->out_in_ratio)

static void d150_0(stage_t * p, fifo_t * output_fifo)
{
    int num_in = fifo_occupancy(&p->fifo) - p->pre_post;
    if (num_in > 0)
    {
        int i;
        //const int max_num_out = M_N_O;
        const int max_num_out = (num_in * p->divisor - p->at.parts.integer + p->step.parts.integer - 1)/p->step.parts.integer;
        const double * const input = ((double *)fifo_read_ptr(&p->fifo) + p->pre);
        double * output = (double*)fifo_reserve(output_fifo, max_num_out);
        div_t divided2;

        for (i = 0; p->at.parts.integer < num_in * p->divisor; ++i, p->at.parts.integer += p->step.parts.integer) {
            div_t divided = div1(p->at.parts.integer, p->divisor);
            const double * const at = input + divided.quot;
            const double * const coef = &p->shared->poly_fir_coefs[d150_l*divided.rem];
            __m128d XMM0, XMM1, XMM2, XMM3, XMM4, XMM5, XMM6;

            ADD2_to3(0);
            ADD2_3to4(4);
            ADD2_4to3(8);
            ADD2_3to4(12);
            ADD2_4to3(16);
            ADD2_3to4(20);
            ADD2_4to3(24);
            ADD2_3to4(28);
            ADD3_4to6(32);
            OUT_6;
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
        const double * const input = ((double *)fifo_read_ptr(&p->fifo) + p->pre);
        double * output = (double*)fifo_reserve(output_fifo, max_num_out);

        for (i = 0; p->at.parts.integer < num_in; ++i, p->at.all += p->step.all)
        {
            const double * const at = input + p->at.parts.integer;
            uint32_t fraction = p->at.parts.fraction;
            int phase = fraction >> (32 - d150_2_b);
            const double * const coef = &p->shared->poly_fir_coefs[d150_l*3*phase];
            __m128d XMM0, XMM1, XMM2, XMM3, XMM4, XMM5, XMM6;

            XMM0 = _mm_setzero_pd();
            XMM1 = _mm_set1_pd((fraction << d150_2_b) * (1 / MULT32));
            XMM2 = XMM1;
            XMM2 = _mm_mul_pd(XMM2, XMM2);

            IADD3( 0);
            IADD3( 2);
            IADD3( 4);
            IADD3( 6);
            IADD3( 8);
            IADD3(10);
            IADD3(12);
            IADD3(14);
            IADD3(16);
            IADD3(18);
            IADD3(20);
            IADD3(22);
            IADD3(24);
            IADD3(26);
            IADD3(28);
            IADD3(30);
            IADD3(32);
            IADD3(34);
            IADD3(36);

            XMM0 = _mm_hadd_pd(XMM0, XMM0);
            _mm_store_sd(output+i, XMM0);
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
        const double * const input = (double*)fifo_read_ptr(&p->fifo) + p->pre;
        double * output = (double*)fifo_reserve(output_fifo, max_num_out);
        div_t divided2;

        for (i = 0; p->at.parts.integer < num_in * p->divisor; ++i, p->at.parts.integer += p->step.parts.integer) {
            div_t divided = div1(p->at.parts.integer, p->divisor);
            const double * const at = input + divided.quot;
            const double * const coef = &p->shared->poly_fir_coefs[u150_l*divided.rem];
            __m128d XMM0, XMM1, XMM2, XMM3, XMM4;

            ADD2_to3(0);
            ADD2_3to4(4);
            ADD2_4to3(8);
            ADD2_3to4(12);
            ADD2_4to3(16);
            ADD2_3to4(20);
            OUT_4;
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
        const double * const input = ((double *)fifo_read_ptr(&p->fifo) + p->pre);
        double * output = (double*)fifo_reserve(output_fifo, max_num_out);

        for (i = 0; p->at.parts.integer < num_in; ++i, p->at.all += p->step.all)
        {
            const double * const at = input + p->at.parts.integer;
            uint32_t fraction = p->at.parts.fraction;
            int phase = fraction >> (32 - u150_2_b);
            const double * const coef = &p->shared->poly_fir_coefs[u150_l*3*phase];
            __m128d XMM0, XMM1, XMM2, XMM3, XMM4, XMM5, XMM6;

            XMM0 = _mm_setzero_pd();
            XMM1 = _mm_set1_pd((fraction << u150_2_b) * (1 / MULT32));
            XMM2 = XMM1;
            XMM2 = _mm_mul_pd(XMM2, XMM2);

            /*--------------------------------------------------------
            sum += (coef[0]*x*x + coef[1]*x + coef[2])*at[0];
            sum += (coef[3]*x*x + coef[4]*x + coef[5])*at[1];
            vvv
            sum += (coef[0]*x*x + coef[2]*x + coef[4])*at[0];
            sum += (coef[1]*x*x + coef[3]*x + coef[5])*at[1];
            vvv
            new1 = old3, new2 = old1, new3 = old4, new4 = old2;
            t = c[6*i + 1], c[6*i + 1] = c[6*i + 3], c[6*i + 3] = c[6*i + 4], c[6*i + 4] = c[6*i + 2], c[6*i + 2] = t;
            --------------------------------------------------------*/

            IADD3( 0);
            IADD3( 2);
            IADD3( 4);
            IADD3( 6);
            IADD3( 8);
            IADD3(10);
            IADD3(12);
            IADD3(14);
            IADD3(16);
            IADD3(18);
            IADD3(20);
            IADD3(22);

            XMM0 = _mm_hadd_pd(XMM0, XMM0);
            _mm_store_sd(output+i, XMM0);
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
  {d150_l,  1, 1.5, 165, {{0, d150_0}, {0, NULL}, {d150_2_b, d150_2}, {0, NULL}}},
  {u150_l, .5, 1.5, 174, {{0, u150_0}, {0, NULL}, {u150_2_b, u150_2}, {0, NULL}}},
};
