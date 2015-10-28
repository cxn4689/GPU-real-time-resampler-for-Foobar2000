/* Effect: change sample rate     Copyright (c) 2008 robs@users.sourceforge.net
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
   This file is a preprocessed and modified version of original rate_filters.h
*/

#pragma warning(disable : 4305)
static const float half_fir_coefs_25[] = {
  4.9866643051942178e-001, 3.1333582318860204e-001, 1.2567743716165585e-003,
 -9.2035726038137103e-002, -1.0507348255277846e-003, 4.2764945027796687e-002,
  7.7661461450703555e-004, -2.0673365323361139e-002, -5.0429677622613805e-004,
  9.4223774565849357e-003, 2.8491539998284476e-004, -3.8562347294894628e-003,
 -1.3803431143314762e-004, 1.3634218103234187e-003, 5.6110366313398705e-005,
 -3.9872042837864422e-004, -1.8501044952475473e-005, 9.0580351350892191e-005,
  4.6764104835321042e-006, -1.4284332593063177e-005, -8.1340436298087893e-007,
  1.1833367010222812e-006, 7.3979325233687461e-008,
};
#pragma warning(default : 4305)


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
    int i;
    float const * input = ((float*)fifo_read_ptr(&p->fifo) + p->pre);
    float * output = (float*)fifo_reserve(output_fifo, num_out);

    for (i = 0; i < num_out; ++i, input += 2)
    {
      int j = 1;
      float sum = input[0] * half_fir_coefs_25[0];
      sum += (input[-j] + input[j]) * half_fir_coefs_25[j], ++j;
      sum += (input[-j] + input[j]) * half_fir_coefs_25[j], ++j;
      sum += (input[-j] + input[j]) * half_fir_coefs_25[j], ++j;
      sum += (input[-j] + input[j]) * half_fir_coefs_25[j], ++j;
      sum += (input[-j] + input[j]) * half_fir_coefs_25[j], ++j;
      sum += (input[-j] + input[j]) * half_fir_coefs_25[j], ++j;
      sum += (input[-j] + input[j]) * half_fir_coefs_25[j], ++j;
      sum += (input[-j] + input[j]) * half_fir_coefs_25[j], ++j;
      sum += (input[-j] + input[j]) * half_fir_coefs_25[j], ++j;
      sum += (input[-j] + input[j]) * half_fir_coefs_25[j], ++j;
      sum += (input[-j] + input[j]) * half_fir_coefs_25[j], ++j;
      sum += (input[-j] + input[j]) * half_fir_coefs_25[j], ++j;
      sum += (input[-j] + input[j]) * half_fir_coefs_25[j], ++j;
      sum += (input[-j] + input[j]) * half_fir_coefs_25[j], ++j;
      sum += (input[-j] + input[j]) * half_fir_coefs_25[j], ++j;
      sum += (input[-j] + input[j]) * half_fir_coefs_25[j], ++j;
      sum += (input[-j] + input[j]) * half_fir_coefs_25[j], ++j;
      sum += (input[-j] + input[j]) * half_fir_coefs_25[j], ++j;
      sum += (input[-j] + input[j]) * half_fir_coefs_25[j], ++j;
      sum += (input[-j] + input[j]) * half_fir_coefs_25[j], ++j;
      sum += (input[-j] + input[j]) * half_fir_coefs_25[j], ++j;
      sum += (input[-j] + input[j]) * half_fir_coefs_25[j], ++j;

      assert(j == array_length(half_fir_coefs_25));
      output[i] = sum;
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
    
      float sum = 0;
      int j = 0;
    
      sum += coef[j] *at[j], ++j;
      sum += coef[j] *at[j], ++j;
      sum += coef[j] *at[j], ++j;
      sum += coef[j] *at[j], ++j;
    
      sum += coef[j] *at[j], ++j;
      sum += coef[j] *at[j], ++j;
      sum += coef[j] *at[j], ++j;
      sum += coef[j] *at[j], ++j;
    
      sum += coef[j] *at[j], ++j;
      sum += coef[j] *at[j], ++j;
      sum += coef[j] *at[j], ++j;
      sum += coef[j] *at[j], ++j;
    
      sum += coef[j] *at[j], ++j;
      sum += coef[j] *at[j], ++j;
      sum += coef[j] *at[j], ++j;
      sum += coef[j] *at[j], ++j;
    
      sum += coef[j] *at[j], ++j;
      sum += coef[j] *at[j], ++j;
      sum += coef[j] *at[j], ++j;
      sum += coef[j] *at[j], ++j;
    
      sum += coef[j] *at[j], ++j;
      sum += coef[j] *at[j], ++j;
      sum += coef[j] *at[j], ++j;
      sum += coef[j] *at[j], ++j;
    
      sum += coef[j] *at[j], ++j;
      sum += coef[j] *at[j], ++j;
      sum += coef[j] *at[j], ++j;
      sum += coef[j] *at[j], ++j;
    
      sum += coef[j] *at[j], ++j;
      sum += coef[j] *at[j], ++j;
      sum += coef[j] *at[j], ++j;
      sum += coef[j] *at[j], ++j;

      assert(j == d120_l);
      output[i] = sum;
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

      float x = (float) (fraction << d120_1_b) * (1 / MULT32);

      float sum = 0;
      int j = 0;
    
      sum += (coef[2*j]*x + coef[2*j + 1])*at[j], ++j;
      sum += (coef[2*j]*x + coef[2*j + 1])*at[j], ++j;
      sum += (coef[2*j]*x + coef[2*j + 1])*at[j], ++j;
      sum += (coef[2*j]*x + coef[2*j + 1])*at[j], ++j;
    
      sum += (coef[2*j]*x + coef[2*j + 1])*at[j], ++j;
      sum += (coef[2*j]*x + coef[2*j + 1])*at[j], ++j;
      sum += (coef[2*j]*x + coef[2*j + 1])*at[j], ++j;
      sum += (coef[2*j]*x + coef[2*j + 1])*at[j], ++j;
    
      sum += (coef[2*j]*x + coef[2*j + 1])*at[j], ++j;
      sum += (coef[2*j]*x + coef[2*j + 1])*at[j], ++j;
      sum += (coef[2*j]*x + coef[2*j + 1])*at[j], ++j;
      sum += (coef[2*j]*x + coef[2*j + 1])*at[j], ++j;
    
      sum += (coef[2*j]*x + coef[2*j + 1])*at[j], ++j;
      sum += (coef[2*j]*x + coef[2*j + 1])*at[j], ++j;
      sum += (coef[2*j]*x + coef[2*j + 1])*at[j], ++j;
      sum += (coef[2*j]*x + coef[2*j + 1])*at[j], ++j;
    
      sum += (coef[2*j]*x + coef[2*j + 1])*at[j], ++j;
      sum += (coef[2*j]*x + coef[2*j + 1])*at[j], ++j;
      sum += (coef[2*j]*x + coef[2*j + 1])*at[j], ++j;
      sum += (coef[2*j]*x + coef[2*j + 1])*at[j], ++j;
    
      sum += (coef[2*j]*x + coef[2*j + 1])*at[j], ++j;
      sum += (coef[2*j]*x + coef[2*j + 1])*at[j], ++j;
      sum += (coef[2*j]*x + coef[2*j + 1])*at[j], ++j;
      sum += (coef[2*j]*x + coef[2*j + 1])*at[j], ++j;
    
      sum += (coef[2*j]*x + coef[2*j + 1])*at[j], ++j;
      sum += (coef[2*j]*x + coef[2*j + 1])*at[j], ++j;
      sum += (coef[2*j]*x + coef[2*j + 1])*at[j], ++j;
      sum += (coef[2*j]*x + coef[2*j + 1])*at[j], ++j;
    
      sum += (coef[2*j]*x + coef[2*j + 1])*at[j], ++j;
      sum += (coef[2*j]*x + coef[2*j + 1])*at[j], ++j;
      sum += (coef[2*j]*x + coef[2*j + 1])*at[j], ++j;
      sum += (coef[2*j]*x + coef[2*j + 1])*at[j], ++j;

      assert(j == d120_l);
      output[i] = sum;
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
    
      float sum = 0;
      int j = 0;
      sum += coef[j] *at[j], ++j;
      sum += coef[j] *at[j], ++j;
      sum += coef[j] *at[j], ++j;
      sum += coef[j] *at[j], ++j;
    
      sum += coef[j] *at[j], ++j;
      sum += coef[j] *at[j], ++j;
      sum += coef[j] *at[j], ++j;
      sum += coef[j] *at[j], ++j;
    
      sum += coef[j] *at[j], ++j;
      sum += coef[j] *at[j], ++j;
      sum += coef[j] *at[j], ++j;
      sum += coef[j] *at[j], ++j;
    
      sum += coef[j] *at[j], ++j;
      sum += coef[j] *at[j], ++j;
      sum += coef[j] *at[j], ++j;
      sum += coef[j] *at[j], ++j;
    
      sum += coef[j] *at[j], ++j;
      sum += coef[j] *at[j], ++j;
      sum += coef[j] *at[j], ++j;
      sum += coef[j] *at[j], ++j;
    
      sum += coef[j] *at[j], ++j;
      sum += coef[j] *at[j], ++j;
      sum += coef[j] *at[j], ++j;
      sum += coef[j] *at[j], ++j;
    
      sum += coef[j] *at[j], ++j;
      sum += coef[j] *at[j], ++j;
      sum += coef[j] *at[j], ++j;
      sum += coef[j] *at[j], ++j;
    
      sum += coef[j] *at[j], ++j;
      sum += coef[j] *at[j], ++j;
      sum += coef[j] *at[j], ++j;
      sum += coef[j] *at[j], ++j;
    
      sum += coef[j] *at[j], ++j;
      sum += coef[j] *at[j], ++j;
      sum += coef[j] *at[j], ++j;
      sum += coef[j] *at[j], ++j;
    
      sum += coef[j] *at[j], ++j;
      sum += coef[j] *at[j], ++j;
      sum += coef[j] *at[j], ++j;
      sum += coef[j] *at[j], ++j;

      assert(j == d150_l);
      output[i] = sum;
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

      float x = (float) (fraction << d150_2_b) * (1 / MULT32);

      float sum = 0;
      int j = 0;
      sum += ((coef[3*j] *x + coef[3*j + 1])*x + coef[3*j + 2])*at[j], ++j;
      sum += ((coef[3*j] *x + coef[3*j + 1])*x + coef[3*j + 2])*at[j], ++j;
      sum += ((coef[3*j] *x + coef[3*j + 1])*x + coef[3*j + 2])*at[j], ++j;
      sum += ((coef[3*j] *x + coef[3*j + 1])*x + coef[3*j + 2])*at[j], ++j;
    
      sum += ((coef[3*j] *x + coef[3*j + 1])*x + coef[3*j + 2])*at[j], ++j;
      sum += ((coef[3*j] *x + coef[3*j + 1])*x + coef[3*j + 2])*at[j], ++j;
      sum += ((coef[3*j] *x + coef[3*j + 1])*x + coef[3*j + 2])*at[j], ++j;
      sum += ((coef[3*j] *x + coef[3*j + 1])*x + coef[3*j + 2])*at[j], ++j;
    
      sum += ((coef[3*j] *x + coef[3*j + 1])*x + coef[3*j + 2])*at[j], ++j;
      sum += ((coef[3*j] *x + coef[3*j + 1])*x + coef[3*j + 2])*at[j], ++j;
      sum += ((coef[3*j] *x + coef[3*j + 1])*x + coef[3*j + 2])*at[j], ++j;
      sum += ((coef[3*j] *x + coef[3*j + 1])*x + coef[3*j + 2])*at[j], ++j;
    
      sum += ((coef[3*j] *x + coef[3*j + 1])*x + coef[3*j + 2])*at[j], ++j;
      sum += ((coef[3*j] *x + coef[3*j + 1])*x + coef[3*j + 2])*at[j], ++j;
      sum += ((coef[3*j] *x + coef[3*j + 1])*x + coef[3*j + 2])*at[j], ++j;
      sum += ((coef[3*j] *x + coef[3*j + 1])*x + coef[3*j + 2])*at[j], ++j;
    
      sum += ((coef[3*j] *x + coef[3*j + 1])*x + coef[3*j + 2])*at[j], ++j;
      sum += ((coef[3*j] *x + coef[3*j + 1])*x + coef[3*j + 2])*at[j], ++j;
      sum += ((coef[3*j] *x + coef[3*j + 1])*x + coef[3*j + 2])*at[j], ++j;
      sum += ((coef[3*j] *x + coef[3*j + 1])*x + coef[3*j + 2])*at[j], ++j;
    
      sum += ((coef[3*j] *x + coef[3*j + 1])*x + coef[3*j + 2])*at[j], ++j;
      sum += ((coef[3*j] *x + coef[3*j + 1])*x + coef[3*j + 2])*at[j], ++j;
      sum += ((coef[3*j] *x + coef[3*j + 1])*x + coef[3*j + 2])*at[j], ++j;
      sum += ((coef[3*j] *x + coef[3*j + 1])*x + coef[3*j + 2])*at[j], ++j;
    
      sum += ((coef[3*j] *x + coef[3*j + 1])*x + coef[3*j + 2])*at[j], ++j;
      sum += ((coef[3*j] *x + coef[3*j + 1])*x + coef[3*j + 2])*at[j], ++j;
      sum += ((coef[3*j] *x + coef[3*j + 1])*x + coef[3*j + 2])*at[j], ++j;
      sum += ((coef[3*j] *x + coef[3*j + 1])*x + coef[3*j + 2])*at[j], ++j;
    
      sum += ((coef[3*j] *x + coef[3*j + 1])*x + coef[3*j + 2])*at[j], ++j;
      sum += ((coef[3*j] *x + coef[3*j + 1])*x + coef[3*j + 2])*at[j], ++j;
      sum += ((coef[3*j] *x + coef[3*j + 1])*x + coef[3*j + 2])*at[j], ++j;
      sum += ((coef[3*j] *x + coef[3*j + 1])*x + coef[3*j + 2])*at[j], ++j;
    
      sum += ((coef[3*j] *x + coef[3*j + 1])*x + coef[3*j + 2])*at[j], ++j;
      sum += ((coef[3*j] *x + coef[3*j + 1])*x + coef[3*j + 2])*at[j], ++j;
      sum += ((coef[3*j] *x + coef[3*j + 1])*x + coef[3*j + 2])*at[j], ++j;
      sum += ((coef[3*j] *x + coef[3*j + 1])*x + coef[3*j + 2])*at[j], ++j;
    
      sum += ((coef[3*j] *x + coef[3*j + 1])*x + coef[3*j + 2])*at[j], ++j;
      sum += ((coef[3*j] *x + coef[3*j + 1])*x + coef[3*j + 2])*at[j], ++j;
      sum += ((coef[3*j] *x + coef[3*j + 1])*x + coef[3*j + 2])*at[j], ++j;
      sum += ((coef[3*j] *x + coef[3*j + 1])*x + coef[3*j + 2])*at[j], ++j;

      assert(j == d150_l);
      output[i] = sum;
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
      float sum = 0;
      int j = 0;
    
      sum += coef[j] *at[j], ++j;
      sum += coef[j] *at[j], ++j;
      sum += coef[j] *at[j], ++j;
      sum += coef[j] *at[j], ++j;
    
      sum += coef[j] *at[j], ++j;
      sum += coef[j] *at[j], ++j;
      sum += coef[j] *at[j], ++j;
      sum += coef[j] *at[j], ++j;
    
      sum += coef[j] *at[j], ++j;
      sum += coef[j] *at[j], ++j;
      sum += coef[j] *at[j], ++j;
      sum += coef[j] *at[j], ++j;
    
      sum += coef[j] *at[j], ++j;
      sum += coef[j] *at[j], ++j;
      sum += coef[j] *at[j], ++j;
      sum += coef[j] *at[j], ++j;

      assert(j == u120_l);
      output[i] = sum;
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

      float x = (float) (fraction << u120_1_b) * (1 / MULT32);

      float sum = 0;
      int j = 0;
    
      sum += (coef[2*j]*x + coef[2*j + 1])*at[j], ++j;
      sum += (coef[2*j]*x + coef[2*j + 1])*at[j], ++j;
      sum += (coef[2*j]*x + coef[2*j + 1])*at[j], ++j;
      sum += (coef[2*j]*x + coef[2*j + 1])*at[j], ++j;
    
      sum += (coef[2*j]*x + coef[2*j + 1])*at[j], ++j;
      sum += (coef[2*j]*x + coef[2*j + 1])*at[j], ++j;
      sum += (coef[2*j]*x + coef[2*j + 1])*at[j], ++j;
      sum += (coef[2*j]*x + coef[2*j + 1])*at[j], ++j;
    
      sum += (coef[2*j]*x + coef[2*j + 1])*at[j], ++j;
      sum += (coef[2*j]*x + coef[2*j + 1])*at[j], ++j;
      sum += (coef[2*j]*x + coef[2*j + 1])*at[j], ++j;
      sum += (coef[2*j]*x + coef[2*j + 1])*at[j], ++j;
    
      sum += (coef[2*j]*x + coef[2*j + 1])*at[j], ++j;
      sum += (coef[2*j]*x + coef[2*j + 1])*at[j], ++j;
      sum += (coef[2*j]*x + coef[2*j + 1])*at[j], ++j;
      sum += (coef[2*j]*x + coef[2*j + 1])*at[j], ++j;

      assert(j == u120_l);
      output[i] = sum;
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

    for (i = 0; p->at.parts.integer < num_in * p->divisor; ++i, p->at.parts.integer += p->step.parts.integer) {
      div_t divided = div1(p->at.parts.integer, p->divisor);
      float const * at = input + divided.quot;
      float const * coef = &p->shared->poly_fir_coefs[u150_l*divided.rem];
      float sum = 0;
      int j = 0;

      sum += coef[j] *at[j], ++j;
      sum += coef[j] *at[j], ++j;
      sum += coef[j] *at[j], ++j;
      sum += coef[j] *at[j], ++j;
    
      sum += coef[j] *at[j], ++j;
      sum += coef[j] *at[j], ++j;
      sum += coef[j] *at[j], ++j;
      sum += coef[j] *at[j], ++j;
    
      sum += coef[j] *at[j], ++j;
      sum += coef[j] *at[j], ++j;
      sum += coef[j] *at[j], ++j;
      sum += coef[j] *at[j], ++j;
    
      sum += coef[j] *at[j], ++j;
      sum += coef[j] *at[j], ++j;
      sum += coef[j] *at[j], ++j;
      sum += coef[j] *at[j], ++j;
    
      sum += coef[j] *at[j], ++j;
      sum += coef[j] *at[j], ++j;
      sum += coef[j] *at[j], ++j;
      sum += coef[j] *at[j], ++j;
    
      assert(j == u150_l);
      output[i] = sum;
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

      float x = (float) (fraction << u150_2_b) * (1 / MULT32);

      float sum = 0;
      int j = 0;

      sum += ((coef[3*j] *x + coef[3*j + 1])*x + coef[3*j + 2])*at[j], ++j;
      sum += ((coef[3*j] *x + coef[3*j + 1])*x + coef[3*j + 2])*at[j], ++j;
      sum += ((coef[3*j] *x + coef[3*j + 1])*x + coef[3*j + 2])*at[j], ++j;
      sum += ((coef[3*j] *x + coef[3*j + 1])*x + coef[3*j + 2])*at[j], ++j;
    
      sum += ((coef[3*j] *x + coef[3*j + 1])*x + coef[3*j + 2])*at[j], ++j;
      sum += ((coef[3*j] *x + coef[3*j + 1])*x + coef[3*j + 2])*at[j], ++j;
      sum += ((coef[3*j] *x + coef[3*j + 1])*x + coef[3*j + 2])*at[j], ++j;
      sum += ((coef[3*j] *x + coef[3*j + 1])*x + coef[3*j + 2])*at[j], ++j;
    
      sum += ((coef[3*j] *x + coef[3*j + 1])*x + coef[3*j + 2])*at[j], ++j;
      sum += ((coef[3*j] *x + coef[3*j + 1])*x + coef[3*j + 2])*at[j], ++j;
      sum += ((coef[3*j] *x + coef[3*j + 1])*x + coef[3*j + 2])*at[j], ++j;
      sum += ((coef[3*j] *x + coef[3*j + 1])*x + coef[3*j + 2])*at[j], ++j;
    
      sum += ((coef[3*j] *x + coef[3*j + 1])*x + coef[3*j + 2])*at[j], ++j;
      sum += ((coef[3*j] *x + coef[3*j + 1])*x + coef[3*j + 2])*at[j], ++j;
      sum += ((coef[3*j] *x + coef[3*j + 1])*x + coef[3*j + 2])*at[j], ++j;
      sum += ((coef[3*j] *x + coef[3*j + 1])*x + coef[3*j + 2])*at[j], ++j;
    
      sum += ((coef[3*j] *x + coef[3*j + 1])*x + coef[3*j + 2])*at[j], ++j;
      sum += ((coef[3*j] *x + coef[3*j + 1])*x + coef[3*j + 2])*at[j], ++j;
      sum += ((coef[3*j] *x + coef[3*j + 1])*x + coef[3*j + 2])*at[j], ++j;
      sum += ((coef[3*j] *x + coef[3*j + 1])*x + coef[3*j + 2])*at[j], ++j;

      assert(j == u150_l);
      output[i] = sum;
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
