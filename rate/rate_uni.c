/* Copyright (c) 2011 lvqcl
 * This library is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation; either version 2.1 of the License, or (at
 * your option) any later version.
 *
 * This library is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser
 * General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this library; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
 */

#include <assert.h>
#include <intrin.h>
#include "rate_i.h"
#include "sox_i.h"

int haveSSE = -1;
int haveSSE3 = -1;

int initialized = 0;

RR_handle* RR_open(const RR_config* config, int nchannels, int* error)
{
	RR_handle* h = NULL;

	assert(initialized);
	if (!initialized)
	{
		if (error) *error = RR_EXTUNINIT;
		return NULL;
	}

	if (config->quality == RR_best)
	{
		if (haveSSE3 == 1)
			h = RR_it_SSE3(config, nchannels);
		else
			h = RR_it_double(config, nchannels);
	}
	else
	{
		if (haveSSE == 1)
			h = RR_it_SSE(config, nchannels);
		else
			h = RR_it_float(config, nchannels);
	}

	if (error)
	{
		if (h == NULL) *error = RR_ENOMEM;
		else *error = RR_OK;
	}

	return h;
}

void RR_close(RR_handle** h)
{
	assert(h);
	if (h == NULL || *h == NULL) return;

	((RR_vtable2*)*h)->close(*h); // <= lsx_free(*h) there
	*h = NULL;
}

const char* RR_strerror(int error)
{
	switch(error)
	{
	case RR_OK:
		return "OK";
	case RR_ENOMEM:
		return "Not enough memory";
	case RR_INTERNAL:
		return "Internal error";
	case RR_NULLHANDLE:
		return "NULL handle";
	case RR_RATEERROR:
		return "Error in rate() functions";
	case RR_EXTUNINIT:
		return "Externals not initialized";
	default:
		return "Other error";
	}
}



static void init_SSEflags()
{
	int info[4];
	__cpuid(info, 1);

	if (info[2] & 1) /* SSE3 */
		haveSSE3 = 1;
	else 
		haveSSE3 = 0;

	if (info[3] & (1 << 25)) /* SSE */
		haveSSE = 1;
	else 
		haveSSE = 0;

	//if (info[3] & (1 << 26)) {...} /* SSE2 */

	assert(haveSSE3 != -1 && haveSSE != -1);
}

static int init_fft4g()
{
	int i;
	int res = 0;

	lsx_rdft = (haveSSE3 == 1) ? lsx_rdft_SSE3 : lsx_rdft_generic;

	for(i=1; i<=16; i++)
	{
		if (lsx_rdft_init(1<<i, &fftx[i]) < 0)
			res = -1;
	}

	if (res == -1)
	{
		for(i=1; i<=16; i++)
			lsx_rdft_close(&fftx[i]);
	}

	return res;
}

static int init_ffmpeg()
{
	int i;
	int res = 0;
	for (i=4; i<=16; i++)
	{
		if (ff_rdft_init(&ff_fwd[i], i,  DFT_R2C, haveSSE == 1) < 0)
			res = -1;
		if (ff_rdft_init(&ff_bkd[i], i, IDFT_C2R, haveSSE == 1) < 0)
			res = -1;
	}

	if (res == -1)
	{
		for (i=4; i<=16; i++)
		{
			ff_rdft_end(&ff_fwd[i]);
			ff_rdft_end(&ff_bkd[i]);
		}
	}

	return res;
}

static void close_fft4g()
{
	int i;

	for(i=1; i<=16; i++)
		lsx_rdft_close(&fftx[i]);
}

static void close_ffmpeg()
{
	int i;

	for (i=4; i<=16; i++)
	{
		ff_rdft_end(&ff_fwd[i]);
		ff_rdft_end(&ff_bkd[i]);
	}
}

int init_ratelib(void (*h)(void)) // not thread-safe
{
	int res = 0;

	initialized = 0;

	if (h==NULL) return -1;
	set_malloc_error_handler(h);
	init_SSEflags();
	if (init_fft4g() < 0) return -1;
	if (init_ffmpeg() < 0) { close_fft4g(); return -1; }

	initialized = 1;
	return 0;
}

void close_ratelib() // not thread-safe
{
	close_fft4g();
	close_ffmpeg();

	initialized = 0;
}

void (*lsx_rdft)(int isgn, double *a, FFTcontext z);
FFTcontext fftx[17];

RDFTContext ff_fwd[17];
RDFTContext ff_bkd[17];
