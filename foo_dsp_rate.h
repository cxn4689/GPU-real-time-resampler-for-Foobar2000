/* Copyright (c) 2008-10 lvqcl.  All rights reserved.
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

#ifndef FOODSPRATE_H
#define FOODSPRATE_H

//#include "rate/ratelib.h"
#include "dsp_config.h"

class dsp_rate : public dsp_impl_base_t<dsp_v2>
{
	RateConfig cfg_;
	unsigned out_rate_;

	audio_sample* buffer_;

	size_t in_samples_accum_, out_samples_gen_accum_;

	unsigned sample_rate_;
	unsigned channel_count_;
	unsigned channel_map_;
	unsigned buf_channels_;

	//RR_handle* rate_state_;

	void init();
	bool set_data(const dsp_preset & p_data);
 
public:
	dsp_rate();
	dsp_rate(const t_dsp_rate_params& params);
	~dsp_rate();

private:
	void reinit(unsigned sample_rate, unsigned channel_count, unsigned channel_map);
	void flushwrite(bool close);

protected:
	virtual void on_endoftrack(abort_callback & p_abort);
	virtual void on_endofplayback(abort_callback & p_abort);
	virtual bool on_chunk(audio_chunk * chunk, abort_callback & p_abort);

public:
	virtual void flush();
	virtual double get_latency();
	virtual bool need_track_change_mark() { return false; }

	//my test linear
	double n_ratio;
	double n_last[2];
	int n_in_sample;
	int n_out_sample;
	//my test linear end
};

#endif
