/* Copyright (c) 2008-11 lvqcl.  All rights reserved.
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

#include "foo_dsp_rate.h"
#include "main.h"
#include "SRC_Stream.h"
#include "cl_interface.h"

extern "C" {
void throw_new_exception()
{
	#pragma warning(disable:4297)
	throw std::bad_alloc();
}

}

const int buf_samples = 8192;

void dsp_rate::init()
{
	t_dsp_rate_params params;
	params.get_rateconfig(cfg_);

	//rate_state_ = NULL;

	sample_rate_ = 0;
	out_rate_ = 0;
	channel_count_ = 0;
	channel_map_ = 0;
	
	buf_channels_ = 0;

	buffer_ = NULL; 
	
	in_samples_accum_ = out_samples_gen_accum_ = 0;

	//my linear
	this->n_ratio=0;
	this->n_last[0]=this->n_last[1]=0;
	this->n_in_sample=0;
	this->n_out_sample=0;
	//my sinc
	__my_reset();
}

dsp_rate::dsp_rate()
{
	init();
	platform_initial();//my opencl
}

dsp_rate::dsp_rate(const t_dsp_rate_params& params)
{
	init();
	params.get_rateconfig(cfg_);
	platform_initial();//my opencl
}

/*dsp_rate::dsp_rate(const dsp_preset& p_data)
{
	init();
	set_data(p_data);
}*/

bool dsp_rate::set_data(const dsp_preset & p_data)
{
	t_dsp_rate_params params;
	if (!params.set_data(p_data)) return false;
	params.get_rateconfig(cfg_);
	return true;
}

dsp_rate::~dsp_rate()
{
	//RR_close(&rate_state_); //rate_state_ == NULL
	if(buffer_!=NULL) delete[] buffer_;
	buffer_ = NULL;
	close_platform();
	__my_clear_all();//fix memory leaking
}

void dsp_rate::on_endoftrack(abort_callback & p_abort) 
{
	init(); 
	if(buffer_!=NULL) delete []buffer_;//my test
	flushwrite(false); 
}

void dsp_rate::on_endofplayback(abort_callback & p_abort) 
{
	init();
	if(buffer_!=NULL) delete []buffer_;//my test
	flushwrite(true);
}

void dsp_rate::reinit(unsigned sample_rate, unsigned channel_count, unsigned channel_map)
{
	int rate_error;

	if (buf_channels_ < channel_count)
	{
		if(buffer_!=NULL) delete[] buffer_;
		buf_channels_ = channel_count;
		buffer_ = new audio_sample[buf_samples*buf_channels_];
	}

	out_rate_ = cfg_.realrate(sample_rate);
	//RR_config c = {sample_rate, out_rate_, cfg_.phase, double(cfg_.passband10)/10.0, cfg_.allow_aliasing?true:false, (RR_quality)(RR_best - cfg_.quality) };

	//rate_state_ = RR_open(&c, channel_count, &rate_error);
	//if (rate_state_ == NULL) { console::error(RR_strerror(rate_error)); throw exception_out_of_resources(); }

	channel_count_ = channel_count; channel_map_ = channel_map; sample_rate_ = sample_rate;
	in_samples_accum_ = out_samples_gen_accum_ = 0;

	//my test
	this->n_ratio=0;
	this->n_last[0]=this->n_last[1]=0;
	this->n_in_sample=0;
	this->n_out_sample=0;
	//my test end
	__my_reset();
	__my_set_parameter(this->sample_rate_,this->out_rate_);
}

//for test beg
double linear_interpolation(double ratio,float nL,float nR,float vL,float vR,float nI){
	double dnL=nL,dnR=nR,dvL=vL,dvR=vR;
	//linear interpolation
	if(ratio>1.0) return 0;//down sampling not support currently
	double dvI=dvL*(dnR-nI)+dvR*(nI-dnL);
	double vI=(double)dvI;
	return vI;
}
//for test end

bool dsp_rate::on_chunk(audio_chunk * chunk, abort_callback & p_abort)
{
	size_t in_samples_used=0, out_samples_gen=0;
	int rate_error;

	unsigned channel_count = chunk->get_channels();
	unsigned channel_map = chunk->get_channel_config();
	unsigned sample_rate = chunk->get_sample_rate();
	t_size sample_count = chunk->get_sample_count();
	audio_sample * current = chunk->get_data();

	//my test
	this->n_ratio=44100.0/48000.0;

	if (/*rate_state_*/1 == NULL) //uninitialized state
	{
		if (cfg_.is_no_resample(sample_rate)) return true;
		reinit(sample_rate, channel_count, channel_map);
		//if (rate_state_ == NULL) { console::error("SoX resampler: cannot initialize rate library"); return true; } //to avoid crashes etc.
	}
	else if ((channel_count_ != channel_count) || (channel_map_ != channel_map) || (sample_rate_ != sample_rate))
	{	// number of channels or samplerate has changed - reinitialize
		flushwrite(true); //here channel_count_, channel_map_ and sample_rate_ must have old values
		//rate_state_ == NULL here
		if (cfg_.is_no_resample(sample_rate)) return true;
		reinit(sample_rate, channel_count, channel_map);
		//if (rate_state_ == NULL) { console::error("SoX resampler: cannot initialize rate library"); return true; } //to avoid crashes etc.
	}

	int beg_index=in_samples_accum_;//my test
	double delta=(double)this->n_out_sample*this->n_ratio-(double)beg_index;//my test

	in_samples_accum_ += sample_count;

	//my interpolation
	int o_size=0;
	int f_in=sample_count;
	//buffer 4096*channel
	do{
		o_size=__my_get_data_and_src(current,f_in*2,this->buffer_,o_size);
		if(f_in!=0) f_in=0;//flush buffer
		if(o_size!=0){
			audio_chunk *out=insert_chunk(/*o_size*2*/);
			out->set_data(buffer_,o_size,2,out_rate_);//out_rate_ == output sample rate
			out_samples_gen_accum_+=o_size;
		}
	}while(o_size!=0);

	/*in_samples_used=f_in;//used
	current += in_samples_used * channel_count_;*///fix memory leaking? audio_sample shoud += to release the audio input buffer

	//out_samples_gen_accum_+=sample_count;

	while (in_samples_accum_ > sample_rate_ && out_samples_gen_accum_ > out_rate_)
	{
		in_samples_accum_ -= sample_rate_;
		out_samples_gen_accum_ -= out_rate_;
	}
	return false;
}

void dsp_rate::flush()
{
	int rate_error;
	if (!0/*rate_state_*/) return;
	rate_error = 0;//RR_reset(rate_state_);
	//if (rate_error) { console::error(RR_strerror(rate_error)); throw exception_out_of_resources(); }
	in_samples_accum_ = out_samples_gen_accum_ = 0;
	return;
}

void dsp_rate::flushwrite(bool close)
{
	size_t out_samples_gen=0;
	int rate_error;

	//if (!rate_state_) return;

	rate_error = 0;//RR_drain(rate_state_);
	//if (rate_error) { console::error(RR_strerror(rate_error)); throw exception_out_of_resources(); }
	while(0)
	{
		//rate_error = RR_flow(rate_state_, NULL, buffer_, 0, buf_samples, NULL, &out_samples_gen);
		//if (rate_error) { console::error(RR_strerror(rate_error)); throw exception_out_of_resources(); }
		if (out_samples_gen == 0) break;
		out_samples_gen_accum_ += out_samples_gen;
			
		audio_chunk * out = insert_chunk(/*out_samples_gen*channel_count_*/);
		out->set_data(buffer_, out_samples_gen, channel_count_, out_rate_, channel_map_);
	}
	if (close){}
		//RR_close(&rate_state_);
	else
	{
		//rate_error = RR_reset(rate_state_);
		//if (rate_error) { console::error(RR_strerror(rate_error)); throw exception_out_of_resources(); }
	}
	in_samples_accum_ = out_samples_gen_accum_ = 0;
	return;
}

double dsp_rate::get_latency()
{
	if (sample_rate_ && out_rate_)
		return double(in_samples_accum_)/double(sample_rate_) - double(out_samples_gen_accum_)/double(out_rate_);
	else return 0;
	// warning: sample_rate_ and out_rate_ can be from previous track (?),
	// but in_samples_accum_ == out_samples_gen_accum_ == 0,  so it's ok.
}
