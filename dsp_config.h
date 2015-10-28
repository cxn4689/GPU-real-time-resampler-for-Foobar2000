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

#ifndef DSPRATECONFIG_H
#define DSPRATECONFIG_H

//#include "../ATLhelpers/ATLhelpers.h"
//#include "../SDK/foobar2000.h"
#include "ATLhelpers.h"
#include "foobar2000.h"

#include "resource.h"



#define RESAMPLER_MAJOR_VERSION  0 /* Major version number */
#define RESAMPLER_MINOR_VERSION  7 /* Minor version number */
#define RESAMPLER_PATCH_VERSION  6/* Patch level */
#define RESAMPLER_BUILD_VERSION  0 /* Build number */

#define RESAMPLER_VERSION (RESAMPLER_MAJOR_VERSION*0x10000 + RESAMPLER_MINOR_VERSION*0x100 + RESAMPLER_PATCH_VERSION)

#ifndef STR
#define STR__(x)  #x
#define STR(x)    STR__(x)
#endif

#if defined(RESAMPLER_ALPHA_VERSION)
#define RESAMPLER_VERSION_STR STR(RESAMPLER_MAJOR_VERSION)"."STR(RESAMPLER_MINOR_VERSION)"."STR(RESAMPLER_PATCH_VERSION)" alpha "STR(RESAMPLER_ALPHA_VERSION)
#elif defined(RESAMPLER_BETA_VERSION)
#define RESAMPLER_VERSION_STR STR(RESAMPLER_MAJOR_VERSION)"."STR(RESAMPLER_MINOR_VERSION)"."STR(RESAMPLER_PATCH_VERSION)" beta "STR(RESAMPLER_BETA_VERSION)
#else
#define RESAMPLER_VERSION_STR STR(RESAMPLER_MAJOR_VERSION)"."STR(RESAMPLER_MINOR_VERSION)"."STR(RESAMPLER_PATCH_VERSION)
#endif

enum QUAL
{
	Qbest = 0,
	Qgood = 1,
	Qfast = 2,
	Qlast = 2
};

const int Pminimum = 0, Plinear = 50, Pmaximum = 100;
const int defPassband10 = 950;
const int minPassband10 = 900;
const int maxPassband10 = 990;
const int tickPassband10 = 10;

enum TargetRates
{
	TRuninit = 0,
	TRerror = -1,

	TRupsample_2 = -2,
	TRupsample_4 = -5,

	TRdnsample_2 = -3,
	TRdnsample_4 = -4
};

class RateConfig
{
public:
	t_int32 outRate;
	t_int32 quality;
	t_int32 allow_aliasing;
	t_int32 passband10;
	t_int32 phase;

	unsigned realrate(unsigned in_samplerate) const
	{
		if (outRate > 0) return outRate;
		switch(outRate)
		{
		case TRupsample_2:
			return in_samplerate*2;
		case TRupsample_4:
			return in_samplerate*4;
		case TRdnsample_2:
			return in_samplerate/2;
		case TRdnsample_4:
			return in_samplerate/4;
		}
		return in_samplerate; //something strange
	}

	bool is_no_resample(unsigned in_samplerate) const
	{
		if (in_samplerate == outRate) return true;
		return false;
	}

	void CloneTo(RateConfig& to) const
	{
		to.outRate = outRate;
		to.quality = quality;
		to.allow_aliasing = allow_aliasing;
		to.passband10 = passband10;
		to.phase = phase;
	}
};

class t_dsp_rate_params
{
	RateConfig cfg_;
public:	
	t_dsp_rate_params();
	t_dsp_rate_params(const RateConfig& cfg):cfg_(cfg){}
	void get_rateconfig(RateConfig& cfg) const;

#if 1
	static const GUID &g_get_guid()
	{	// {52EEE681-6ADx-4378-B1B1-480CA9FF3014}
		static const GUID guid = { 0x52eee681, 0x6ade, 0x4378, { 0xb1, 0xb1, 0x48, 0x0c, 0xa9, 0xff, 0x30, 0x14 } };
		return guid;
	}
#else
	static const GUID &g_get_guid() //ALPHA
	{	// {61B9849B-595F-4852-9394-E7D3213CF234}
		static const GUID guid = { 0x61b9849b, 0x595f, 0x4852, { 0x93, 0x94, 0xe7, 0xd3, 0x21, 0x3c, 0xf2, 0x34 } };
		return guid;
	}
#endif

	//store data from preset
	bool set_data(const dsp_preset& p_data);
	//put data to preset
	bool get_data(dsp_preset& p_data) const;

	//inline t_int32 outRate(...) const { return cfg_.outRate; }
	inline t_int32 quality() const { return cfg_.quality; }
	inline t_int32 aliasing() const { return cfg_.allow_aliasing; }
	inline t_int32 passband10() const { return cfg_.passband10; }
	inline t_int32 phase() const { return cfg_.phase; }
	const TCHAR* toutRateStr(TCHAR*) const;

	void tset_outRate(const TCHAR* rate);
	inline void set_outRate(t_int32 rate)
	{
		if (rate < audio_chunk::sample_rate_min) rate = audio_chunk::sample_rate_min;
		else if (rate > audio_chunk::sample_rate_max) rate = audio_chunk::sample_rate_max;
		cfg_.outRate = rate;
	}
	inline void set_quality(t_int32 q){ cfg_.quality = q; }
	inline void set_aliasing(t_int32 a){ cfg_.allow_aliasing = a; }
	inline void set_passband10(t_int32 pb){ cfg_.passband10 = pb; }
	inline void set_phase(t_int32 ph){ cfg_.phase = ph; }
};

class dialog_dsp_rate : public CDialogImpl<dialog_dsp_rate>
{
public:
	enum
	{
		IDD = IDD_CONFIG
	};
 
	BEGIN_MSG_MAP_EX(dialog_dsp_rate)
		MSG_WM_INITDIALOG(OnInitDialog)
		MSG_WM_HSCROLL(OnHScroll)
		MSG_WM_COMMAND(OnCommand)
	END_MSG_MAP()

private:
	t_dsp_rate_params& params_;

	void update_phresponseText(TCHAR* buf, int sz, CStatic phresponseText, int val);

public:
	dialog_dsp_rate(t_dsp_rate_params& p_params) : params_(p_params) {}

protected:
	BOOL OnInitDialog(CWindow wndFocus, LPARAM lInitParam);
	void OnHScroll(UINT nSBCode, UINT nPos, CScrollBar pScrollBar);
	void OnCommand(UINT uNotifyCode, int nID, CWindow wndCtl);
};

class statics
{
	static advconfig_integer_factory*  p_int_soxr_passband;
	static advconfig_checkbox_factory* p_chk_soxr_aliasing;
	static advconfig_integer_factory*  p_int_soxr_phase;
	static advconfig_integer_factory*  p_int_soxr_priority;

public:
	static bool aliasing() { return *p_chk_soxr_aliasing; }
	static t_uint64 passband() { return *p_int_soxr_passband; }
	static t_uint64 phase() { return *p_int_soxr_phase; }
	static t_uint64 priority() { return *p_int_soxr_priority; }
};

#endif
