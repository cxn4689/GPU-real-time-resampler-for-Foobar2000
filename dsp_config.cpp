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

#define _CRT_SECURE_NO_WARNINGS 1
#include "dsp_config.h"

const int STRMAXLEN = 20;

static const long samplerates[] = {8000, 11025, 16000, 22050, 24000, 32000, 44100, 48000, 64000, 88200, 96000, 176400, 192000};

static const TCHAR* srates_z[] = {0, 0,
					_T("Upsample x2"), _T("Downsample x2"),
					_T("Downsample x4"), _T("Upsample x4")};



static const TCHAR* samplerates_s[] = {_T("8000"), _T("11025"), _T("16000"), _T("22050"), _T("24000"), _T("32000"),
					_T("44100"), _T("48000"), _T("64000"), _T("88200"), _T("96000"), _T("176400"), _T("192000"),
					srates_z[2], srates_z[5], srates_z[3], srates_z[4]};
static const TCHAR* qualities[] = {_T("Best"), _T("Good"), _T("Normal")};

t_dsp_rate_params::t_dsp_rate_params()
{
	cfg_.outRate        = 44100;
	cfg_.quality        = Qbest;

	cfg_.allow_aliasing = t_int32(statics::aliasing());
	cfg_.passband10     = t_int32(statics::passband()*10);
	cfg_.phase          = t_int32(statics::phase());
}

void t_dsp_rate_params::get_rateconfig(RateConfig& cfg) const
{
	cfg_.CloneTo(cfg);
}

const int NParam = 5;
bool t_dsp_rate_params::set_data(const dsp_preset& p_data)
{
	t_int32 temp[NParam];
	if (p_data.get_owner() != g_get_guid()) return false;
	if (p_data.get_data_size() != sizeof(temp))
	{
		//from previous versions? -- in future
		return false;
	}
	memcpy(temp, p_data.get_data(), sizeof(temp));

	for(int i=0; i<NParam; i++)
		byte_order::order_le_to_native_t(temp[i]);
		
	cfg_.outRate        = temp[0];
	cfg_.quality        = temp[1];
	cfg_.allow_aliasing = temp[2];
	cfg_.passband10     = temp[3];
	cfg_.phase          = temp[4];

	if (cfg_.quality > Qlast)  // if someone uses prev. versions and Low Quality
		cfg_.quality = Qlast;

	if (cfg_.passband10 == 0 || cfg_.passband10 == 1) //bool steepness -- from previous version (0.3.2.  ha!)
	{
		cfg_.passband10 = cfg_.passband10?maxPassband10:defPassband10;
		cfg_.phase *= 25;
	}
	return true;
}

bool t_dsp_rate_params::get_data(dsp_preset& p_data) const
{
	t_int32 temp[NParam];
	temp[0] = cfg_.outRate;
	temp[1] = cfg_.quality;
	temp[2] = cfg_.allow_aliasing;
	temp[3] = cfg_.passband10;
	temp[4] = cfg_.phase;

	p_data.set_owner(g_get_guid());
	for(int i=0; i<NParam; i++)
		byte_order::order_native_to_le_t(temp[i]);
	p_data.set_data(temp, sizeof(temp));
	return true;
}

void t_dsp_rate_params::tset_outRate(const TCHAR* rate)
{
	for(int i = -TRupsample_2; i<= -TRupsample_4; i++)
	{
		if (_tcscmp(rate, srates_z[i]) == 0) {
		cfg_.outRate = -i; return; }
	}

	set_outRate(_ttoi(rate));
}

const TCHAR* t_dsp_rate_params::toutRateStr(TCHAR* buf) const
{
	if (cfg_.outRate < 0) //cfg_.outRate = [-5...-2].
	{
		assert(cfg_.outRate >= -5 && cfg_.outRate <= -2);
		_tcsncpy_s(buf, 30, srates_z[-cfg_.outRate], _TRUNCATE);
		return buf;
	}

	_itot_s(cfg_.outRate, buf, STRMAXLEN-1, 10);
	return buf;
}


void dialog_dsp_rate::update_phresponseText(TCHAR* buf, int sz, CStatic phresponseText, int val)
{
	if (val == Pminimum) {
		phresponseText.SetWindowTextW(_T(" 0 % (minimum)"));
	}
	else if (val == Plinear) {
		phresponseText.SetWindowTextW(_T("50 % (linear)"));
	}
	else if (val == Pmaximum) {
		phresponseText.SetWindowTextW(_T("100 % (maximum)"));
	}
	else {
		_sntprintf(buf, sz, _T("%2d %%"), val);
		phresponseText.SetWindowTextW(buf);
	}
}


BOOL dialog_dsp_rate::OnInitDialog(CWindow wndFocus, LPARAM lInitParam)
{
	TCHAR buf[30];

	CComboBox rateCombo(GetDlgItem(IDC_RATE));
	CComboBox qualCombo(GetDlgItem(IDC_QUALITY));
	CButton aliasCheckbox(GetDlgItem(IDC_ALIASING));
	CTrackBarCtrl passbandSlider(GetDlgItem(IDC_PASSBAND));
	CTrackBarCtrl phresponseSlider(GetDlgItem(IDC_PHRESPONSE));
	
	CStatic passbandText(GetDlgItem(IDC_BWTEXT));
	CStatic phresponseText(GetDlgItem(IDC_PHASETEXT));

	rateCombo.LimitText(12);
	for(int i=0; i<sizeof(samplerates_s)/sizeof(samplerates_s[0]); i++)
		rateCombo.AddString((samplerates_s[i]));
	params_.toutRateStr(buf);
	for(int i=0; i<sizeof(samplerates_s)/sizeof(samplerates_s[0]); i++)
	{
		if (!_tcscmp(samplerates_s[i], buf)) { rateCombo.SetCurSel(i); break; }
	}
	rateCombo.SetWindowTextW(buf);

	for(int i=0; i<sizeof(qualities)/sizeof(qualities[0]); i++)
		qualCombo.AddString((qualities[i]));
	qualCombo.SetCurSel(params_.quality());

	aliasCheckbox.SetCheck(params_.aliasing()?BST_CHECKED:BST_UNCHECKED);

	passbandSlider.SetRange(minPassband10, maxPassband10);
	passbandSlider.SetPos(params_.passband10());
	passbandSlider.SetTicFreq(tickPassband10);
	passbandSlider.SetPageSize(tickPassband10);

	phresponseSlider.SetRange(Pminimum, Plinear);
	phresponseSlider.SetPos(params_.phase());
	phresponseSlider.SetTicFreq(5);
	phresponseSlider.SetPageSize(5);

	_sntprintf(buf, sizeof(buf)/sizeof(TCHAR), _T("%2.1f %%"), (float(params_.passband10())/10.0));
	passbandText.SetWindowTextW(buf);

	update_phresponseText(buf, sizeof(buf)/sizeof(TCHAR), phresponseText, params_.phase());
				
	return 0;
}

void dialog_dsp_rate::OnHScroll(UINT nSBCode, UINT nPos, CScrollBar pScrollBar)
{
	TCHAR buf[30];

	CTrackBarCtrl passbandSlider(GetDlgItem(IDC_PASSBAND));
	CTrackBarCtrl phresponseSlider(GetDlgItem(IDC_PHRESPONSE));
	
	if (HWND(pScrollBar) == HWND(passbandSlider))
	{
		CStatic passbandText(GetDlgItem(IDC_BWTEXT));
		_sntprintf(buf, sizeof(buf)/sizeof(TCHAR), _T("%2.1f %%"), float(passbandSlider.GetPos()/10.0));
		passbandText.SetWindowTextW(buf);
	}
	if (HWND(pScrollBar) == HWND(phresponseSlider))
	{
		CStatic phresponseText(GetDlgItem(IDC_PHASETEXT));
		update_phresponseText(buf, sizeof(buf)/sizeof(TCHAR), phresponseText, phresponseSlider.GetPos());
	}
}

void dialog_dsp_rate::OnCommand(UINT uNotifyCode, int nID, CWindow wndCtl)
{
	switch (nID)
	{
	case IDOK:
		{
			TCHAR str[STRMAXLEN];
			CComboBox rateCombo(GetDlgItem(IDC_RATE));
			CComboBox qualCombo(GetDlgItem(IDC_QUALITY));
			CButton aliasCheckbox(GetDlgItem(IDC_ALIASING));
			CTrackBarCtrl passbandSlider(GetDlgItem(IDC_PASSBAND));
			CTrackBarCtrl phresponseSlider(GetDlgItem(IDC_PHRESPONSE));

			rateCombo.GetWindowTextW(str, STRMAXLEN);
			params_.tset_outRate(str);

			params_.set_quality(qualCombo.GetCurSel());
			params_.set_aliasing(aliasCheckbox.GetCheck());
			params_.set_passband10(passbandSlider.GetPos());
			params_.set_phase(phresponseSlider.GetPos());

			// Data (potentially) changed
			EndDialog(IDOK);
		}
		break;

		case IDCANCEL:
		{
			// Data not changed
			EndDialog(0);
		}
		break;

		/*case IDC_DEFAULTS:
		{
			TCHAR buf[30];
			//CComboBox rateCombo(GetDlgItem(IDC_RATE));
			CComboBox qualCombo(GetDlgItem(IDC_QUALITY));
			CButton aliasCheckbox(GetDlgItem(IDC_ALIASING));
			CTrackBarCtrl passbandSlider(GetDlgItem(IDC_PASSBAND));
			CTrackBarCtrl phresponseSlider(GetDlgItem(IDC_PHRESPONSE));
			CStatic passbandText(GetDlgItem(IDC_BWTEXT));
			CStatic phresponseText(GetDlgItem(IDC_PHASETEXT));
			
			qualCombo.SetCurSel(Qbest);

			t_dsp_rate_params params;

			aliasCheckbox.SetCheck(params.aliasing());

			passbandSlider.SetPos(params.passband10());
			_sntprintf(buf, sizeof(buf)/sizeof(TCHAR), _T("%2.1f %%"), float(params.passband10()/10.0));
			passbandText.SetWindowTextW(buf);

			phresponseSlider.SetPos(params.phase());
			update_phresponseText(buf, sizeof(buf)/sizeof(TCHAR), phresponseText, (int)params.phase());
		}
		break;*/
	}
}
