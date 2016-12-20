#include "stdafx.h"
#include<memory.h>
#include<math.h>
#include "trunk_sinc.h"


//float sinc_table [_table_size] ;//sinc_pi

#define pi 3.14159265358979323846264338327950288 

double delta = 1.0 / prec_;
double windows_(double t,double nT){//bn sequence length, blackman window,2*n
	if(t<=nT&&t>=-nT) return 0.42+0.5*cos(pi*t/nT)+0.08*cos(2*pi*t/nT);
	else return 0;
}
double sinc_(double t,double ratio){ //ratio=ori/des
	if(ratio<=1.0){
		if(fabs(t)<(double)pow((double)10,(double)-30)) return 1.0;
		double x=sin(t*pi);
		double y=(pi*t);
		if (x/y>1.0) return 1.0;
		return x/y;
	}
	else{//pre filter to prevent alias
		double max=1.0/ratio;
		if(fabs(t)<(double)pow((double)10,(double)-30)) return max;
		double x=sin(t*pi*max);
		double y=(pi*t);
		if (x/y>max) return max;
		else return x/y;
	}
}
void generate_table(int sample_half,float factor/*org_sr/des_sr*/,float *sinc_table)//30 points per zero crossing
{//zero crossing number
	//normalized Ts = 1, Fs = 1;
	//interval 1/pi/30
	double delta = 1.0  / prec_;
	int n = (int)((double)sample_half * prec_); //length
	if(n >= _table_size) return;
	double index = 0.0;
	for(int i = 0; i < n; i++)
	{
		sinc_table[i] = sinc_(index,factor) * windows_(index,sample_half);
		index += delta;
	}
}

float lookup(float t,float *sinc_table)//normallized, real t
{
	//normalized Ts = 1, Fs = 1;
	//interval 1/pi/30
	/*unsigned int exp_plus = 7 <<23;
	*t_ += exp_plus;//t*128*/
	//unsigned int *t_ = reinterpret_cast<unsigned int *>(&t);
	//*t_ &= 0x7fffffff;//convert to positive number
	t = fabs(t);//only get from half
	t *= prec_;
	int index = (int)(t);//get index
	if(index > _table_size) return 0.0;
	float sap = sinc_table[index];
	//unsigned exp_ = (*t_) & 0xff;
	//unsigned fra_ = (*t_) & 0x7fffff;
	float dt = (t - (float)index);//find a quick replace algorithm fraction x.M - x.0
	float dt2 = sinc_table[index+1] - sap;
	return sap + dt*dt2;

}