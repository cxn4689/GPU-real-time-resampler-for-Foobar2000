#ifdef _intel_cpu_
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#endif

#define _data_type_ float
#define buf_fix_size 4096

int linear_interpolation_GPU(double ratio,int nL,int nR,int vL,int vR,double nI){
	double dnL=nL,dnR=nR,dvL=vL,dvR=vR;
	//linear interpolation
	if(ratio>1.0) return 0;//down sampling not support currently
	double dvI=dvL*(dnR-nI)+dvR*(nI-dnL);
	int vI=(int)dvI;
	return vI;
}
#define pi 3.14159265358979323846264338327950288f

double windows_GPU(double t,double nT){//bn sequence length, blackman window,2*n
	if(t<=nT&&t>=-nT) return 0.42+0.5*cos(pi*t/nT)+0.08*cos(2*pi*t/nT);
	else return 0;
}

double sinc_GPU(double t,double fs,double ratio){ //ratio=ori/des
	if(ratio<=1.0){
		if(fabs(t)<(double)pow((double)10,(double)-30)) return 1.0;
		double x=sin(fs*t*pi);
		double y=(pi*fs*t);
		if (x/y>1.0) return 1.0;
		return x/y;
	}
	else{//pre filter to prevent alias
		double max=1.0/ratio;
		if(fabs(t)<(double)pow((double)10,(double)-30)) return max;
		double x=sin(fs*t*pi*max);
		double y=(pi*t*fs);
		if (x/y>max) return max;
		else return x/y;
	}
}
#define prec_ 512
#define _table_size 65537


inline double lookup(double t,__global _data_type_ *sinc_table,_data_type_ ratio)//normallized, real t
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
	double sap = (double)sinc_table[index];
	//unsigned exp_ = (*t_) & 0xff;
	//unsigned fra_ = (*t_) & 0x7fffff;
	double dt = (t - (double)index);//find a quick replace algorithm fraction x.M - x.0
	double dt2 = (double)sinc_table[index+1] - sap;
	double ret_v = sap + dt*dt2;
	double _1_ratio = (double)(1.0/ratio);
	ret_v = (ratio > 1.0 && ret_v > _1_ratio)? _1_ratio : ret_v;
	return ret_v;

}

__kernel void interpolation(__global _data_type_ *buffer_process,int buf_size,__global _data_type_ *out_buf,int out_size,
int beg_index,int out_samples,int chn_num,int ori_sr,int des_sr,int w_hn, __global _data_type_ *sinc_table_)
{
	__global _data_type_ *sinc_table = (__global _data_type_ *)sinc_table_;
	int ix=get_global_id(0);//destination frame id
	if(ix>out_size) return;
	double ratio=(double)ori_sr/(double)des_sr;
	int ori_sample_rate=ori_sr;
	int des_sample_rate=des_sr;
	double rs_p_begin=(double)(out_samples+ix)*ratio;//resample point
	double ti=rs_p_begin;//(double)(out_samples+ix)/(double)des_sr;
	int nL=(int)rs_p_begin;
	int nR=nL+1;//caculating left and right
	//cacaulating actuall index
	int a_nL=nL-beg_index;
	int a_nR=nR-beg_index;
	if(a_nL>buf_fix_size-w_hn){return;}
	_data_type_ vL=(a_nL>=0&&a_nL<=buf_fix_size-1)?buffer_process[a_nL*chn_num]:0;
	_data_type_ vR=(a_nR<=buf_fix_size-1&&a_nR>=0)?buffer_process[a_nR*chn_num]:0;
	/*for(int i=0;i<10;i++)
	{
		if(buffer_process[i*2]!=0){
			printf("%i:%i\n",i,buffer_process[i*2]);
		}
	}*/
	//interpolation
	{
		double nIVL=0;
		//double ti=rs_p_begin;//insert point(Fs')
		for(int i=0;i<w_hn;i++){//right wing
			//window size 32
			int nId=a_nR+i;
			_data_type_ vL=(nId>=0&&nId<=buf_fix_size-1)?buffer_process[nId*chn_num]:0;
			int fabs_Index=nR+i;
			int index_f=fabs_Index;
			nIVL+=(double)((double)vL*lookup(ti-(double)index_f,sinc_table,ratio));//sinc_GPU(fabs(ti-(double)index_f),1,ratio)*windows_GPU(fabs(ti-(double)index_f),(double)w_hn));
		}
		for(int i=0;i<w_hn;i++){//right wing
			//window size 32
			int nId=a_nL-i;
			_data_type_ vL=(nId>=0&&nId<=buf_fix_size-1)?buffer_process[nId*chn_num]:0;
			int fabs_Index=nL-i;
			int index_f=fabs_Index;
			nIVL+=(double)((double)vL*lookup(ti-(double)index_f,sinc_table,ratio));//sinc_GPU(fabs(ti-(double)index_f),1,ratio)*windows_GPU(fabs(ti-(double)index_f),(double)w_hn));
		}
		out_buf[ix*chn_num]=(_data_type_)nIVL;
	}//left channel
	{
		double nIVL=0;
		//double ti=rs_p_begin;//insert point(Fs')
		for(int i=0;i<w_hn;i++){//right wing
			//window size 32
			int nId=a_nR+i;
			_data_type_ vL=(nId>=0&&nId<=buf_fix_size-1)?buffer_process[nId*chn_num+1]:0;
			int fabs_Index=nR+i;
			int index_f=fabs_Index;
			nIVL+=(double)((double)vL*lookup(ti-(double)index_f,sinc_table,ratio));//sinc_GPU(fabs(ti-(double)index_f),1,ratio)*windows_GPU(fabs(ti-(double)index_f),(double)w_hn));
		}
		for(int i=0;i<w_hn;i++){//right wing
			//window size 32
			int nId=a_nL-i;
			_data_type_ vL=(nId>=0&&nId<=buf_fix_size-1)?buffer_process[nId*chn_num+1]:0;
			int fabs_Index=nL-i;
			int index_f=fabs_Index;
			nIVL+=(double)((double)vL*lookup(ti-(double)index_f,sinc_table,ratio));//sinc_GPU(fabs(ti-(double)index_f),1,ratio)*windows_GPU(fabs(ti-(double)index_f),(double)w_hn));
		}
		out_buf[ix*chn_num+1]=(_data_type_)nIVL;
	}//right channel
	//sinc
	//out_sample_count++;//temp counter output
	//break;
	return;
}
