#ifndef _data_type_
#define _data_type_ float
#endif

int platform_initial();
void close_platform();
int launch_cl_kernel(_data_type_ *buffer_process,int buf_size,_data_type_ *out_buf,
	int out_size,int beg_index,int out_samples,int chn_num,int ori_sr,int des_sr,int w_hn,float *sinc_table);