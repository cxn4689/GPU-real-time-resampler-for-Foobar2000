#include "stdafx.h"
#include <stdlib.h>
#include <memory.h>
#include <math.h>
#include "fifo.h"
#include "SRC_Stream.h"
#include <stdio.h>
#include <omp.h>
#include "trunk_sinc.h"

#define _CUDA_

/*
float double 混用造成精度损失！！！！！！
*/

static _data_type_ *buffer_process=0;//for conversion use size fixed to be 512*2
static int buf_frame_in_use=0;//how many buf(frame) is in used,actually size is doubled

static int ori_sr=44100,des_sr=48000;

static fifo_t fifo_in;
static int fifo_init_flag=0;

static int in_samples=0;
static int out_samples=0;//input&output frame counter

static double ratio=44100.0/48000.0;//fix ratio

static int beg_index=-64;

#define buf_fix_size 2048 //frame, actual size= multiple of channel number
#define w_hn 64   //half window size
#define channel_num 2 //channel number
#define window_size 128 //full window size
#define pi 3.1415926 //Pi
#define chn_num 2 //channel number

static double ori_sample_rate=44100.0;
static double des_sample_rate=48000.0;

//trunk sinc
static float *sinc_table = 0;//sinc table

int __my_set_parameter(int in_org_sr,int in_des_sr)
{
	in_samples=0;
	out_samples=0;
	beg_index=-w_hn;
	buf_frame_in_use=w_hn;
	ori_sr=in_org_sr;
	des_sr=in_des_sr;
	ori_sample_rate=(double)ori_sr;
	des_sample_rate=(double)des_sr;
	ratio=ori_sample_rate/des_sample_rate;
	//generate sinc table here
	generate_table(w_hn,ratio,sinc_table);
	return 0;
}


//sinc math constant
double windows(double t,double nT){//bn sequence length*sample period, blackman window
	if(t<=nT&&t>=-nT) return 0.42+0.5*cos(pi*t/nT)+0.08*cos(2*pi*t/nT);
	else return 0;
}
double sinc(double t,double fs){
	if(t==0) return 1;
	double x=sin(fs*t*pi);
	double y=(pi*fs*t);
	return x/y;
}
//end of sinc

//CUDA
extern int launch_cl_kernel(_data_type_ *buffer_process,int buf_size,_data_type_ *out_buf,int out_size,int beg_index,
	int out_samples,int chn_num_x,int ori_sr,int des_sr,int w_hn_x,float *sinc_table);
//end of CUDA

int __my_initial_buffer()
{
	if(buffer_process!=0) free(buffer_process);//fix buffer leaking problem
	buffer_process=(_data_type_ *)malloc(sizeof(_data_type_)*buf_fix_size*chn_num);
	if(fifo_init_flag==0){
		fifo_create(&fifo_in,sizeof(_data_type_));
		fifo_init_flag=1;
	}
	in_samples=0;
	out_samples=0;
	beg_index=-w_hn;
	buf_frame_in_use=w_hn;
	memset(buffer_process,0,buf_fix_size*chn_num*sizeof(_data_type_));
	//
	if(sinc_table == 0) sinc_table = new float[_table_size];//fix currently
	memset(sinc_table,0,_table_size*sizeof(float));
	return (int)buffer_process;
}

int __my_get_data_and_src(_data_type_ *in_buf,int in_size,_data_type_ *out_buf,int out_size){//input intem, out_size=out frame,in_size= channel_number*frame= number of items
	//initial buffer
	if(in_size!=0){//fixed memory leaking
		in_samples+=in_size;
		fifo_write(&fifo_in,in_size,in_buf);
	}
	if(buf_frame_in_use<buf_fix_size){
		/*if(out_samples==0){
			buf_frame_in_use=32;
		}*/
		//get fifo size
		int fifo_size=fifo_occupancy(&fifo_in);
		if(fifo_size>(buf_fix_size-buf_frame_in_use)*chn_num)
		{
			fifo_read(&fifo_in,(buf_fix_size-buf_frame_in_use)*chn_num,(void *)((int)buffer_process+buf_frame_in_use*chn_num*sizeof(_data_type_)));
			buf_frame_in_use=buf_fix_size;
		}
		else{
			fifo_read(&fifo_in,fifo_size,(void *)((int)buffer_process+buf_frame_in_use*chn_num*sizeof(_data_type_)));
			buf_frame_in_use+=fifo_size/chn_num;
		}
		if(buf_frame_in_use<buf_fix_size){
			return 0;
		}
	}
	//caculating output sample's frame number
	double A_=((double)out_samples*ratio)-(double)beg_index;
	int out_pred_num=(int)(((double)(buf_fix_size-w_hn+1)-A_)/ratio)+1;
	//begin caculation
	int out_sample_count=0;

#ifdef _CUDA_
	launch_cl_kernel(buffer_process,buf_fix_size,out_buf,out_pred_num,beg_index,out_samples,chn_num,ori_sr,des_sr,w_hn,sinc_table);
#endif
#ifndef _CUDA_
#pragma omp parallel for
	for(int ix=0;ix<out_pred_num;ix++)
	{
		double rs_p_begin=(double)(out_samples+ix)*ratio;//resample point
		double ti=(double)(out_samples+ix)/(double)des_sample_rate;
		int nL=(int)rs_p_begin;
		int nR=nL+1;//caculating left and right
		//cacaulating actuall index
		int a_nL=nL-beg_index;
		int a_nR=nR-beg_index;
		if(a_nL>buf_fix_size-w_hn){break;}
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
				int abs_Index=nR+i;
				int index_f=abs_Index;
				nIVL+=(double)vL*sinc(abs(ti-(double)index_f/(double)ori_sample_rate),ori_sample_rate)*windows(abs(ti-(double)index_f/(double)ori_sample_rate),(double)w_hn/(double)ori_sample_rate);
			}
			for(int i=0;i<w_hn;i++){//right wing
				//window size 32
				int nId=a_nL-i;
				_data_type_ vL=(nId>=0&&nId<=buf_fix_size-1)?buffer_process[nId*chn_num]:0;
				int abs_Index=nL-i;
				int index_f=abs_Index;
				nIVL+=(double)vL*sinc(abs(ti-(double)index_f/(double)ori_sample_rate),ori_sample_rate)*windows(abs(ti-(double)index_f/(double)ori_sample_rate),(double)w_hn/(double)ori_sample_rate);
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
				int abs_Index=nR+i;
				int index_f=abs_Index;
				nIVL+=(double)vL*sinc(abs(ti-(double)index_f/(double)ori_sample_rate),ori_sample_rate)*windows(abs(ti-(double)index_f/(double)ori_sample_rate),(double)w_hn/(double)ori_sample_rate);
			}
			for(int i=0;i<w_hn;i++){//right wing
				//window size 32
				int nId=a_nL-i;
				_data_type_ vL=(nId>=0&&nId<=buf_fix_size-1)?buffer_process[nId*chn_num+1]:0;
				int abs_Index=nL-i;
				int index_f=abs_Index;
				nIVL+=(double)vL*sinc(abs(ti-(double)index_f/(double)ori_sample_rate),ori_sample_rate)*windows(abs(ti-(double)index_f/(double)ori_sample_rate),(double)w_hn/(double)ori_sample_rate);
			}
			out_buf[ix*chn_num+1]=(_data_type_)nIVL;
		}//right channel
		//sinc
		//out_sample_count++;//temp counter output
		//break;
	}//end of while
#endif
	out_sample_count=out_pred_num;
	//adjust out_sample
	out_samples+=out_sample_count;
	//adjust buffer
	int reserve_index=buf_fix_size-window_size;
	memcpy(buffer_process,(char *)buffer_process+reserve_index*chn_num*sizeof(_data_type_),window_size*chn_num*sizeof(_data_type_));
	buf_frame_in_use=window_size;
	//adjust begin index
	beg_index+=buf_fix_size-window_size;
	//check prediction output_num
	if(out_sample_count!=out_pred_num&&out_sample_count!=0){
		printf("prediction fail!\n");
		exit(0);
	}
	return out_sample_count;
}

int __my_flush_and_reset(_data_type_ *out_buf,_data_type_ out_size){
	__my_reset();
	return 0;
}

int __my_clear_all()//on dsp destruction, do this once
{
	if(buffer_process!=0) free(buffer_process);
	if(fifo_init_flag) fifo_delete(&fifo_in);
	if(sinc_table) free(sinc_table);
	fifo_init_flag=0;
	buffer_process=0;
	sinc_table = 0;
	in_samples=0;
	out_samples=0;
	beg_index=-w_hn;
	buf_frame_in_use=w_hn;
	return 0;
}

int __my_reset(){
	if(buffer_process==0) buffer_process=(_data_type_ *)malloc(sizeof(_data_type_)*buf_fix_size*chn_num);
	if(fifo_init_flag==0){
		fifo_create(&fifo_in,sizeof(_data_type_));
		fifo_init_flag=1;
	}
	else{
		fifo_delete(&fifo_in);
		fifo_create(&fifo_in,sizeof(_data_type_));
		fifo_init_flag=1;
	}
	fifo_clear(&fifo_in);
	in_samples=0;
	out_samples=0;
	beg_index=-w_hn;
	buf_frame_in_use=w_hn;
	memset(buffer_process,0,buf_fix_size*chn_num*sizeof(_data_type_));
	//
	if(sinc_table == 0) sinc_table = new float[_table_size];//fix currently
	memset(sinc_table,0,_table_size*sizeof(float));
	//generate sinc table here
	generate_table(w_hn,ratio,sinc_table);
	return 0;
}


int __my_get_data_and_src__backup(_data_type_ *in_buf,int in_size,_data_type_ *out_buf,int out_size){//input intem, out_size=out frame,in_size= channel_number*frame= number of items
	//initial buffer
	in_samples+=in_size;
	if(in_size!=0) fifo_write(&fifo_in,in_size,in_buf);
	if(buf_frame_in_use<buf_fix_size){
		/*if(out_samples==0){
			buf_frame_in_use=32;
		}*/
		//get fifo size
		int fifo_size=fifo_occupancy(&fifo_in);
		if(fifo_size>(buf_fix_size-buf_frame_in_use)*chn_num)
		{
			fifo_read(&fifo_in,(buf_fix_size-buf_frame_in_use)*chn_num,(void *)((int)buffer_process+buf_frame_in_use*chn_num*sizeof(_data_type_)));
			buf_frame_in_use=buf_fix_size;
		}
		else{
			fifo_read(&fifo_in,fifo_size,(void *)((int)buffer_process+buf_frame_in_use*chn_num*sizeof(_data_type_)));
			buf_frame_in_use+=fifo_size/chn_num;
		}
		if(buf_frame_in_use<buf_fix_size){
			return 0;
		}
	}
	//begin caculation
	int out_sample_count=0;
	while(1){
		double rs_p_begin=(double)(out_samples+out_sample_count)*ratio;//resample point
		double ti=(double)(out_samples+out_sample_count)/(double)des_sample_rate;
		int nL=(int)rs_p_begin;//floor function all the same in this project
		int nR=nL+1;//caculating left and right
		//cacaulating actuall index
		int a_nL=nL-beg_index;
		int a_nR=nR-beg_index;
		if(a_nL>buf_fix_size-w_hn){break;}
		int vL=(a_nL>=0&&a_nL<=buf_fix_size-1)?buffer_process[a_nL*chn_num]:0;
		int vR=(a_nR<=buf_fix_size-1&&a_nR>=0)?buffer_process[a_nR*chn_num]:0;
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
				int vL=(nId>=0&&nId<=buf_fix_size-1)?buffer_process[nId*chn_num]:0;
				int abs_Index=nR+i;
				int index_f=abs_Index;
				nIVL+=(double)vL*sinc(abs(ti-(double)index_f/(double)ori_sample_rate),ori_sample_rate)*windows(abs(ti-(double)index_f/(double)ori_sample_rate),(double)w_hn/(double)ori_sample_rate);
			}
			for(int i=0;i<w_hn;i++){//right wing
				//window size 32
				int nId=a_nL-i;
				int vL=(nId>=0&&nId<=buf_fix_size-1)?buffer_process[nId*chn_num]:0;
				int abs_Index=nL-i;
				int index_f=abs_Index;
				nIVL+=(double)vL*sinc(abs(ti-(double)index_f/(double)ori_sample_rate),ori_sample_rate)*windows(abs(ti-(double)index_f/(double)ori_sample_rate),(double)w_hn/(double)ori_sample_rate);
			}
			out_buf[out_sample_count*chn_num]=(_data_type_)nIVL;
		}//left channel
		{
			double nIVL=0;
			//double ti=rs_p_begin;//insert point(Fs')
			for(int i=0;i<w_hn;i++){//right wing
				//window size 32
				int nId=a_nR+i;
				int vL=(nId>=0&&nId<=buf_fix_size-1)?buffer_process[nId*chn_num+1]:0;
				int abs_Index=nR+i;
				int index_f=abs_Index;
				nIVL+=(double)vL*sinc(abs(ti-(double)index_f/(double)ori_sample_rate),ori_sample_rate)*windows(abs(ti-(double)index_f/(double)ori_sample_rate),(double)w_hn/(double)ori_sample_rate);
			}
			for(int i=0;i<w_hn;i++){//right wing
				//window size 32
				int nId=a_nL-i;
				int vL=(nId>=0&&nId<=buf_fix_size-1)?buffer_process[nId*chn_num+1]:0;
				int abs_Index=nL-i;
				int index_f=abs_Index;
				nIVL+=(double)vL*sinc(abs(ti-(double)index_f/(double)ori_sample_rate),ori_sample_rate)*windows(abs(ti-(double)index_f/(double)ori_sample_rate),(double)w_hn/(double)ori_sample_rate);
			}
			out_buf[out_sample_count*chn_num+1]=(_data_type_)nIVL;
		}//right channel
		//sinc
		out_sample_count++;//temp counter output
		//break;
	}//end of while
	//adjust out_sample
	out_samples+=out_sample_count;
	//adjust buffer
	int reserve_index=buf_fix_size-window_size;
	memcpy(buffer_process,(char *)buffer_process+reserve_index*chn_num*sizeof(int),window_size*chn_num*sizeof(int));
	buf_frame_in_use=window_size;
	//adjust begin index
	beg_index+=buf_fix_size-window_size;
	return out_sample_count;
}