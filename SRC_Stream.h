#include "fifo.h"
#define _data_type_ float


int __my_initial_buffer(); //on dsp creation, do this once

int __my_get_data_and_src(_data_type_ *in_buf,int in_size,_data_type_ *out_buf,int out_size);

int __my_flush_and_reset(_data_type_ *out_buf,int out_size); //on end of file or stop, flush the buffer

int __my_reset();//on end of file or stop, then reset

int __my_clear_all();//on dsp destruction, do this once

int __my_set_parameter(int org_sr,int des_sr);