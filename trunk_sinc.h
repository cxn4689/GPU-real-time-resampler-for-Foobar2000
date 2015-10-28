#define _table_size 65538

#define prec_ 512.0


void generate_table(int sample_half,float factor/*org_sr/des_sr*/,float *sinc_table);//30 points per zero crossing
float lookup(float t,float *sinc_table);//normallized, real t