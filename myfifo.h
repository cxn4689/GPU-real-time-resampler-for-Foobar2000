#include <memory.h>
#include <stdlib.h>

typedef struct {
  char* data;
  size_t allocation;   /* Number of bytes allocated for data. */
  size_t item_size;    /* Size of each item in data */
  size_t begin;        /* Offset of the first byte to read. */
  size_t end;          /* 1 + Offset of the last byte byte to read. */
} fifo_t;

int fifo_create(fifo_t *f,int initial_size,int item_size);
int fifo_read(fifo_t *f, void *des_buf,int size);
int fifo_write(fifo_t *f, void *src_buf,int size);
int get_fifo_size(fifo_t *f);
int fifo_free(fifo_t *f);

int fifo_create(fifo_t *f,int initial_size,int item_size){
	if(initial_size<=0) return 0;
	int alloc_size=item_size*initial_size;
	f->allocation=alloc_size;
	f->item_size=item_size;
	f->begin=0;
	f->end=0;//pos end+1,beg==end means empty
	f->data=(char *)malloc(alloc_size);
	if(f->data==0) return 0;
	return 1;
}

int fifo_read(fifo_t *f, void *des_buf,int size)
{
	if(size>f->end-f->begin) return 0;
	memcpy(des_buf,f->data+f->begin*f->item_size,size*f->item_size);
	if(f->end-f->begin==size){
		f->begin=0;
	}
	return 0;
}
