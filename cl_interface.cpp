#include "stdafx.h"
#include <CL/cl.h>
#include <stdio.h>
#include <stdlib.h>
#include <tchar.h>
#include <memory.h>
#include<string.h>
#include "cl_interface.h"
#include "trunk_sinc.h"
#include <Windows.h>
#include<string>
//cl
static cl_context __context=0;
static cl_command_queue __cmd_queue=0;
static cl_device_id *__devices=0;
static cl_program __program=0;
static cl_kernel __kernel=0;
static cl_mem __buffer[3];
//cl
//get dll path use
#include "dll_global_var.h"

char *ReadSources(const char *fileName);

cl_platform_id GetIntelOCLPlatform_intel()
{
    cl_platform_id pPlatforms[10] = { 0 };
    char pPlatformName[128] = { 0 };

    cl_uint uiPlatformsCount = 0;
    cl_int err = clGetPlatformIDs(10, pPlatforms, &uiPlatformsCount);
    for (cl_uint ui = 0; ui < uiPlatformsCount; ++ui)
    {
        err = clGetPlatformInfo(pPlatforms[ui], CL_PLATFORM_NAME, 128 * sizeof(char), pPlatformName, NULL);
        if ( err != CL_SUCCESS )
        {
            printf("ERROR: Failed to retreive platform vendor name.\n", ui);
            return NULL;
        }

        if (!strcmp(pPlatformName, "Intel(R) OpenCL"))
            return pPlatforms[ui];
    }

    return NULL;
}//intel

cl_int oclGetPlatformID_nvidia(cl_platform_id* clSelectedPlatformID)
{
    char chBuffer[1024];
    cl_uint num_platforms; 
    cl_platform_id* clPlatformIDs;
    cl_int ciErrNum;
    *clSelectedPlatformID = NULL;

    // Get OpenCL platform count
    ciErrNum = clGetPlatformIDs (0, NULL, &num_platforms);
    if (ciErrNum != CL_SUCCESS)
    {
        return -1000;
    }
    else 
    {
        if(num_platforms == 0)
        {
            return -2000;
        }
        else 
        {
            // if there's a platform or more, make space for ID's
            if ((clPlatformIDs = (cl_platform_id*)malloc(num_platforms * sizeof(cl_platform_id))) == NULL)
            {
                return -3000;
            }

            // get platform info for each platform and trap the NVIDIA platform if found
            ciErrNum = clGetPlatformIDs (num_platforms, clPlatformIDs, NULL);
            for(cl_uint i = 0; i < num_platforms; ++i)
            {
                ciErrNum = clGetPlatformInfo (clPlatformIDs[i], CL_PLATFORM_NAME, 1024, &chBuffer, NULL);
                if(ciErrNum == CL_SUCCESS)
                {
                    if(strstr(chBuffer, "NVIDIA") != NULL)
                    {
                        *clSelectedPlatformID = clPlatformIDs[i];
                        break;
                    }
                }
            }

            // default to zeroeth platform if NVIDIA not found
            if(*clSelectedPlatformID == NULL)
            {
                *clSelectedPlatformID = clPlatformIDs[0];
            }

            free(clPlatformIDs);
        }
    }

    return CL_SUCCESS;
}

int platform_initial()
{
	cl_int err;
	size_t cb=0;
	/*cl_platform_id intel_platform_id = GetIntelOCLPlatform();
    if( intel_platform_id == NULL )
    {
        printf("ERROR: Failed to find Intel OpenCL platform.\n");
        return -1;
    }
    cl_context_properties context_properties[3] = {CL_CONTEXT_PLATFORM, (cl_context_properties)intel_platform_id, NULL };*/
	//nvidia
	cl_platform_id platform_id=0;
	platform_id=GetIntelOCLPlatform_intel();
    if( platform_id == NULL )
    {
        printf("ERROR: Failed to find OpenCL platform.\n");
        return -1;
    }
    cl_context_properties context_properties[3] = {CL_CONTEXT_PLATFORM, (cl_context_properties)platform_id, NULL };
	//nvidia
	__context=clCreateContextFromType(context_properties, CL_DEVICE_TYPE_CPU, NULL, NULL, NULL);
	clGetContextInfo(__context,CL_CONTEXT_DEVICES,0,NULL,&cb);
	__devices=(cl_device_id *)malloc(cb);
	clGetContextInfo(__context,CL_CONTEXT_DEVICES,cb,__devices,NULL);
	__cmd_queue=clCreateCommandQueue(__context, __devices[0], 0, NULL);
	/*get path from dll*/
	char _cl_path[1000];
	HINSTANCE gdll_ptr = get_dll_hdl();
	GetModuleFileName((HMODULE)gdll_ptr,(LPWSTR)_cl_path,sizeof(_cl_path));
	int errCode;
    int size = WideCharToMultiByte(CP_ACP, 0,(LPWSTR)_cl_path, -1, 0, 0, NULL, NULL); 
	char* buf__ = new char[size]; 
	WideCharToMultiByte(CP_ACP, 0, (LPWSTR)_cl_path, -1, buf__, size, NULL, NULL); //unicode to multi-bytes "windows.h"
	std::string _cl_cs(buf__);
	delete []buf__;
	int pos_f = _cl_cs.find("foo_dsp_resampler.dll");
	std::string _full_cl_path = _cl_cs.substr(0,pos_f) + "kernel.cl";
	//
	char *sources=ReadSources(_full_cl_path.c_str());
	//char *sources=ReadSources("D:\\Software_Installed\\Foobar2000\\user-components\\kernel.cl");//replaced by upon codes
	__program=clCreateProgramWithSource(__context, 1, (const char**)&sources, NULL, NULL);
	err=clBuildProgram(__program, 0, NULL, NULL, NULL, NULL);
	{
		cl_build_status build_status;
		clGetProgramBuildInfo(__program, __devices[0], CL_PROGRAM_BUILD_STATUS, sizeof(cl_build_status), &build_status, NULL);
		char *build_log;                                                        
		size_t ret_val_size;
		clGetProgramBuildInfo(__program, __devices[0], CL_PROGRAM_BUILD_LOG, 0, NULL, &ret_val_size);
		build_log = new char[ret_val_size+1];                                                        
		clGetProgramBuildInfo(__program, __devices[0], CL_PROGRAM_BUILD_LOG, ret_val_size, build_log, NULL);   
		build_log[ret_val_size] = '\0';
		printf("%s\n",build_log);
		delete []build_log;
	}
	__kernel = clCreateKernel(__program, "interpolation", NULL);
	free(sources);
	//memory object
	__buffer[0]=__buffer[1]=__buffer[2]=0;
	/*buffer[0]=clCreateBuffer(context,CL_MEM_READ_ONLY,sizeof(cl_int)*512*2,NULL,&err);
	buffer[1]=clCreateBuffer(context,CL_MEM_WRITE_ONLY,sizeof(cl_int)*1024*2,NULL,&err);*/
	//end of memory
	return err;
}

void close_platform()
{
    if( __kernel ) {clReleaseKernel( __kernel ); __kernel = NULL;}
    if( __program ) {clReleaseProgram( __program );__program = NULL;}
    if( __cmd_queue ) {clReleaseCommandQueue( __cmd_queue ); __cmd_queue = NULL;}
    if( __context ) {clReleaseContext( __context );__context = NULL;}
	/*if(buffer[0]) clReleaseMemObject(buffer[0]);
	if(buffer[1]) clReleaseMemObject(buffer[1]);*/
	__buffer[0]=__buffer[1]=__buffer[2]=0;
}

char *ReadSources(const char *fileName);

int launch_cl_kernel(_data_type_ *buffer_process,int buf_size,_data_type_ *out_buf,
	int out_size,int beg_index,int out_samples,int chn_num,int ori_sr,int des_sr,int w_hn,float *sinc_table)
{
	//
	size_t global[]={0};
	global[0]=out_size;
	/*if(out_size>256){
		global[0]=out_size/256+1;
	}
	else global[0]=1;
	size_t local_size[]={0};
	if(out_size>256){
		local_size[0]=256;
	}
	else local_size[0]=out_size;*/
	//
	cl_int err;
	//initial ok
	/*_kernel void interpolation(__global _data_type_ *buffer_process,int buf_size,__global _data_type_ *out_buf,int out_size,
	int beg_index,int out_samples,int chn_num,int ori_sr,int des_sr,int w_hn)*/
	__buffer[0]=clCreateBuffer(__context,CL_MEM_READ_ONLY,sizeof(_data_type_)*buf_size*chn_num,NULL,&err);
	__buffer[1]=clCreateBuffer(__context,CL_MEM_WRITE_ONLY,sizeof(_data_type_)*out_size*chn_num,NULL,&err);
	__buffer[2]=clCreateBuffer(__context,CL_MEM_WRITE_ONLY,sizeof(_data_type_)*_table_size,NULL,&err);//fix 2000
	//initial buffer
	err=clEnqueueWriteBuffer(__cmd_queue,__buffer[0],CL_TRUE,0,sizeof(_data_type_)*buf_size*chn_num,buffer_process,0,NULL,NULL);
	err=clEnqueueWriteBuffer(__cmd_queue,__buffer[2],CL_TRUE,0,sizeof(_data_type_)*_table_size,sinc_table,0,NULL,NULL);//fix 2000
	//set parameter
	err=clSetKernelArg(__kernel,
					   0,
					   sizeof(cl_mem),
					   (void *)&(__buffer[0]));
	err=clSetKernelArg(__kernel,
					   1,
					   sizeof(cl_int),
					   (void *)&buf_size);
	err=clSetKernelArg(__kernel,
					   2,
					   sizeof(cl_mem),
					   (void *)&(__buffer[1]));
	err=clSetKernelArg(__kernel,
					   3,
					   sizeof(cl_int),
					   (void *)&out_size);
	err=clSetKernelArg(__kernel,
					   4,
					   sizeof(cl_int),
					   (void *)&beg_index);
	err=clSetKernelArg(__kernel,
					   5,
					   sizeof(cl_int),
					   (void *)&out_samples);
	err=clSetKernelArg(__kernel,
					   6,
					   sizeof(cl_int),
					   (void *)&chn_num);
	err=clSetKernelArg(__kernel,
					   7,
					   sizeof(cl_int),
					   (void *)&ori_sr);
	err=clSetKernelArg(__kernel,
					   8,
					   sizeof(cl_int),
					   (void *)&des_sr);
	err=clSetKernelArg(__kernel,
					   9,
					   sizeof(cl_int),
					   (void *)&w_hn);
	err=clSetKernelArg(__kernel,
					   10,
					   sizeof(cl_mem),
					   (void *)&(__buffer[2]));
	//launch kernel
	cl_event event_0;
	err=clEnqueueNDRangeKernel(__cmd_queue,__kernel,1,NULL,global,0,0,NULL,&event_0);
	if(err!=CL_SUCCESS){printf("executing kernel error!\n");return 0;}
	err=clWaitForEvents(1,&event_0);
	if(err!=CL_SUCCESS){printf("executing kernel event waiting error!\n");return 0;}
	//release event
	clReleaseEvent(event_0);
	//read memory
	cl_event event_1;
	err=clEnqueueReadBuffer(__cmd_queue,__buffer[1],CL_TRUE,0,sizeof(_data_type_)*out_size*chn_num,out_buf,0,NULL,&event_1);
	err=clWaitForEvents(1,&event_1);

	//release event
	clReleaseEvent(event_1);

	clReleaseMemObject(__buffer[0]);
	clReleaseMemObject(__buffer[1]);
	clReleaseMemObject(__buffer[2]);

	return 1;
}

char *ReadSources(const char *fileName)
{
    FILE *file = fopen(fileName, "rb");
    if (!file)
    {
        printf("ERROR: Failed to open file '%s'\n", fileName);
        return NULL;
    }

    if (fseek(file, 0, SEEK_END))
    {
        printf("ERROR: Failed to seek file '%s'\n", fileName);
		fclose(file);
        return NULL;
    }

    long size = ftell(file);
    if (size == 0)
    {
        printf("ERROR: Failed to check position on file '%s'\n", fileName);
		fclose(file);
        return NULL;
    }

    rewind(file);

    char *src = (char *)malloc(sizeof(char) * size + 1);
    if (!src)
    {
        printf("ERROR: Failed to allocate memory for file '%s'\n", fileName);
		fclose(file);
        return NULL;
    }

    printf("Reading file '%s' (size %ld bytes)\n", fileName, size);
    size_t res = fread(src, 1, sizeof(char) * size, file);
    if (res != sizeof(char) * size)
    {
        printf("ERROR: Failed to read file '%s'\n", fileName);
		fclose(file);
		free(src);
        return NULL;
    }

    src[size] = '\0'; /* NULL terminated */
    fclose(file);

    return src;
}