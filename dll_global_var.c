#include<Windows.h>

static HINSTANCE _global_hinstDll = 0;//for get path use

void set_dll_hdl(HINSTANCE ptr)
{
	_global_hinstDll = ptr;
}
HINSTANCE get_dll_hdl()
{
	return _global_hinstDll;
}