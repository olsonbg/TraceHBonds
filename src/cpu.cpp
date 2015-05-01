#include "cpu.h"

#ifdef _WIN32
#include <windows.h>
#else
#include <unistd.h>
#endif

unsigned int NumberOfCPUs()
{
#ifdef WIN32
	SYSTEM_INFO sysinfo;
	GetSystemInfo(&sysinfo);
	return sysinfo.dwNumberOfProcessors;
#elif __linux
	return sysconf( _SC_NPROCESSORS_ONLN );
#else
#error "OS not supported!"
	return 0;
#endif // WIN32, __linux
}

