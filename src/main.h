#include <unistd.h>
#ifdef _WIN32
#include <windows.h>
#endif


#ifdef __linux
	const unsigned int numCPUs = sysconf( _SC_NPROCESSORS_ONLN );
#elif _WIN32
	SYSTEM_INFO sysinfo;
	GetSystemInfo(&sysinfo);
	const unsigned int numCPUs = sysinfo.dwNumberOfProcessors;
#else
#error "OS not supported!"
#endif // __linux

const unsigned int NumThreads = numCPUs;
