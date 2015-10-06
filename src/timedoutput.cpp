#include "timedoutput.h"

timedOutput::timedOutput( unsigned int m, double diff, std::string s_prefix ) {
	max = m;
	maxSet = true;
	starttime = time(NULL);
	timeout = diff;
	prefix = s_prefix;
}

timedOutput::timedOutput( double diff, std::string s_prefix ) {
	max = 0;
	maxSet = false;
	starttime = time(NULL);
	timeout = diff;
	prefix = s_prefix;
}

void
timedOutput::maximum(unsigned int m) {
	max = m;
	maxSet = true;
}

void
timedOutput::message(std::string msg) {
	prefix = msg;
}

void
timedOutput::print(unsigned int c) {
	if ( (difftime(time(NULL), starttime) > timeout ) || ( c == max ) )
	{
		std::cout << prefix << c;
		if ( maxSet ) {
			std::cout << "/" << max; }
		std::cout << "\r" << std::flush;
		starttime = time(NULL);
	}
}

void
timedOutput::print(unsigned int c, std::string extra) {
	if ( (difftime(time(NULL), starttime) > timeout ) || ( c == max ) )
	{
		std::cout << prefix << c;
		if ( maxSet ) {
			std::cout << "/" << max; }
		std::cout << " (" << extra << ")\r" << std::flush;
		starttime = time(NULL);
	}
}

void
timedOutput::print(unsigned int c, std::string extra, unsigned int c2) {
	if ( (difftime(time(NULL), starttime) > timeout ) || ( c == max ) )
	{
		std::cout << prefix << c;
		if ( maxSet ) {
			std::cout << "/" << max; }
		std::cout << " (" << extra << c2 << ")\r" << std::flush;
		starttime = time(NULL);
	}
}
