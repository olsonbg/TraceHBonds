#ifndef _timedoutput_h
#define _timedoutput_h

#include <iostream>
#include <string>
#include <time.h>

class timedOutput
{
	private:
		unsigned int max;
		bool maxSet;
		time_t starttime;
		double timeout;
		std::string prefix;

	public:
		timedOutput( unsigned int m, double diff, std::string s_prefix );
		timedOutput( double diff, std::string s_prefix );
		timedOutput( unsigned int m );
		void maximum( unsigned int m );
		void message( std::string msg );
		void print(unsigned int c);
		void print(unsigned int c, std::string extra);
		void print(unsigned int c, std::string extra, unsigned int c2);
};
#endif // _timedoutput_h
