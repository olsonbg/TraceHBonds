#ifndef _OutputFormat_h
#define _OutputFormat_h

#include <iostream>
#include <iomanip>

class OFmt
{
	int myWidth;
	int myPrecision;
	public:
		OFmt( int width, int precision )
			: myWidth( width )
			  , myPrecision( precision )
	{
	}

	friend std::ostream&
		operator<<( std::ostream& dest, OFmt const& fmt )
		{
			dest.setf( std::ios_base::fixed, std::ios_base::floatfield );
			dest.precision( fmt.myPrecision );
			dest.width( fmt.myWidth );
			return dest;
		}
};
#endif
