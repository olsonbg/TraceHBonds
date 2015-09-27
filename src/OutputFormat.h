/**
 * \file
 * \brief Set width and precision of numbers for stream output
 *
 **/
#ifndef _OutputFormat_h
#define _OutputFormat_h

#include <iostream>
#include <iomanip>

/**
 * Set width and precision of numbers for stream output
 */
class OFmt
{
	int myWidth;
	int myPrecision;
	public:
		/**
		 * Constructor
		 *
		 * \param[in] width      Set width of number for output stream
		 * \param[in] precision  Set precision of number for output stream.
		 */
		OFmt( int width, int precision )
			: myWidth( width )
			  , myPrecision( precision )
	{
	}

	/**
	 * Send formatted number to stream
	 */
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
