/**
 * \file
 * \brief   Flags for command line
 * \author  Brian G. Olson (olsonbg@gmail.com)
 * \date    14 October 2015
 */
#ifndef _flags_h
#define _flags_h

extern bool THB_VERBOSE;

/**
 * Flags for command line options
 */
enum {
	VERBOSE       = 0x01, /**< Be verbose                         */
	LIFETIME      = 0x02, /**< Calculate lifetime correlations    */
	SIZE_HIST     = 0x04, /**< Generate size histograms           */
	LENGTHS       = 0x08, /**< Generate hydrogen bond length list */
	POVRAY        = 0x10, /**< Ouput in PovRay format             */
	NEIGHBOR_HIST = 0x20, /**< Generge neighbor histograms        */
	ANGLES        = 0x40, /**< Generate hydrogen bond angle list  */
	JSON          = 0x80, /**< Ouput in JSON format               */
	/** Do all possible calculations */
	ALL           = (SIZE_HIST|LENGTHS|NEIGHBOR_HIST|LIFETIME|ANGLES),
};

#endif // _flags_h
