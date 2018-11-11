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
enum class Flags : unsigned int {
	ZERO          = 0,
	ONE           = 1 << 0,
	VERBOSE       = 1 << 1, /**< Be verbose                                */
	LIFETIME      = 1 << 2, /**< Calculate lifetime correlations           */
	SIZE_HIST     = 1 << 3, /**< Generate size histograms                  */
	LENGTHS       = 1 << 4, /**< Generate hydrogen bond length list        */
	POVRAY        = 1 << 5, /**< Output in PovRay format                   */
	NEIGHBOR_HIST = 1 << 6, /**< Generge neighbor histograms               */
	ANGLES        = 1 << 7, /**< Generate hydrogen bond angle list         */
	JSON          = 1 << 8, /**< Output in JSON format                     */
	INCELL        = 1 << 9, /**< All chains start inside PBC cell          */
	JSONALL       = 1 <<10, /**< Include all atoms in output, implies JSON, and not INCELL */
	LIST          = 1 <<11, /**< Generate list of hydrogen bonds          */

	/** Do all possible calculations */
	ALL           = (SIZE_HIST|LENGTHS|NEIGHBOR_HIST|LIFETIME|ANGLES|LIST),
	COUNT         = 14
};

unsigned int operator& (unsigned int l, Flags r);
unsigned int operator| (unsigned int l, Flags r);

Flags operator& (Flags l, Flags r);
Flags operator| (Flags l, Flags r);
#endif // _flags_h
