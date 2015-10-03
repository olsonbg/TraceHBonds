/**
 * \file
 * \author Brian G. Olson
 * \date   30 April 2015
 * \brief  Typedefs for commonly used vectors
 *
 **/
#ifndef _VectorTypes_h
#define _VectorTypes_h

#include <vector>

/** Vector of unsigned integers */
typedef std::vector<unsigned int> vui;

/** Vector of vector of unsigned integers */
typedef std::vector< vui > vvui;

/** Vector of doubles */
typedef std::vector< double > vd;

/** Vector of vector of doubles */
typedef std::vector< vd > vvd;

#endif // _VectorTypes_h
