/**
 * \file
 * \date  14 June 2018
 * \brief Open a possible compressed file for reading
 */
#ifndef _OpenFile_h
#define _OpenFile_h
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#ifdef USE_LZMA
#include "lzma.h"
#endif

/**
 * Opens a file for reading
 *
 * Determines the type of file and sets up the appropriate
 * boost::iostreams::filtering_stream.
 *
 * \param[in]      filename Name of file
 * \param[in,out]  in       boost filtering_stream
 * \param[in,out]  ifp      ifstream
 *
 * \return \c TRUE if the file could be open, \c FALSE otherwise
 */
bool openfile(const char *filename,
              boost::iostreams::filtering_stream<boost::iostreams::input> *in,
              std::ifstream *ifp);
#endif
