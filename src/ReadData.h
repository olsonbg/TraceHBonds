#ifndef _ReadData_h
#define _ReadData_h
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <vector>
#include <math.h>
#include "ListOfHBonds.h"
#include "TraceHBonds.h"

int ReadData( const char *filename,
              std::vector<struct HydrogenBond *> *hb,
              std::vector<struct thbAtom *> *atom,
              struct PBC *Cell);

#endif

