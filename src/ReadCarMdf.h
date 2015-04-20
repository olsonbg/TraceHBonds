#ifndef _ReadCarMdf_h
#define _ReadCarMdf_h
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <vector>
#include <math.h>
#include "ListOfHBonds.h"
#include "TraceHBonds.h"

int ReadCarMdf( const char *filename,
                std::vector<struct thbAtom *> *atom,
                struct PBC *Cell);

int ReadCar( const char *filename,
            std::vector<struct thbAtom *> *atom,
            struct PBC *Cell );

bool ReadMdf( const char *filename,
             std::vector<struct thbAtom *> *atom);
#endif

