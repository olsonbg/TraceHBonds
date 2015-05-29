#ifndef _ReadCarMdf_h
#define _ReadCarMdf_h
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <vector>
#include <math.h>
#include "queue.h"
#include "WorkerThreads.h"
#include "cpu.h"
#include "ListOfHBonds.h"
#include "TraceHBonds.h"

bool openfile(const char *filename,
              boost::iostreams::filtering_stream<boost::iostreams::input> *in,
              std::ifstream *ifp);

bool ReadCar(boost::iostreams::filtering_stream<boost::iostreams::input> *in,
             std::vector<struct thbAtom *> *atom,
             struct PBC *Cell, std::vector<Point> *Coordinates);

bool ConnectionsMDF(const char *filename,
                    std::vector<struct thbAtom *> *atom);

bool PositionsCAR(const char *filename,
                  std::vector<struct thbAtom *> *atom,
                  struct PBC *Cell,
                  std::vector<struct thbAtom *> *hydrogens,
                  std::vector<struct thbAtom *> *acceptors,
                  double rCutoff, double angleCutoff);

bool ReadMdf(boost::iostreams::filtering_stream<boost::iostreams::input> *in,
             std::vector<struct thbAtom *> *atom);
#endif

