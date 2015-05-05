#ifndef _TraceHBonds_h
#define _TraceHBonds_h
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <unistd.h>
#include <math.h>
#include <string.h>
#include <vector>
#include <iostream>
#include <set>
#include "ReadCarMdf.h"
#include "ListOfHBonds.h"
#include "Point.h"

extern bool THB_VERBOSE;

#ifdef DEBUG
#define DEBUG_MSG(str) do { if (DEBUG) std::cout << "DEBUG: " << str << "\n"; } while ( false )
#else
#define DEBUG_MSG(str) do { } while ( false )
#endif

#ifdef __linux
	#define VERBOSE_MSG(str)  do { if (THB_VERBOSE) std::cout << str << "\n"; } while ( false )
	#define VERBOSE_CMSG(str) do { if (THB_VERBOSE) std::cout << str; } while ( false )
	#define VERBOSE_RMSG(str) do { if (THB_VERBOSE) std::cout << str << "\r" << std::flush; } while ( false )
#elif _WIN32
	#define VERBOSE_MSG(str)  do { if (THB_VERBOSE) std::cout << str << "\n" << std::flush; } while ( false )
	#define VERBOSE_CMSG(str) do { if (THB_VERBOSE) std::cout << str << std::flush; } while ( false )
	#define VERBOSE_RMSG(str) do { if (THB_VERBOSE) std::cout << str << "\r" << std::flush; } while ( false )
#endif

//

struct HydrogenBondMatching
{
	std::set<std::string>Hydrogens;
	std::set<std::string>Acceptors;
};

struct HydrogenBondIterator_s
{
	std::vector<struct HydrogenBond *>::iterator begin;
	std::vector<struct HydrogenBond *>::iterator end;
};

void HBs( std::vector<struct HydrogenBond *> *hb,
          Point cell,
          std::vector<struct thbAtom *>*hydrogens,
          std::vector<struct thbAtom *>*acceptors,
          double TrjIdx, double rCutoff, double angleCutoff);

int doArcFile(char *ifilename,
              char *ofPrefix, char *ofSuffix,
              struct HydrogenBondMatching *match,
              double rCutoff, double angleCutoff,
              int NumBins, bool POVRAY);

bool SameAtom( struct thbAtom *A,
               struct thbAtom *B);

bool Trace( ListOfHBonds **HBonds,
            std::vector< struct HydrogenBondIterator_s > *TrjIdx_iter,
            std::vector<struct HydrogenBond *>::iterator iter_hbmain);

std::vector< std::vector<bool> >
Lifetime( std::vector<struct HydrogenBondIterator_s > *TrjIdx_iter,
          std::vector<struct HydrogenBond *>::iterator iter_hbmain,
          unsigned int NumFrames);

void RemoveDuplicates( std::vector<struct HydrogenBond *> *hb,
                       std::vector< struct HydrogenBondIterator_s > *);

void RemoveDuplicatesThread( struct HydrogenBondIterator_s );

template<class T> void DeleteVectorPointers( T v );

std::vector<struct HydrogenBondIterator_s>
TrajectoryIndexIterator( std::vector<struct HydrogenBond *> *hb,
                         unsigned int NumFramesInTrajectory);

#endif
