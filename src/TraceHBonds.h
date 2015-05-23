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
	std::vector<std::string>Hydrogens;
	std::vector<std::string>Acceptors;
};

struct HydrogenBondIterator_s
{
	std::vector<struct HydrogenBond *>::iterator begin;
	std::vector<struct HydrogenBond *>::iterator end;
};


int doArcFile(char *ifilename,
              char *ofPrefix, char *ofSuffix,
              struct HydrogenBondMatching *match,
              double rCutoff, double angleCutoff,
              int NumBins, bool POVRAY);



std::vector< std::vector<bool> >
Lifetime( std::vector<struct HydrogenBondIterator_s > *TrjIdx_iter );

template<class T> void DeleteVectorPointers( std::vector<T*> v );

// std::vector<struct HydrogenBondIterator_s>
// TrajectoryIndexIterator( std::vector<struct HydrogenBond *> *hb,
//                          unsigned int NumFramesInTrajectory);

void
TrajectoryIndexIterator( std::vector< struct HydrogenBondIterator_s > *TrjIdx_iter,
                         std::vector< struct HydrogenBond *> *hb);
#endif
