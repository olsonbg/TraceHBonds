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

// Macros

typedef std::vector<unsigned int> vui;
typedef std::vector< vui > vvui;

//

struct HydrogenBondMatching
{
	std::set<std::string>Hydrogens;
	std::set<std::string>Acceptors;
};

void HBs( std::vector<struct HydrogenBond *> *hb,
		  std::vector<double>cell,
		  std::vector<struct thbAtom *>*hydrogens,
		  std::vector<struct thbAtom *>*acceptors,
		  double TrjIdx, double rCutoff, double angleCutoff,
		  unsigned int ThreadID=0, unsigned int Threads=1);

int doArcFile(char *ifilename,
              char *ofPrefix, char *ofSuffix,
              struct HydrogenBondMatching *match,
              double rCutoff, double angleCutoff,
              int NumBins, bool POVRAY);

int makeHistograms( std::ostream *out,
                    std::vector<ListOfHBonds *> HBStrings,
                    std::string CC, unsigned int NumBins,
                    struct PBC *Cell, unsigned int TrjIdx,
                    bool POVRAY);

template<class T> bool alloc_vector(std::vector< std::vector<T> > *v,
                                    T val,
                                    unsigned int nelem,
                                    unsigned int melem);

template<class T> bool alloc_vector( std::vector<T> *v,
                                     T val,
                                     unsigned int nelem);

bool SameAtom( struct thbAtom *A,
               struct thbAtom *B);

bool Trace( ListOfHBonds **HBonds,
            std::vector< std::vector<struct HydrogenBond *>::iterator > *TrjIdx_iter,
            std::vector<struct HydrogenBond *>::iterator iter_hbmain);

void RemoveDuplicates( std::vector<struct HydrogenBond *> *hb,
            std::vector< std::vector<struct HydrogenBond *>::iterator > *);


template<class T> void DeleteVectorPointers( T v );

std::vector< std::vector<struct HydrogenBond *>::iterator >
TrajectoryIndexIterator( std::vector<struct HydrogenBond *> *hb);

#endif
