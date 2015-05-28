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
#include "MessageDefines.h"
#include "ReadCarMdf.h"
#include "ListOfHBonds.h"
#include "Lifetime.h"
#include "Point.h"

extern bool THB_VERBOSE;

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



template<class T> void DeleteVectorPointers( std::vector<T*> v );

// std::vector<struct HydrogenBondIterator_s>
// TrajectoryIndexIterator( std::vector<struct HydrogenBond *> *hb,
//                          unsigned int NumFramesInTrajectory);

void
TrajectoryIndexIterator( std::vector< struct HydrogenBondIterator_s > *TrjIdx_iter,
                         std::vector< struct HydrogenBond *> *hb);
#endif
