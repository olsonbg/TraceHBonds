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
#include "ReadData.h"
#include "ListOfHBonds.h"

int doAllFiles(char *progname,
               char *fPrefix , char *fSuffix, int first, int last,
               char *ofPrefix, char *ofSuffix,
               int NumBins, bool POVRAY );

int doFrame(const char *ifile, const char *ofile,
            unsigned int NumBins, bool POVRAY);

int makeHistograms( std::ostream *out,
                     std::vector<ListOfHBonds *> HBStrings,
                     std::string CC, unsigned int NumBins, 
                     bool POVRAY, struct PBC *Cell);

bool alloc_vector(std::vector<unsigned int> *v,
                  unsigned int val,
                  unsigned int nelem,
                  unsigned int melem);

bool alloc_vector(std::vector< std::vector<unsigned int> > *v,
                  unsigned int val,
                  unsigned int nelem);

bool SameAtom( struct thbAtom *A,
               struct thbAtom *B);

bool Trace( ListOfHBonds **HBonds,
            std::vector<struct HydrogenBond *> hb,
            unsigned int current);

void RemoveDuplicates( std::vector<struct HydrogenBond *> *hb );

// void DeleteVectorPointers( std::vector<struct HydrogenBond *> v);
// void DeleteVectorPointers( std::vector<struct thbAtom *> v);
template<class T> void DeleteVectorPointers( T v );

#endif
