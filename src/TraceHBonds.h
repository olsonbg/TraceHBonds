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
#include "ReadCarMdf.h"
#include "ListOfHBonds.h"

// Macros

typedef std::vector<unsigned int> vui;
typedef std::vector< vui > vvui;

//

int doAllFiles(char *progname,
               char *fPrefix , char *fSuffix, int first, int last,
               char *ofPrefix, char *ofSuffix,
               int NumBins, bool POVRAY );

int doFrame(const char *ifile, const char *ofile,
            unsigned int NumBins, bool POVRAY);

int makeHistograms( std::ostream *out,
                     std::vector<ListOfHBonds *> HBStrings,
                     std::string CC, unsigned int NumBins,
                     struct PBC *Cell, unsigned int TrjIdx,
                     bool POVRAY);

bool alloc_vector(std::vector< std::vector<unsigned int> > *v,
                  unsigned int val,
                  unsigned int nelem,
                  unsigned int melem);

template<class T> bool alloc_vector( std::vector<T> *v,
                                     T val,
                                     unsigned int nelem);

bool alloc_vector(struct thbAtom *v,
                  double val,
                  unsigned int nelem);

bool alloc_vector(struct PBC *v,
                  double val,
                  unsigned int nelem);

bool SameAtom( struct thbAtom *A,
               struct thbAtom *B);

bool Trace( ListOfHBonds **HBonds,
            std::vector< std::vector<struct HydrogenBond *>::iterator > *TrjIdx_iter,
            std::vector<struct HydrogenBond *>::iterator iter_hbmain);

void RemoveDuplicates( std::vector<struct HydrogenBond *> *hb,
            std::vector< std::vector<struct HydrogenBond *>::iterator > *);


template<class T> void DeleteVectorPointers( T v );

#endif
