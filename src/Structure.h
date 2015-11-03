#ifndef _Structure_h
#define _Structure_h

#include <fstream>
#include <sstream>
#include "ListOfHBonds.h"
#include "MessageDefines.h"


void
saveStructure(unsigned int NumFramesInTrajectory,
              char *ofPrefix, char *ofSuffix,
              struct PBC *Cell,
              std::vector<struct thbMolecule *> *molecules);

#endif // _Structure_h
