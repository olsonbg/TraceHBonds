#ifndef _Lifetime_h
#define _Lifetime_h
#include "TraceHBonds.h"

extern bool THB_VERBOSE;

typedef std::vector<struct HydrogenBond *> HBVec;
typedef std::vector<struct HydrogenBondIterator_s> HBVecIter;

std::vector< std::vector<bool> >
Lifetime( std::vector<struct HydrogenBondIterator_s > *TrjIdx_iter );

#endif // _Lifetime_h

