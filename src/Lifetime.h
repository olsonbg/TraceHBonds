#ifndef _Lifetime_h
#define _Lifetime_h
#include "TraceHBonds.h"
#include "queue.h"
#include "WorkerThreads.h"
#include "cpu.h"

extern bool THB_VERBOSE;

typedef std::vector<struct HydrogenBond *> HBVec;
typedef std::vector<struct HydrogenBondIterator_s> HBVecIter;

void
Lifetime( std::vector< std::vector<bool> >*b,
          std::vector<struct HydrogenBondIterator_s > *TrjIdx_iter );

void
LifetimeThread(std::vector< std::vector<bool> >*b,  HBVecIter *TrjIdx_iter,
               unsigned int NumThreads=1, unsigned int ThreadID=0);
#endif // _Lifetime_h

