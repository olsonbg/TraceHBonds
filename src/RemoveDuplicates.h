#ifndef _RemoveDuplicates_h
#define _RemoveDuplicates_h

#include "TraceHBonds.h"
#include "queue.h"
#include "WorkerThreads.h"
#include "cpu.h"

extern bool THB_VERBOSE;

#ifdef PTHREADS
extern Queue<struct worker_data_s> inQueue;
extern Queue<struct worker_data_s> outQueue;
#endif

typedef std::vector<struct HydrogenBond *> HBVec;
typedef std::vector<struct HydrogenBondIterator_s> HBVecIter;

struct sort_Neighbors {
	bool operator()(const HBVec::iterator &left, const HBVec::iterator &right)
	{
		return (*left)->length < (*right)->length;
	}
};

bool deleteMarked( struct HydrogenBond *hb );
void removeMarked( HBVec *hb );
void RemoveDuplicates( std::vector<struct HydrogenBond *> *hb,
                       std::vector< struct HydrogenBondIterator_s > *);

void RemoveDuplicatesThread( struct HydrogenBondIterator_s );

bool SameAtom( struct thbAtom *A,
               struct thbAtom *B);

#endif // _RemoveDuplicates_h

