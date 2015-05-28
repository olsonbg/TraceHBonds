#ifndef _WorkerThreads_h
#define _WorkerThreads_h

#ifdef PTHREADS
#include <vector>
#include <list>
#include <queue>
#include <pthread.h>
#include "Point.h"
#include "TraceHBonds.h"
#include "ListOfHBonds.h"
#include "Trace.h"
#include "HydrogenBonds.h"
#include "RemoveDuplicates.h"
#include "Histograms.h"


//Shared data.
const unsigned int THREAD_JOB_HBS           =  1; // HBs().
const unsigned int THREAD_JOB_RMDUPS        =  2; // RemoveDuplicatesThread().
const unsigned int THREAD_JOB_TRACE         =  3; // TraceThread().
const unsigned int THREAD_JOB_CORR          =  4; // CorrelationsThread().
const unsigned int THREAD_JOB_LIFETIME      =  5; // LifetimeThread().
const unsigned int THREAD_JOB_POSITIONS_CAR =  6; // PositionsCAR().
const unsigned int THREAD_JOB_PAUSE         = 90; // Pause thread.
const unsigned int THREAD_JOB_EXIT          = 99; // Exit thread.

//The queues.

struct worker_data_s
{
	unsigned int jobtype; // The type of job to run (e.g. THREAD_JOB_HBS)

	// The next two values are Used for loops which are split between more than
	// one thread.  e.g. for(i=jobnum; i < end ; i += num_threads)
	// 
	// Common to all jobs
	// 
	unsigned int jobnum; // 0..num_threads-1
	unsigned int num_threads;
	//
	// HBs
	// 
	std::vector<struct HydrogenBond *> *hb;
	Point cell;
	std::vector<struct thbAtom *>*hydrogens;
	std::vector<struct thbAtom *>*acceptors;
	unsigned int TrjIdx;
	double rCutoff;
	double angleCutoff;
	//
	// Trace and RemoveDuplicates
	// 
	struct HydrogenBondIterator_s *HBit;
	std::vector<ListOfHBonds *> *HBStrings;
	//
	// Lifetime
	// 
	std::vector< std::vector<unsigned int> > *vvuiC;
	std::vector< std::vector<unsigned int> > *vvuiI;
	std::vector< double > *vdC;
	std::vector< double > *vdI;
	std::vector< std::vector<bool> > *b;
	std::vector<struct HydrogenBondIterator_s> *TrjIdx_iter;
	//
	// ReadCarMdf
	//
	const char *filename;
	std::vector<struct thbAtom *> *atom;
	struct PBC *Cell;
};

struct thread_detail_s {
	unsigned int num;
};

#endif // PTHREADS
#endif // _WorkerThreads.h
