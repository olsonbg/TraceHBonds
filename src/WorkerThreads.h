#ifndef _WorkerThreads_h
#define _WorkerThreads_h

#ifdef PTHREADS
#include <vector>
#include <list>
#include <queue>
#include <pthread.h>
#include "Point.h"


//Shared data.
const unsigned int THREAD_JOB_HBS   =  1; // Run HBs().
const unsigned int THREAD_JOB_PAUSE = 90; // Pause thread.
const unsigned int THREAD_JOB_EXIT  = 99; // Exit thread.

//The queues.

struct worker_data_s
{
	unsigned int jobtype; // The type of job to run (e.g. THREAD_JOB_HBS)

	// The next two valuse are Used for loops which are split between more than
	// one thread.  e.g. for(i=jobnum; i < end < i += num_threads)
	unsigned int jobnum; // 0..num_threads-1
	unsigned int num_threads;
	std::vector<struct HydrogenBond *> *hb;
	Point cell;
	std::vector<struct thbAtom *>*hydrogens;
	std::vector<struct thbAtom *>*acceptors;
	unsigned int TrjIdx;
	double rCutoff;
	double angleCutoff;
};

struct thread_detail_s {
	unsigned int num;
};

#endif // PTHREADS
#endif // _WorkerThreads.h
