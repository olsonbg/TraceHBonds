#ifndef _WorkerThreads_h
#define _WorkerThreads_h

#ifdef PTHREADS
#include <vector>
#include <list>
#include <pthread.h>
#include <sys/signal.h>


//Shared data.
#define THREAD_JOB_HBS     1 // Run HBs().
#define THREAD_JOB_PAUSE  90 // Pause thread.
#define THREAD_JOB_EXIT   99 // Exit thread.

//The queues.

struct worker_data_s
{
	unsigned int jobtype; // The type of job to run (e.g. THREAD_JOB_HBS)

	// The next two valuse are Used for loops which are split between more than
	// one thread.  e.g. for(i=jobnum; i < end < i += num_threads)
	unsigned int jobnum; // 0..num_threads-1
	unsigned int num_threads;
	std::vector<struct HydrogenBond *> *hb;
	std::vector<double>cell;
	std::vector<struct thbAtom *>*hydrogens;
	std::vector<struct thbAtom *>*acceptors;
	unsigned int TrjIdx;
	double rCutoff;
	double angleCutoff;
};

struct thread_detail_s {
	unsigned int num;
};

void *workerThread(void *);
unsigned int NumWorkerThreads();
void CheckWorkerPause(void);
void PauseWorkerThreads(void);
void ContinueWorkerThreads(void);
bool StartWorkerThreads(unsigned int NUM_THREADS);
bool StopWorkerThreads(void);

#endif
#endif
