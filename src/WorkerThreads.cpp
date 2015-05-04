#include "WorkerThreads.h"
#include "queue.h"
#include "Thread.h"
#include "TraceHBonds.h"

#ifdef PTHREADS

Queue<struct worker_data_s> inQueue;
Queue<struct worker_data_s> outQueue;

void *MyThread::run() 
{
	struct worker_data_s wd;
	while ( 1 )
	{
		wd = inQueue.pop();

		switch (wd.jobtype)
		{
			case THREAD_JOB_HBS:
				HBs( wd.hb,
				     wd.cell,
				     wd.hydrogens,
				     wd.acceptors,
				     wd.TrjIdx,
				     wd.rCutoff,
				     wd.angleCutoff,
				     wd.jobnum,
				     wd.num_threads );
				break;
			case THREAD_JOB_EXIT:
				return NULL;
				break;
			default:
				break;
		}

		// Put the job in the outQueue
		outQueue.push(wd);
	}
	return NULL;
}

#endif // PTHREADS