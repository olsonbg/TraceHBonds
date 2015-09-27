/**
 * \file
 * \brief Defines inQueue and outQueue
 */
#include "queue.h"
#include "Thread.h"
#include "WorkerThreads.h"

#ifdef PTHREADS

/**
 * \anchor inQueue
 * Queue for starting jobs.
 */
Queue<struct worker_data_s> inQueue;

/**
 * \anchor outQueue
 * Queue of completed jobs.
 */
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
				     wd.angleCutoff );
				break;
			case THREAD_JOB_HBS2:
				HBs( wd.hb,
				     wd.cell,
				     wd.hydrogens,
				     wd.acceptors,
				     wd.coordinates,
				     wd.TrjIdx,
				     wd.rCutoff,
				     wd.angleCutoff );
				break;
			case THREAD_JOB_RMDUPS:
				RemoveDuplicatesThread( *wd.HBit );
				break;
			case THREAD_JOB_TRACE:
				TraceThread( wd.HBStrings, wd.HBit );
				break;
			case THREAD_JOB_CORR:
				CorrelationsThread( wd.vdC,wd.vdI,
				                    wd.vvuiC, wd.vvuiI,
				                    wd.num_threads, wd.jobnum );
				break;
			case THREAD_JOB_CORR_TABLE:
				CorrelationsTableThread( wd.b,
				                         wd.vvuiC, wd.vvuiI,
				                         wd.numHBs, wd.fcutoff,
				                         wd.num_threads, wd.jobnum );
				break;
			case THREAD_JOB_LIFETIME:
				LifetimeThread( wd.b, wd.TrjIdx_iter,
				                wd.num_threads, wd.jobnum);
				break;
			case THREAD_JOB_POSITIONS_CAR:
				PositionsCAR( wd.filename, wd.atom, wd.Cell,
				              wd.hydrogens, wd.acceptors,
				              wd.rCutoff, wd.angleCutoff,
				              wd.saveMemory);
				wd.TrjIdx = wd.Cell->frames;
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
