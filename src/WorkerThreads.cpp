#include "WorkerThreads.h"
#include "TraceHBonds.h"

#ifdef PTHREADS

std::list< pthread_t *> threadIdList;
pthread_mutex_t inQueueLock;
pthread_mutex_t runningQueueLock;
pthread_mutex_t outQueueLock;
pthread_mutex_t pauseLock;
pthread_cond_t resumeCond;

//The value is the vector index of worker_data;
std::list< unsigned int > inQueue;
std::list< unsigned int > outQueue;
bool WorkerPauseFlag = false;

std::vector< struct worker_data_s > worker_data;
// Each thread will put their results in one of the elements of this.
std::vector< std::vector<struct HydrogenBond *> *> worker_hb;

void *workerThread(void *threadarg)
{
	// struct thread_detail_s *my_data;
	// my_data = (struct thread_detail_s *) threadarg;
	// unsigned int thread_id = my_data->num;

	while ( 1 )
	{
		CheckWorkerPause();
		pthread_mutex_lock( &inQueueLock );
		unsigned int job;
		bool GotSomeWork=false;
		if (inQueue.size() != 0)
		{
			job = inQueue.front();
			inQueue.pop_front();
			GotSomeWork = true;
		}
		pthread_mutex_unlock( &inQueueLock );

		if ( GotSomeWork )
		{
			switch (worker_data.at(job).jobtype)
			{
				case THREAD_JOB_HBS:
					HBs( worker_data.at(job).hb,
						 worker_data.at(job).cell,
						 worker_data.at(job).hydrogens,
						 worker_data.at(job).acceptors,
						 worker_data.at(job).TrjIdx,
						 worker_data.at(job).rCutoff,
						 worker_data.at(job).angleCutoff,
						 worker_data.at(job).jobnum,
						 worker_data.at(job).num_threads );
					break;
				default:
					break;
			}

			// Put the job in the outQueue
			pthread_mutex_lock( &outQueueLock );
			outQueue.push_back( job );
			pthread_mutex_unlock( &outQueueLock );

			// This should not put anything on the outQueue.
			if ( worker_data.at(job).jobtype == THREAD_JOB_EXIT )
				pthread_exit(NULL);
		}
	}
}

void PauseWorkerThreads()
{
	pthread_mutex_lock(&pauseLock);
	WorkerPauseFlag = true;
	pthread_mutex_unlock(&pauseLock);
}

void ContinueWorkerThreads()
{
	pthread_mutex_lock(&pauseLock);
	WorkerPauseFlag = false;
	pthread_cond_broadcast(&resumeCond);
	pthread_mutex_unlock(&pauseLock);
}

void CheckWorkerPause()
{
	pthread_mutex_lock(&pauseLock);
	while ( WorkerPauseFlag )
		pthread_cond_wait(&resumeCond, &pauseLock);
	pthread_mutex_unlock(&pauseLock);
}

bool StartWorkerThreads(unsigned int NUM_THREADS)
{
	pthread_mutex_init( &inQueueLock     , NULL );
	pthread_mutex_init( &runningQueueLock, NULL );
	pthread_mutex_init( &outQueueLock    , NULL );
	pthread_mutex_init( &pauseLock       , NULL );
	pthread_cond_init( &resumeCond, NULL );

	// start worker threads
	std::list< struct thread_detail_s > thread_detail;

	for (unsigned int j=0; j < NUM_THREADS; ++j)
	{
		pthread_t *tId( new pthread_t );
		threadIdList.push_back(tId);
		struct thread_detail_s Y;
		Y.num = j;
		thread_detail.push_back(Y);

		struct worker_data_s data;
		worker_data.push_back(data);

		int rc;
		rc = pthread_create( tId, NULL, workerThread,
		                     (void *)(&(thread_detail.back()  )) );
		if (rc)
		{
			std::cout << "Error creating threads: " << strerror(rc) << "\n";
			return(false);
		}
	}
	return(true);
}

unsigned int NumWorkerThreads()
{
	return (threadIdList.size());
}

bool StopWorkerThreads()
{
	// Tell the threads to exit by sending THREAD_JOB_EXIT.
	// Make sure the are not paused first.
	ContinueWorkerThreads();

	for (unsigned int j=0; j < threadIdList.size(); ++j)
	{
		pthread_mutex_lock(&inQueueLock);
		worker_data.at(j).jobtype = THREAD_JOB_EXIT;
		inQueue.push_back(j);

		pthread_mutex_unlock(&inQueueLock);
	}
	// Wait queue to empty.
	unsigned int Done_count=0;
	while (Done_count != threadIdList.size() )
	{
		pthread_mutex_lock(&outQueueLock);
		Done_count = outQueue.size();
		pthread_mutex_unlock(&outQueueLock);
	}

	// kill all worker threads
	std::list< pthread_t *>:: iterator i;
	i = threadIdList.begin();
	while ( i!=threadIdList.end() )
	{
		int rc;
		rc = pthread_join( *(*i), NULL);

		if ( rc )
			std::cout << "Error joining thread: " << strerror(rc) << "\n";

		delete (*i);
		i = threadIdList.erase(i);
	}
	pthread_mutex_destroy( &inQueueLock );
	pthread_mutex_destroy( &outQueueLock );

	return(true);
}
#endif
