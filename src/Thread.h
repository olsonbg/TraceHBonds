/**
 * \file
 * \brief Classes for making threads.
 */
#ifndef _Thread_h
#define _Thread_h

#ifdef PTHREADS

#include <pthread.h>

/**
 * Thread class
 */
class Thread
{
	private:
		pthread_t tid;  /**< Thread id                */
		bool running;   /**< Running state of thread  */
		bool detached;  /**< Detached state of thread */
	public:
		/**
		 * Constructor
		 */
		Thread();
		/**
		 * Destructor detached thread with pthread_detach, if that fails, calles pthread_cancel.
		 */
		virtual ~Thread();

		/**
		 * Start a thread with pthread_create()
		 *
		 * \return \c TRUE on pthread_create() success, result of pthread_create() otherwise.
		 */
		int start();
		/**
		 * Call pthread_join() on thread
		 *
		 * \return Result of pthread_join()
		 */
		int join();
		/**
		 * Detatch thread: pthread_detach()
		 *
		 * \return Result of pthread_detach()
		 */
		int detach();
		/**
		 * \return Thread id, pthread_t
		 */
		pthread_t self();

		/**
		 * Varies depending on job, see WorkerThreads.h
		 */
		virtual void *run() = 0;
};

/**
 * Setup function to run on thread creation.
 */
class MyThread : public Thread {
	public:
		/**
		 * Function to run on thread creation. Varies depending on job type,
		 * see \ref JobTypes in WorkerThreads.h
		 */
		void *run();
};

#endif // PTHREADS
#endif // _Thread_h
