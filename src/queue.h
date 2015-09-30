/**
 * \file
 * \author Brian G. Olson
 * \date   27 April 2015
 * \brief  Setup queue of jobs for threads to work on
 *
 */
#ifndef _queue_h
#define _queue_h

#ifdef PTHREADS
#include <queue>
#include <pthread.h>


/**
 * \class MutexCondition
 * \brief MutexCondition
 */
class MutexCondition {
	protected:
		pthread_mutex_t m; /**< Mutex for locking           */
		pthread_cond_t  c; /**< Cond for wait and broadcast */

	public:
		/**
		 * Constructor
		 */
		MutexCondition();
		/**
		 * destructor
		 */
		~MutexCondition();
		/**
		 * lock
		 */
		void lock();
		/**
		 * unlock
		 */
		void unlock();
		/**
		 * wait
		 */
		void wait();
		/**
		 * broadcast
		 */
		void broadcast();
};

/**
 * Mutex for locking threads
 */
class lock {
	MutexCondition &mc;
	public:
	/**
	 * Constructor, sets mutex lock
	 */
	lock(MutexCondition &mc);
	/**
	 * Destructor, removes mutex lock (unlocks)
	 */
	~lock();
};

/**
 * Simple queuing system for threads
 */
template <typename T>
class Queue {
	private:
		MutexCondition QueueMC;
		std::queue<T> queue;
	public:
		Queue();
		~Queue();
		/**
		 * Put a job in the queue. Used for both \ref inQueue and \ref outQueue.
		 */
		void push(T work);

		/**
		 * Take a job off the queue. Used for both \ref inQueue and \ref outQueue.
		 */
		T pop();
};

#endif // PTHREADS
#endif // _queue_h
