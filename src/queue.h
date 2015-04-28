#ifndef _queue_h
#define _queue_h

#ifdef PTHREADS
#include <queue>
#include <pthread.h>


class MutexCondition {
	protected:
		pthread_mutex_t m;
		pthread_cond_t  c;
	public:
		MutexCondition();
		~MutexCondition();
		void lock();
		void unlock();
		void wait();
		void broadcast();
};

class lock {
	MutexCondition &mc;
	public:
	lock(MutexCondition &mc);
	~lock();
};

template <typename T>
class Queue {
	private:
		MutexCondition QueueMC;
		std::queue<T> queue;
	public:
		Queue();
		~Queue();
		void push(T work);
		T pop();
};

#endif // PTHREADS
#endif // _queue_h
