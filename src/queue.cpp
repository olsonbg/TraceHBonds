#include "queue.h"
#include "WorkerThreads.h"
#include<iostream>

#ifdef PTHREADS

MutexCondition::MutexCondition()
{
	pthread_mutex_init(&m, NULL);
	pthread_cond_init (&c, NULL);
}
MutexCondition::~MutexCondition()
{
	pthread_cond_destroy (&c);
	pthread_mutex_destroy(&m);
}
void MutexCondition::lock()      { pthread_mutex_lock    (&m);     }
void MutexCondition::unlock()    { pthread_mutex_unlock  (&m);     }
void MutexCondition::wait()      { pthread_cond_wait     (&c, &m); }
void MutexCondition::broadcast() { pthread_cond_broadcast(&c);     }

lock::lock(MutexCondition &mc) : mc(mc) { mc.lock(); }
lock::~lock() { mc.unlock(); }

template <typename T>
Queue<T>::Queue(void) { }

template <typename T>
Queue<T>::~Queue(void) { }

template <typename T>
void Queue<T>::push(T work)
{
	lock Lock(QueueMC);
	queue.push(work);
	QueueMC.broadcast();
}

template <typename T>
T Queue<T>::pop()
{
	T work;
	lock Lock(QueueMC);

	while ( queue.empty() ) {
		QueueMC.wait(); }

	work = queue.front();
	queue.pop();

	return(work);
}

// Explicit template instantiation
template class Queue<struct worker_data_s>;
#endif // PTHREADS
