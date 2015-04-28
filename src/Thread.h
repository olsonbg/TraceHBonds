#ifndef _Thread_h
#define _Thread_h

#ifdef PTHREADS

#include <pthread.h>

class Thread
{
	private:
		pthread_t tid;
		bool running;
		bool detached;
	public:
		Thread();
		virtual ~Thread();

		int start();
		int join();
		int detach();
		pthread_t self();
		virtual void *run() = 0;
};

class MyThread : public Thread {
	public:
		void *run();
};

#endif // PTHREADS
#endif // _Thread_h
