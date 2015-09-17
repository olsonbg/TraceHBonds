#ifdef PTHREADS

#include "Thread.h"

Thread::Thread() : tid(0), running(false), detached(false) {}
Thread::~Thread()
{
	if ( running && !detached ) {
		pthread_detach(tid); }
	else if ( running ) {
		pthread_cancel(tid); }
}

static void *runThread(void *arg);
int Thread::start()
{
	int result = pthread_create( &tid, NULL, runThread, this);
	if ( result == 0 ) {
		running = true; }

	return(result);
}

// run() is defined in WorkerThreads.cpp
static void *runThread(void *arg)
{
	return ((Thread*)arg)->run();
}

int Thread::join()
{
	int result = -1;
	if ( running )
	{
		result = pthread_join(tid, NULL);
		if ( result == 0 ) {
			running = false; }
	}
	return result;
}

int Thread::detach()
{
	int result = -1;
	if ( running && !detached )
	{
		result = pthread_detach(tid);
		if ( result == 0 ) {
			detached = 1; }
	}
	return result;
}

pthread_t Thread::self() {
	return tid;
}
#endif // PTHREADS
