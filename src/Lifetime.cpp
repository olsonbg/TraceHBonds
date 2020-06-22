#include "Lifetime.h"

void
Lifetime(std::vector< std::vector<bool> >*b,  HBVecIter *TrjIdx_iter )
{
	HBVec::iterator iter_hb;
	HBVec::iterator iter_begin;
	HBVec::iterator iter_end;


	unsigned int MaxNumHBsInFrame=0;

	for(unsigned int i=0; i < TrjIdx_iter->size(); ++i) {
		iter_begin = TrjIdx_iter->at(i).begin;
		iter_end   = TrjIdx_iter->at(i).end;

		unsigned int NumHBsInFrame=0;
		for(iter_hb = iter_begin ; iter_hb < iter_end; ++iter_hb ) {
			++NumHBsInFrame; }

		if (NumHBsInFrame > MaxNumHBsInFrame) {
			MaxNumHBsInFrame = NumHBsInFrame; }
	}

	VERBOSE_MSG("\tMaximum number of Hydrogen bonds in a frame: " << MaxNumHBsInFrame);


	if (MaxNumHBsInFrame == 0)
		return;

	unsigned int NumFrames = TrjIdx_iter->size();
	iter_begin = TrjIdx_iter->at( 0 ).begin;
	iter_end   = TrjIdx_iter->at( 0 ).end;

	// Initialize all elements to false.
	b->assign(MaxNumHBsInFrame, std::vector<bool>(NumFrames,false));

#ifdef PTHREADS
	for( unsigned int jobnum=0; jobnum < NumberOfCPUs(); ++jobnum)
	{
		struct worker_data_s wd;
		wd.jobtype = THREAD_JOB_LIFETIME;
		wd.jobnum = jobnum;
		wd.num_threads = NumberOfCPUs();
		wd.TrjIdx_iter = TrjIdx_iter;

		wd.b = new std::vector< std::vector<bool> >(MaxNumHBsInFrame,
		                                            std::vector<bool>(NumFrames,false));

		inQueue.push(wd);
	}

	// Get the results back from the worker threads.
	struct worker_data_s wd[NumberOfCPUs()];
	for( unsigned int jobnum=0; jobnum < NumberOfCPUs(); ++jobnum) {
		wd[jobnum] = outQueue.pop(); }

	// Merge the results, if any b->at(i).at(j) is true, set it.
	for ( unsigned int i=0; i < MaxNumHBsInFrame; ++i ) {
		for ( unsigned int j=0; j < NumFrames; ++j ) {
			for ( unsigned int t=0; t < NumberOfCPUs(); ++t ) {
				if ( wd[t].b->at(i).at(j) ) {
					b->at(i).at(j) = true;
					break;
				}
			}
		}
	}

	for ( unsigned int i=0; i < NumberOfCPUs(); ++i ) {
		delete wd[i].b; }
#else
	LifetimeThread( b, TrjIdx_iter );
#endif // PTHREADS
}

void
LifetimeThread(std::vector< std::vector<bool> >*b,  HBVecIter *TrjIdx_iter,
               unsigned int NumThreads, unsigned int ThreadID)
{
	HBVec::iterator iter_hb;
	HBVec::iterator iter_begin;
	HBVec::iterator iter_end;

	unsigned int NumFrames = TrjIdx_iter->size();
	iter_begin = TrjIdx_iter->at( 0 ).begin;
	iter_end   = TrjIdx_iter->at( 0 ).end;

	HBVec::iterator iter_hbmain = iter_begin;

	unsigned int NumHBsInFrameZero = b->size();

	for( unsigned int h=ThreadID; h < NumHBsInFrameZero; h += NumThreads)
	{
		b->at(h).at(0) = true;
		for( unsigned int f=1; f < NumFrames; ++f )
		{
			// The range of hydrogen bonds in frame f.
			iter_begin = TrjIdx_iter->at( f ).begin;
			iter_end   = TrjIdx_iter->at( f ).end;

			for(iter_hb = iter_begin ; iter_hb < iter_end; ++iter_hb )
			{
				if ( ((*iter_hb)->hydrogen == (*(iter_hbmain+h))->hydrogen) &&
				     ((*iter_hb)->acceptor == (*(iter_hbmain+h))->acceptor) )
				{
					// This matched the hydrogen bond in the first frame.
					b->at(h).at(f) = true;
					break;
				}
			}
		}
	}
}
