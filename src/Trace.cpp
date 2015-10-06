#include "TraceHBonds.h"
#include "queue.h"
#include "WorkerThreads.h"
#include "cpu.h"
#include "timedoutput.h"
#include "Trace.h"

extern bool THB_VERBOSE;

#ifdef PTHREADS
extern Queue<struct worker_data_s> inQueue;
extern Queue<struct worker_data_s> outQueue;
#endif

typedef std::vector<struct HydrogenBond *> HBVec;
typedef std::vector<struct HydrogenBondIterator_s> HBVecIter;

/*
 * Return values:
 *
 * true   : A new string is found.
 * false  : Not a new string.
 *
 */
bool TraceThread( std::vector<ListOfHBonds *> *HBStrings,
                  struct HydrogenBondIterator_s *HBit)
{
	// DonorO --- Hydrogen ... AcceptorO
	// ... Denotes the Hydrogen bond.

	// The Range of HydrogenBonds to search.
	HBVec::iterator iter_begin = HBit->begin;
	HBVec::iterator iter_end   = HBit->end;

	HBVec::iterator iter_hbmain = iter_begin;
	for( ; iter_hbmain < iter_end; ++iter_hbmain)
	{

		// If this hydrogen bond has already been assigned to a chain, skip it
		if ( ( (*iter_hbmain)->Next != NULL) ||
		     ( (*iter_hbmain)->Previous != NULL) )
			continue;


		// Starting a new chain.
		ListOfHBonds *HBonds = new ListOfHBonds();
		(HBonds)->AddAtBegin(*iter_hbmain);

		bool StillLooking = true;
		while ( StillLooking )
		{
			bool FoundOne = false;
			HBVec::iterator iter_hb;
			for(iter_hb = iter_begin ; iter_hb < iter_end; ++iter_hb )
			{
				// If this hydrogen bond has already been assigned to a
				// chain, skip it
				if ( ( (*iter_hb)->Next != NULL) ||
				     ( (*iter_hb)->Previous != NULL) )
					continue;

				// Check that this hydrogen bond is not in the chain already.
				// Checking for only the H is sufficient.
				if ( !(HBonds)->Find(*iter_hb)  )
				{
					if ( (HBonds)->linksAtBegin(*iter_hb) )
					{
						// Found a new link at the beginning of the chain.
						(HBonds)->AddAtBegin(*iter_hb);
						FoundOne = true;
					}
					else if ( (HBonds)->linksAtEnd(*iter_hb) )
					{
						// Found a new link at the end of the chain.
						(HBonds)->AddAtEnd(*iter_hb);
						FoundOne = true;
					}
				}
			}
			StillLooking = FoundOne;
		}
		HBStrings->push_back(HBonds);
	}

	return(true);
}

void Trace( std::vector<ListOfHBonds *> *HBStrings,
            HBVecIter *TrjIdx_iter )
{
	timedOutput msg(TrjIdx_iter->size()+1, 1.0, "Tracing HB strings: frame ");

	for( unsigned int f=0; f < TrjIdx_iter->size(); ++f )
	{
#ifdef PTHREADS
		struct worker_data_s wd;
		wd.jobtype = THREAD_JOB_TRACE;
		wd.jobnum = f;
		wd.num_threads = NumberOfCPUs();
		wd.HBit = &(TrjIdx_iter->at(f));
		wd.HBStrings = new std::vector<ListOfHBonds *>;
		wd.HBStrings->reserve(5000);

		inQueue.push(wd);
#else
		if ( THB_VERBOSE ) { msg.print(f+1); }

		HBVec::iterator iter_hb = TrjIdx_iter->at(f).begin;
		TraceThread( HBStrings, &TrjIdx_iter->at(f) );
#endif // PTHREADS
	}
#ifdef PTHREADS
	// Get the results back from the worker threads.
	for( unsigned int f=0; f < TrjIdx_iter->size(); ++f )
	{
		if ( THB_VERBOSE ) { msg.print(f+1); }

		struct worker_data_s wd = outQueue.pop();
		HBStrings->reserve( TrjIdx_iter->size()*wd.HBStrings->size() );
		HBStrings->insert( HBStrings->end(),
		                   wd.HBStrings->begin(),
		                   wd.HBStrings->end() );

		wd.HBStrings->clear();
		delete wd.HBStrings;
	}
#endif // PTHREADS
	VERBOSE_MSG("Tracing HB strings: frame " << TrjIdx_iter->size() <<"/"<< TrjIdx_iter->size() << ".");
}
