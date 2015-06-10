#include "RemoveDuplicates.h"

bool SameAtom( struct thbAtom *A,
               struct thbAtom *B) {
	return A == B; }

// Used with remove_if()
bool deleteMarked( struct HydrogenBond *hb ) {
	if ( hb->markedDuplicate )
	{
		delete hb;
		return true;
	}
	return false;
}

void removeMarked( HBVec *hb )
{
	// TODO: Why does remove_if not work when cross compiling for windows? It
	// causes the program to crash. I'll have to look into this, but for now
	// use a very simple alternative (which is not very efficient).
#ifdef _WIN32
	HBVec::iterator iter_hb = hb->begin();
	for( iter_hb = hb->end()-1; iter_hb >= hb->begin(); --iter_hb)
	{
		if ( (*iter_hb)->markedDuplicate )
		{
			delete *iter_hb;
			iter_hb = hb->erase(iter_hb);
		}
	}
#else
	HBVec::iterator pos;

	pos = remove_if(hb->begin(), hb->end(), deleteMarked);

	if ( pos != hb->end() )
	{
		// HBVec::iterator toDelete;
		// for( toDelete=pos; toDelete != hb->end(); toDelete++ ) {
		//     delete *toDelete; }

		hb->erase( pos, hb->end() );
	}
#endif
}

void RemoveDuplicatesThread( struct HydrogenBondIterator_s HBit )
{
	HBVec::iterator iter_hbmain;
	HBVec::iterator iter_hb;

	/*
	 * Look for acceptor duplicates
	 */

	for( iter_hbmain = HBit.begin; iter_hbmain < HBit.end-1; ++iter_hbmain )
	{
		// If this is already marked as a duplicate, skip it.
		if ( (*iter_hbmain)->markedDuplicate )
			continue;

		// Store all hydrogen bonds which have the same acceptor atom as
		// iter_hbmain.
		std::vector< HBVec::iterator > Neighbors;
		Neighbors.push_back( iter_hbmain );

		for( iter_hb = iter_hbmain+1; iter_hb < HBit.end; ++iter_hb )
		{
			// If this is already marked as a duplicate, skip it.
			if ( (*iter_hb)->markedDuplicate )
				continue;

			if ( SameAtom( (*iter_hbmain)->acceptor, (*iter_hb)->acceptor) ) {
				Neighbors.push_back( iter_hb ); }
		}

		// Sort the Neighbors vector by length, smallest first.
		std::sort( Neighbors.begin(), Neighbors.end(), sort_Neighbors());

		// Mark all but the shortest length hydrogen bond as duplicate
		unsigned int HydrogenBondMax = (*iter_hbmain)->acceptor->HydrogenBondMax;
		for( unsigned int i=HydrogenBondMax; i < Neighbors.size(); ++i) {
				(*Neighbors[i])->markedDuplicate = true; }
	}

	//
	// Look for H duplicates
	//

	for(iter_hbmain = HBit.begin; iter_hbmain < HBit.end-1; ++iter_hbmain )
	{
		// If this is already marked as a duplicate, skip it.
		if ( (*iter_hbmain)->markedDuplicate )
			continue;

		// Store all hydrogen bonds which have the same hydrogen atom as
		// iter_hbmain.
		std::vector< HBVec::iterator > Neighbors;
		Neighbors.push_back( iter_hbmain );

		for( iter_hb = iter_hbmain+1 ; iter_hb < HBit.end; ++iter_hb )
		{
			// If this is already marked as a duplicate, skip it.
			if ( (*iter_hb)->markedDuplicate )
				continue;

			if ( SameAtom( (*iter_hbmain)->hydrogen, (*iter_hb)->hydrogen) ) {
				Neighbors.push_back( iter_hb ); }
		}

		// Sort the Neighbors vector by length, smallest first.
		std::sort( Neighbors.begin(), Neighbors.end(), sort_Neighbors());

		// Mark all but the shortest length hydrogen bond as duplicate
		unsigned int HydrogenBondMax = (*iter_hbmain)->hydrogen->HydrogenBondMax;
		for( unsigned int i=HydrogenBondMax; i < Neighbors.size(); ++i) {
				(*Neighbors[i])->markedDuplicate = true; }
	}

	return;
}

/*
 * Possible to have more than one:
 *
 *   - acceptor Oxygen (aO) with a single Hydrogen (H)
 *   - Hydrogen (H) with a single acceptor Oxygen (aO)
 *
 * Find the all duplicates and keep the shortest Hydrogen Bond length.
 */

void RemoveDuplicates( HBVec *hb,
                       HBVecIter *TrjIdx_iter)
{
	unsigned int emptyFrames=0;

	for(unsigned int i=0; i != TrjIdx_iter->size(); ++i) {
		if ( TrjIdx_iter->at(i).begin == TrjIdx_iter->at(i).end ) {
			// This frame is empty, skip it.
			emptyFrames++;
			continue; }
#ifdef PTHREADS
		struct worker_data_s wd;
		wd.jobtype = THREAD_JOB_RMDUPS;
		wd.jobnum = i;
		wd.num_threads = NumberOfCPUs();
		wd.HBit = &(TrjIdx_iter->at(i));
		inQueue.push(wd);
#else
		if (  ((i+1)%50==0) || ((i+1)==TrjIdx_iter->size())  )
			VERBOSE_RMSG("Processing frame " << i+1 <<"/"<< TrjIdx_iter->size() << ".");

		RemoveDuplicatesThread(TrjIdx_iter->at(i));
#endif
	}

#ifdef PTHREADS
	// Get the results back from the worker threads.
	for(unsigned int i=0; i != TrjIdx_iter->size()-emptyFrames; ++i) {
		if (  ((i+1)%50==0) || ((i+1)==TrjIdx_iter->size())  )
			VERBOSE_RMSG("Processing frame " << i+1 <<"/"<< TrjIdx_iter->size() << ".");

		// Do not need to do anything with result of outQueue.pop,
		// just have to wait for it to complete.
		outQueue.pop();
	}
#endif // PTHREADS
	VERBOSE_MSG("Processing frame " << TrjIdx_iter->size() <<"/"<< TrjIdx_iter->size() << ".");

	VERBOSE_MSG("Removing Marked");
	removeMarked(hb);
	return;

}
