#include "help.h"
#include "OutputFormat.h"
#include "TraceHBonds.h"
#include "queue.h"
#include "WorkerThreads.h"
#include "cpu.h"
#include "Trace.h"
#include "RemoveDuplicates.h"
#include "sizehist.h"
#include "correlation.h"
#include "neighborhist.h"
#include "deleteVectorPointers.h"

extern bool THB_VERBOSE;

#ifdef PTHREADS
extern Queue<struct worker_data_s> inQueue;
extern Queue<struct worker_data_s> outQueue;
#endif

typedef std::vector<struct HydrogenBond *> HBVec;
typedef std::vector<struct HydrogenBondIterator_s> HBVecIter;

// Used with remove_if()
bool deleteAtom( struct thbAtom *atom ) {
	// delete atom->Molecule;
	delete atom;
	return true;
}

int doArcFile(char *ifilename,
              char *ofPrefix, char *ofSuffix,
              struct HydrogenBondMatching *match,
              double rCutoff, double angleCutoff,
              int NumBins, unsigned char flags)
{
	HBVec hb;
	std::vector<struct thbAtom *> atom;
	unsigned int NumFramesInTrajectory = 0;

	// If neither NEIGHBOR_HIST or SIZE_HIST are specifies, we can try to
	// save some memory by storing the atom coordinates only as long as we
	// need the.
	bool SaveMemory = true;
	if ( flags&(NEIGHBOR_HIST|SIZE_HIST) )
		SaveMemory = false;
	VERBOSE_MSG("SaveMemory: " << SaveMemory);

	atom.reserve(50000);
	DEBUG_MSG("Capacity/size of atom: " << atom.capacity() << "/" << atom.size());
	// Get Atoms and their connections.
	if ( ! ConnectionsMDF( ifilename, &atom ) ) {
		return 1; }

	// Find the Hydrogens and Acceptors.
	std::vector<struct thbAtom *> hydrogens;
	std::vector<struct thbAtom *> acceptors;
	hydrogens.reserve(atom.size()/2);
	acceptors.reserve(atom.size()/2);

	VERBOSE_MSG("");
	getHydrogenBondElements( &atom, &hydrogens, &acceptors, match );

	if ( (flags & (LIFETIME|SIZE_HIST|NEIGHBOR_HIST|LENGTHS)) == 0 )
	{
		DeleteVectorPointers( atom ); atom.clear();
		return(0);
	}

	// VERBOSE_MSG("Finding hydrogen bonds with:\n\n\tRc    < " << rCutoff << " Angstroms, and \n\tangle > " << angleCutoff << " degrees.\n");
	VERBOSE_MSG("Max Distance (A): " << rCutoff);
	VERBOSE_MSG("Min Angle  (deg): " << angleCutoff);
	VERBOSE_MSG("");

	struct PBC *Cell;
	Cell = new struct PBC;

	// Get Atom positions and cell dimensions for each frame.
#ifdef PTHREADS
	struct worker_data_s wd;
	wd.jobtype     = THREAD_JOB_POSITIONS_CAR;
	wd.jobnum      = 1;
	wd.num_threads = NumberOfCPUs();
	wd.filename    = ifilename;
	wd.atom        = &atom;
	wd.hydrogens   = &hydrogens;
	wd.acceptors   = &acceptors;
	wd.Cell        = Cell;
	wd.rCutoff     = rCutoff;
	wd.angleCutoff = angleCutoff;
	wd.saveMemory  = SaveMemory;

	inQueue.push(wd);

	// Get the result back from the worker threads.
	unsigned int FramesProcessed=0;
	time_t timer = time(NULL);
	while ( 1 )
	{
		struct worker_data_s wdOut = outQueue.pop();
		if ( (wdOut.jobtype == THREAD_JOB_HBS2) ||
		     (wdOut.jobtype == THREAD_JOB_HBS ) )
		{
			/** \todo Use a more appropriate value than 5000. */
			hb.reserve(5000*wdOut.hb->size() );
			if ( NumFramesInTrajectory != 0 )
				hb.reserve(NumFramesInTrajectory*wdOut.hb->size() );

			hb.insert(hb.end(),
			          wdOut.hb->begin(),
			          wdOut.hb->end() );
			wdOut.hb->clear();
			delete wdOut.hb;

			if ( wdOut.jobtype == THREAD_JOB_HBS2 ) {
				wdOut.coordinates->clear();
				delete wdOut.coordinates;
			}
			FramesProcessed++;

			// Only show this message if we are done reading all frames.
			if ( THB_VERBOSE && (NumFramesInTrajectory != 0) &&
				 ((difftime(time(NULL), timer) > 1.0) || (FramesProcessed == NumFramesInTrajectory))  )
			{
				VERBOSE_CMSG("Processing frame " << FramesProcessed );
				VERBOSE_CMSG("/" << NumFramesInTrajectory);
				VERBOSE_RMSG(". Hydrogen-acceptor pairs found: " << hb.size() << ".");
				timer = time(NULL);
			}

			if ( FramesProcessed == NumFramesInTrajectory )
				break;
		}
		if ( wdOut.jobtype == THREAD_JOB_POSITIONS_CAR )
		{
			NumFramesInTrajectory = wdOut.TrjIdx;
			if ( FramesProcessed == NumFramesInTrajectory )
				break;
		}
	}
	VERBOSE_MSG("");
#else
	PositionsCAR( ifilename, &atom, Cell, SaveMemory );
#endif //PTHREADS
	// std::vector<double>A, B, C;

	NumFramesInTrajectory = Cell->frames;


	if ( hb.size() == 0 )
	{
		BRIEF_MSG("No hydrogen bonds found.");

		//Cleanup
		delete Cell;
		DeleteVectorPointers( atom ); atom.clear();
		DeleteVectorPointers( hb ); hb.clear();

		return(1);
	}
	// Reserve space for hydrogen bonds to minimize reallocations.
	// This is just an emperical guess for reserve size.
	// hb.reserve((atom.size()/4)*NumFramesInTrajectory);
	DEBUG_MSG("\tCapacity/size of hb: " << hb.capacity() << "/" << hb.size());
	// Now  determine the hydrogen bonds
	// AtomNeighbors( &hb, Cell, &hydrogens, &acceptors,
	//                rCutoff, angleCutoff );
	DEBUG_MSG("Capacity/size of hb: " << hb.capacity() << "/" << hb.size());
	if ( hb.size() != hb.capacity() )
		HBVec(hb).swap(hb);
	DEBUG_MSG("\tCapacity/size of hb: " << hb.capacity() << "/" << hb.size());

	HBVecIter TrjIdx_iter(NumFramesInTrajectory);

	TrajectoryIndexIterator( &TrjIdx_iter, &hb );

	VERBOSE_MSG("Looking for smallest hydrogen-acceptor bond lengths in all frames...");

	RemoveDuplicates ( &hb, &TrjIdx_iter );
	VERBOSE_MSG("Hydrogen bonds: " << hb.size() << ".");
	if ( hb.size() != hb.capacity() )
		HBVec(hb).swap(hb);
	DEBUG_MSG("\tCapacity/size of hb: " << hb.capacity() << "/" << hb.size());

	// Update TrjIdx_iter after removing elements.
	TrajectoryIndexIterator( &TrjIdx_iter, &hb );


	// Hydrogen bond lifetime correlations.
	if ( flags & LIFETIME )
	{
		VERBOSE_MSG("Lifetime of a hydrogen bond.");
		VERBOSE_MSG("\tGenerating truth table.");
		std::vector< std::vector<bool> >correlationData;
		Lifetime(&correlationData, &TrjIdx_iter);

		if ( correlationData.size() != 0) {
			std::stringstream ofilename;
			ofilename << ofPrefix << "-lifetimes" << ofSuffix;

			std::ofstream out;
			out.open(ofilename.str().c_str(),std::ios::out);

			if ( out.is_open() ) {
				Correlations(&out, &correlationData); }
			out.close();
			VERBOSE_MSG("\t\tSaved " << ofilename.str() << ".");
		}
	}

	// Save Hydrogen bond length H...Acceptor, and angle.
	if ( flags & LENGTHS )
	{
		std::stringstream ofilename;

		ofilename << ofPrefix << "-lengths" << ofSuffix;

		std::ofstream out(ofilename.str().c_str(),std::ios::out);

		if ( out.is_open() )
		{
			VERBOSE_MSG("\tSaving Hydrogen bond lengths");

			for( unsigned int i=0; i < hb.size(); i++ ) {
				out << hb.at(i)->length << "\n"; }

			VERBOSE_MSG("\t\tSaved " << ofilename.str() << ".");
		}
	}

	if ( flags & ANGLES )
	{
		std::stringstream ofilename;

		ofilename << ofPrefix << "-angles" << ofSuffix;

		std::ofstream out(ofilename.str().c_str(),std::ios::out);

		if ( out.is_open() )
		{
			VERBOSE_MSG("\tSaving Hydrogen bond angles");

			for( unsigned int i=0; i < hb.size(); i++ ) {
				out << hb.at(i)->angle << "\n"; }

			VERBOSE_MSG("\t\tSaved " << ofilename.str() << ".");
		}
	}


	if ( flags & (SIZE_HIST|NEIGHBOR_HIST) ) {

		// Each element of the vector points to a string of hbonds.
		// ListOfHBonds is a strings of hbonds.
		std::vector<ListOfHBonds *>HBStrings;
		HBStrings.reserve(1000*NumFramesInTrajectory);

		//Find all the strings.
		VERBOSE_RMSG("Tracing HB strings.");
		Trace( &HBStrings, &TrjIdx_iter );

		if ( flags & SIZE_HIST ) {
			sizehist(NumFramesInTrajectory,
			         ofPrefix, ofSuffix,
			         &hb,
			         Cell,
			         NumBins,
			         flags,
			         &HBStrings);
		}

		if ( flags & NEIGHBOR_HIST ) {
			neighborhist(NumFramesInTrajectory,
			             ofPrefix, ofSuffix,
			             Cell,
			             &HBStrings);
		}
		//Cleanup
		DeleteVectorPointers( HBStrings ); HBStrings.clear();
	}

	//Cleanup
	delete Cell;
	DeleteVectorPointers( atom ); atom.clear();
	DeleteVectorPointers( hb ); hb.clear();

	return(0);
}

/*
 * Generate Iterators which point to hydrogen bonds in specific frames.
 * TrjIdx_iter.at(5).begin points to first hydrogen bond (\p hb) in frame 5
 * (frame numbers start at 0).  TrjIdx_iter.at(5).end points to just past the
 * last hydrogen bond in frame 5.  The hydrogen bonds in \p hb are grouped by
 * trajectory index number, however the order of the groups may not be in
 * sequence because they are obtained from threads.
 */
void
TrajectoryIndexIterator( HBVecIter *TrjIdx_iter, HBVec *hb)
{
	HBVec::iterator it_hb;

	// Initialize each to hb->end()
	// If there are no hydrogen bonds in a frame, both
	// TrjIdx_iter[].begin and TrjIdx_iter[].end will be hb->end().
	// Can simply check if TrjIdx_iter.at(i)->begin == TrjIdx_iter.at(i)->end
	// to determine if a frame is void of hydrogen bonds.
	for(unsigned int i=0; i < TrjIdx_iter->size(); ++i)
	{
		TrjIdx_iter->at(i).begin = hb->end();
		TrjIdx_iter->at(i).end   = hb->end();
	}

	unsigned int curIdx=(*hb->begin())->TrajIdx;
	TrjIdx_iter->at(curIdx).begin = hb->begin();

	for(it_hb = hb->begin(); it_hb != hb->end(); ++it_hb)
	{
		if ( (*it_hb)->TrajIdx != curIdx )
		{
			unsigned int newIdx = (*it_hb)->TrajIdx;
			TrjIdx_iter->at(newIdx).begin = it_hb;
			TrjIdx_iter->at(curIdx).end = it_hb;
			curIdx = newIdx;
		}
	}
	TrjIdx_iter->at(curIdx).end = hb->end();
}


