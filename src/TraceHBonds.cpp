#include "Print.h"
#include "OutputFormat.h"
#include "TraceHBonds.h"
#include "queue.h"
#include "WorkerThreads.h"
#include "Histograms.h"
#include "NeighborPrint.h"
#include "cpu.h"
#include "Trace.h"
#include "RemoveDuplicates.h"

extern bool THB_VERBOSE;

#ifdef PTHREADS
extern Queue<struct worker_data_s> inQueue;
extern Queue<struct worker_data_s> outQueue;
#endif

typedef std::vector<struct HydrogenBond *> HBVec;
typedef std::vector<struct HydrogenBondIterator_s> HBVecIter;

// Used with remove_if()
template <typename T> bool deleteVectorPointers( T* v ) {
	delete v;
	return true;
}

template<class T> void DeleteVectorPointers( std::vector<T *> v )
{
	remove_if(v.begin(),v.end(),deleteVectorPointers<T>);
	// for(unsigned int i =0; i < v.size(); ++i)
	//     delete v[i];
}

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

	struct PBC *Cell;
	Cell = new struct PBC;

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
	ConnectionsMDF( ifilename, &atom );

	// Find the Hydrogens and Acceptors.
	std::vector<struct thbAtom *> hydrogens;
	std::vector<struct thbAtom *> acceptors;
	hydrogens.reserve(atom.size()/2);
	acceptors.reserve(atom.size()/2);

	getHydrogenBondElements( &atom, &hydrogens, &acceptors, match );

	if ( (flags & (LIFETIME|SIZE_HIST|NEIGHBOR_HIST|LENGTHS)) == 0 )
	{
		DeleteVectorPointers( atom ); atom.clear();
		return(0);
	}


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
	while ( 1 )
	{
		struct worker_data_s wdOut = outQueue.pop();
		if ( (wdOut.jobtype == THREAD_JOB_HBS2) ||
		     (wdOut.jobtype == THREAD_JOB_HBS ) )
		{
			// TODO: Use a more appropriate value than 5000.
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
			if ( (NumFramesInTrajectory != 0) &&
				 ((FramesProcessed%10==0) || (FramesProcessed == NumFramesInTrajectory))  )
			{
				VERBOSE_CMSG("Processing frame " << FramesProcessed );
				VERBOSE_CMSG("/" << NumFramesInTrajectory);
				VERBOSE_RMSG(". Hydrogen-acceptor pairs found: " << hb.size() << ".");
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

	struct HydrogenBondIterator_s HBit;
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

	unsigned int TrjIdx;

	// Each element of the vector points to a string of hbonds.
	// ListOfHBonds is a strings of hbonds.
	std::vector<ListOfHBonds *>HBStrings;
	HBStrings.reserve(1000*NumFramesInTrajectory);

	//Find all the strings.
	VERBOSE_RMSG("Tracing HB strings.");
	Trace( &HBStrings, &TrjIdx_iter );

	// Hydrogen bond lifetime correlations.
	if ( flags & LIFETIME )
	{
		VERBOSE_MSG("Lifetime of a hydrogen bond.");
		VERBOSE_MSG("\tGenerating truth table.");
		std::vector< std::vector<bool> >correlationData;
		Lifetime(&correlationData, &TrjIdx_iter);
		VERBOSE_MSG("\tFinished truth table.");

		std::ofstream out;
		out.open("Correlations.txt",std::ios::out);
		VERBOSE_MSG("\tGenerating Correlations.");
		if ( out.is_open() ) {
			Correlations(&out, &correlationData); }

		out.close();
	}

	// Save Hydrogen bond length H...Acceptor, and angle.
	if ( flags & LENGTHS )
	{
		std::ofstream out;
		out.open("Lengths-Angles.txt",std::ios::out);
		if ( out.is_open() )
		{
			for( unsigned int i=0; i < hb.size(); i++ ) {
				out << hb.at(i)->length << "\n"; }
			out << "\n";
			for( unsigned int i=0; i < hb.size(); i++ ) {
				out << hb.at(i)->angle << "\n"; }
		}

		out.close();
	}


	if ( flags & SIZE_HIST ) {
		// Make histograms.
		VERBOSE_MSG("Generating size histograms.");
		std::vector<struct Histograms_s> Histograms;
		for( TrjIdx = 0 ; TrjIdx < NumFramesInTrajectory; ++TrjIdx ) {
			Histograms.push_back( makeHistograms(HBStrings, TrjIdx) ); }

		// Save the histograms.
		const char *CC1 = "#";
		const char *CC2 = "//";
		std::string CC;

		// Povray uses a different comment string.
		if ( flags & POVRAY )
			CC = CC2;
		else
			CC = CC1;

		for( TrjIdx = 0 ; TrjIdx < NumFramesInTrajectory; ++TrjIdx )
		{
			std::stringstream ofilename;
			ofilename << ofPrefix << TrjIdx+1 << ofSuffix;


			std::ofstream out;
			out.open(ofilename.str().c_str(),std::ios::out);
			if ( out.is_open() )
			{
				if ( ((TrjIdx+1)%50==0) || ((TrjIdx+1)==NumFramesInTrajectory) )
					BRIEF_RMSG("Saving size histograms, frame " << TrjIdx+1 << "/" << NumFramesInTrajectory << " (" << ofilename.str() << ")");
				// Header
				out << CC
					<< " PBC "
					<< Cell->p.at(TrjIdx).x()      << " "
					<< Cell->p.at(TrjIdx).y()      << " "
					<< Cell->p.at(TrjIdx).z()      << " "
					<< Cell->angles.at(TrjIdx).x() << " "
					<< Cell->angles.at(TrjIdx).y() << " "
					<< Cell->angles.at(TrjIdx).z()
					<< "\n";

				out <<CC<< " Donor Oxygen atoms    : " << hb.size() << "\n";
				out <<CC<< " Hydrogen atoms        : " << hb.size() << "\n";
				out <<CC<< " Acceptor Oxygen atoms : " << hb.size() << "\n";

				// print the histograms and chains.
				prntHistograms( &out, HBStrings, &Histograms.at(TrjIdx), CC, NumBins, Cell, TrjIdx, flags & POVRAY);

				out.close();
			} else {
				BRIEF_MSG("ERROR: Can not save " <<ofilename.str() << "!" );
			}
		}
		BRIEF_MSG("");

		if ( flags & NEIGHBOR_HIST )
		{
			VERBOSE_MSG("Generating neighbor histograms.");
			for( TrjIdx = 0 ; TrjIdx < NumFramesInTrajectory; ++TrjIdx ) {
				getNeighbors( &(Histograms.at(TrjIdx)), HBStrings, Cell );}

			std::ofstream out;
			out.open("Neighbors.txt",std::ios::out);
			if ( out.is_open() )
			{
				BRIEF_MSG("Saving neighbor histograms: " << "Neighbors.txt" << ".");
				Print_AllFrames(&out, &Histograms);
				Print_CombineFrames(&out, &Histograms);
				Print_CombineNeighbors(&out, &Histograms);
			} else {
				BRIEF_MSG("ERROR: Can not save to " << "Neighbors.txt" << "!"); }

			out.close();
		}
	}

	//Cleanup
	delete Cell;
	DeleteVectorPointers( atom ); atom.clear();
	DeleteVectorPointers( hb ); hb.clear();
	DeleteVectorPointers( HBStrings ); HBStrings.clear();

	return(0);
}

// Save Iterators which point to just past the end of a Trajectory Index.
// TrjIdx_iter.at(1).begin points to first element of TrjIdx 1.
// TrjIdx_iter.at(1).end points to just past the last element of TrjIdx 1 The
// hbs are grouped by trajectory index number, however the order of the groups
// may not be in sequence.

void
TrajectoryIndexIterator( HBVecIter *TrjIdx_iter, HBVec *hb)
{
	// struct HydrogenBondIterator_s HBit;
	// HBit.begin = hb->begin();
	// HBVecIterTrjIdx_iter(N+1,HBit);

	
	unsigned int curIdx=(*hb->begin())->TrajIdx;
	TrjIdx_iter->at(curIdx).begin = hb->begin();

	HBVec::iterator it_hb;
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


