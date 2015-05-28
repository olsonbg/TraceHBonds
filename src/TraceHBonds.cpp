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
              int NumBins, bool POVRAY)
{
	HBVec hb;
	std::vector<struct thbAtom *> atom;

	std::map<std::string,double> times;
	time_t t_start;

	struct PBC *Cell;
	Cell = new struct PBC;

	atom.reserve(50000);
	DEBUG_MSG("Capacity/size of atom: " << atom.capacity() << "/" << atom.size());
	t_start = time(NULL);
#ifdef PTHREADS
	struct worker_data_s wd;
	wd.jobtype = THREAD_JOB_READCARMDF;
	wd.jobnum = 1;
	wd.num_threads = NumberOfCPUs();
	wd.filename = ifilename;
	wd.atom = &atom;
	wd.Cell = Cell;

	inQueue.push(wd);

	// Get the result back from the worker threads.
	outQueue.pop();
#else
	ReadCarMdf( ifilename, &atom, Cell );
#endif //PTHREADS
	times["reading data"]  = difftime(t_start, time(NULL));
	// Free up some space, if possible.
	if ( atom.size() != atom.capacity() )
		std::vector<struct thbAtom *>(atom).swap(atom);
	DEBUG_MSG("Capacity/size of atom: " << atom.capacity() << "/" << atom.size());
	// std::vector<double>A, B, C;

	unsigned int NumFramesInTrajectory = 0;
	NumFramesInTrajectory = Cell->frames;

	VERBOSE_MSG("Total frames: " << NumFramesInTrajectory);

	// Reserve space for hydrogen bonds to minimize reallocations.
	// This is just an emperical guess for reserve size.
	hb.reserve((atom.size()/4)*NumFramesInTrajectory);
	DEBUG_MSG("\tCapacity/size of hb: " << hb.capacity() << "/" << hb.size());
	// Now  determine the hydrogen bonds
	t_start = time(NULL);
	AtomNeighbors( &hb, &atom, Cell, match, rCutoff, angleCutoff );
	times["finding pairs"] = difftime(t_start, time(NULL));
	DEBUG_MSG("Capacity/size of hb: " << hb.capacity() << "/" << hb.size());
	if ( hb.size() != hb.capacity() )
		HBVec(hb).swap(hb);
	DEBUG_MSG("\tCapacity/size of hb: " << hb.capacity() << "/" << hb.size());

	struct HydrogenBondIterator_s HBit;
	HBVecIter TrjIdx_iter(NumFramesInTrajectory);

	TrajectoryIndexIterator( &TrjIdx_iter, &hb );

	// for(unsigned int i=0; i != TrjIdx_iter.size(); ++i)
	// {
	//     VERBOSE_MSG("TrjIdx_iter.begin:" << std::distance(hb.begin(), TrjIdx_iter.at(i).begin) );
	//     VERBOSE_MSG("TrjIdx_iter.end  :" << std::distance(TrjIdx_iter.at(i).begin, TrjIdx_iter.at(i).end) );
	// }

	VERBOSE_MSG("Looking for smallest hydrogen-acceptor bond lengths in all frames...");

	t_start = time(NULL);
	RemoveDuplicates ( &hb, &TrjIdx_iter );
	times["finding hydrogen bonds"] = difftime(t_start, time(NULL));
	VERBOSE_MSG("Hydrogen bonds:          " << hb.size() << ".");
	if ( hb.size() != hb.capacity() )
		HBVec(hb).swap(hb);
	DEBUG_MSG("\tCapacity/size of hb: " << hb.capacity() << "/" << hb.size());

	// Update TrjIdx_iter after removing elements.
	TrajectoryIndexIterator( &TrjIdx_iter, &hb );

	// for(unsigned int i=0; i != TrjIdx_iter.size(); ++i)
	// {
	//     VERBOSE_MSG("TrjIdx_iter.begin:" << std::distance(hb.begin(), TrjIdx_iter.at(i).begin) );
	//     VERBOSE_MSG("TrjIdx_iter.end  :" << std::distance(TrjIdx_iter.at(i).begin, TrjIdx_iter.at(i).end) );
	// }

	unsigned int TrjIdx;

	// Each element of the vector points to a string of hbonds.
	// ListOfHBonds is a strings of hbonds.
	std::vector<ListOfHBonds *>HBStrings;
	HBStrings.reserve(1000*NumFramesInTrajectory);

	//Find all the strings.
	VERBOSE_MSG("Tracing HB strings.");

	t_start = time(NULL);
	Trace( &HBStrings, &TrjIdx_iter );
	times["finding hydrogen bond strings"] = difftime(t_start, time(NULL));




	// Hydrogen bond lifetime correlations.
	if ( 1 )
	{
		VERBOSE_MSG("Lifetime of a hydrogen bond.");
		VERBOSE_MSG("\tGenerating truth table.");
		std::vector< std::vector<bool> >correlationData;
		Lifetime(&correlationData, &TrjIdx_iter);
		VERBOSE_MSG("\t\tFinished truth table.");

		std::ofstream out;
		out.open("Correlations.txt",std::ios::out);
		VERBOSE_MSG("\tGenerating Correlations.");
		if ( out.is_open() ) {
			Correlations(&out, &correlationData); }

		out.close();
	}

	// Save Hydrogen bond length H...Acceptor, and angle.
	if ( 0 )
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


	if ( 0 ) {
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
		if ( POVRAY )
			CC = CC2;
		else
			CC = CC1;

		t_start = time(NULL);
		for( TrjIdx = 0 ; TrjIdx < NumFramesInTrajectory; ++TrjIdx )
		{
			if (  ((TrjIdx+1)%50==0) || ((TrjIdx+1)==NumFramesInTrajectory)  )
				VERBOSE_RMSG("Saving size histograms, frame " << TrjIdx+1 << "/" << NumFramesInTrajectory);

			std::stringstream ofilename;
			ofilename << ofPrefix << TrjIdx+1 << ofSuffix;

			std::ofstream out;
			out.open(ofilename.str().c_str(),std::ios::out);
			if ( out.is_open() )
			{
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
				prntHistograms( &out, HBStrings, &Histograms.at(TrjIdx), CC, NumBins, Cell, TrjIdx, POVRAY);

				out.close();
			}
		}

		VERBOSE_MSG("Saving neighbor histograms.");
		if ( 0 )
		{
			VERBOSE_MSG("Generating neighbor histograms.");
			for( TrjIdx = 0 ; TrjIdx < NumFramesInTrajectory; ++TrjIdx ) {
				getNeighbors( &(Histograms.at(TrjIdx)), HBStrings, Cell );}

			std::ofstream out;
			out.open("Neighbors.txt",std::ios::out);
			if ( out.is_open() )
			{
				Print_AllFrames(&out, &Histograms);
				Print_CombineFrames(&out, &Histograms);
				Print_CombineNeighbors(&out, &Histograms);
			}

			out.close();
		}
	}
	times["Saving files"] = difftime(t_start, time(NULL));

	VERBOSE_MSG("\n\nTime spend:");
	std::map<std::string,double>::iterator time_it;
	for( time_it = times.begin(); time_it!=times.end(); ++time_it)
	{
		double t = -1.0*time_it->second;

		VERBOSE_CMSG("  " << time_it->first << " - ");
		if ( t > 60.0 )
			VERBOSE_MSG(t/60.0 << " min.");
		else
			VERBOSE_MSG(t << " sec.");
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


