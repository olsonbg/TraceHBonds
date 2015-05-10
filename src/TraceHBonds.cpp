#include "Print.h"
#include "OutputFormat.h"
#include "TraceHBonds.h"
#include "queue.h"
#include "WorkerThreads.h"
#include "Histograms.h"
#include "NeighborPrint.h"
#include "cpu.h"

extern bool THB_VERBOSE;

#ifdef PTHREADS
extern Queue<struct worker_data_s> inQueue;
extern Queue<struct worker_data_s> outQueue;
#endif

void getHydrogenBondElements( std::vector<struct thbAtom *> *atom,
                              std::vector<struct thbAtom *> *hydrogendonors,
                              std::vector<struct thbAtom *> *acceptors,
                              struct HydrogenBondMatching *match)
{
	std::vector<struct thbAtom *>::iterator it_a1;

	for( it_a1 = atom->begin(); it_a1 < atom->end(); ++it_a1)
	{
		if ( match->Hydrogens.find((*it_a1)->ForceField) != match->Hydrogens.end())
			hydrogendonors->push_back(*it_a1);
		else if( match->Acceptors.find((*it_a1)->ForceField) != match->Acceptors.end())
			acceptors->push_back(*it_a1);
	}
}

void HBs( std::vector<struct HydrogenBond *> *hb,
          Point cell,
          std::vector<struct thbAtom *>*hydrogens,
          std::vector<struct thbAtom *>*acceptors,
          double TrjIdx,
          double rCutoff, double angleCutoff)
{

	std::vector<struct thbAtom *>::iterator it_h;
	std::vector<struct thbAtom *>::iterator it_a;

	Point r;

	double rCutoff2 = pow(rCutoff,2.0);
	double r2;

	Point a,b,c;

	for( it_h = hydrogens->begin(); it_h < hydrogens->end(); ++it_h)
	{
		a = (*it_h)->p.at(TrjIdx);

		for( it_a = acceptors->begin(); it_a < acceptors->end(); ++it_a)
		{
			// Make sure this acceptor is not covalently bonded to the hydrogen
			// we are looking at.  The angle check would catch this, but this
			// will skip a few calculations.
			if ( (*it_a) == (*it_h)->ConnectedAtom.at(0) )
				continue;

			b = (*it_a)->p.at(TrjIdx);

			r = a.minimumImage( b, cell );
			r2 = r.magnitudeSquared();

			if ( r2 < rCutoff2)
			{
				// Distance cutoff is good, now check the angle.
				c = (*it_h)->ConnectedAtom.at(0)->p.at(TrjIdx);

				double angle = a.angle(b,c);
				if ( angle > angleCutoff )
				{
					struct HydrogenBond *NewHB;
					NewHB = new struct HydrogenBond;

					NewHB->length   = sqrt(r2);
					NewHB->angle    = angle;
					NewHB->hydrogen = *it_h;
					NewHB->acceptor = *it_a;
					NewHB->donor    = (*it_h)->ConnectedAtom.at(0);
					NewHB->TrajIdx  = TrjIdx;

					NewHB->acceptorDonorDistance=c.minimumImageDistance(b,cell);

					hb->push_back(NewHB);
				}
			}
		}
	}
}

void AtomNeighbors( std::vector<struct HydrogenBond *> *hb,
                    std::vector<struct thbAtom *> *atom,
                    struct PBC *Cell, struct HydrogenBondMatching *match,
                    double rCutoff, double angleCutoff )
{
	std::vector<struct thbAtom *> hydrogens;
	std::vector<struct thbAtom *> acceptors;

	hydrogens.reserve(5000);
	acceptors.reserve(5000);

	getHydrogenBondElements( atom, &hydrogens, &acceptors, match );

	if (THB_VERBOSE)
	{
		std::set<std::string>::iterator it;

		VERBOSE_CMSG("Total hydrogen donors: " << hydrogens.size() << " [ ");
		for(it=match->Hydrogens.begin(); it != match->Hydrogens.end();++it)
			VERBOSE_CMSG(*it << " ");
		VERBOSE_MSG("]");

		VERBOSE_CMSG("Total acceptors      : " << acceptors.size() << " [ ");
		for(it=match->Acceptors.begin(); it != match->Acceptors.end();++it)
			VERBOSE_CMSG(*it << " ");
		VERBOSE_MSG("]");

		VERBOSE_MSG("Finding hydrogen bonds with Rc < " << rCutoff << " Angstroms, and angle > " << angleCutoff << " degrees.");
	}

	unsigned int NumFramesInTrajectory = 0;
	NumFramesInTrajectory = Cell->frames;

	Point cell;
	for( unsigned int TrjIdx=0; TrjIdx < NumFramesInTrajectory; ++TrjIdx)
	{
		cell = Cell->p.at(TrjIdx);

#ifdef PTHREADS
		struct worker_data_s wd;
		wd.jobtype = THREAD_JOB_HBS;
		wd.jobnum = TrjIdx;
		wd.num_threads = NumberOfCPUs();
		wd.cell = cell;
		wd.hydrogens = &hydrogens;
		wd.acceptors = &acceptors;
		wd.TrjIdx = TrjIdx;
		wd.rCutoff = rCutoff;
		wd.angleCutoff = angleCutoff;

		wd.hb = new std::vector<struct HydrogenBond *>;
		wd.hb->reserve(5000);

		inQueue.push(wd);

#else
		HBs( hb, cell, &hydrogens, &acceptors, TrjIdx, rCutoff, angleCutoff);
#endif
	}

#ifdef PTHREADS
	// Get the results back from the worker threads.
	for( unsigned int TrjIdx=0; TrjIdx < NumFramesInTrajectory; ++TrjIdx)
	{
		if (  ((TrjIdx+1)%10==0) || ((TrjIdx+1)==NumFramesInTrajectory)  )
			VERBOSE_RMSG("Processing frame " << TrjIdx+1 <<"/"<< Cell->frames << ". Hydrogen-acceptor pairs found: " << hb->size() << ".");

		struct worker_data_s wd = outQueue.pop();
		hb->reserve( hb->size() + wd.hb->size() );
		hb->insert(hb->end(),
		           wd.hb->begin(),
		           wd.hb->end() );

		wd.hb->clear();
		delete wd.hb;
	}
	VERBOSE_MSG("Processing frame " << Cell->frames <<"/"<< Cell->frames << ". Hydrogen-acceptor pairs found: " << hb->size() << ".");

#endif // PTHREADS
}

int doArcFile(char *ifilename,
              char *ofPrefix, char *ofSuffix,
              struct HydrogenBondMatching *match,
              double rCutoff, double angleCutoff,
              int NumBins, bool POVRAY)
{
	std::vector<struct HydrogenBond *> hb;
	std::vector<struct thbAtom *> atom;

	std::map<std::string,double> times;
	time_t t_start;

	struct PBC *Cell;
	Cell = new struct PBC;

	atom.reserve(50000);
	t_start = time(NULL);
	ReadCarMdf( ifilename, &atom, Cell );
	times["reading data"]  = difftime(t_start, time(NULL));
	std::vector<double>A, B, C;

	unsigned int NumFramesInTrajectory = 0;
	NumFramesInTrajectory = Cell->frames;

	VERBOSE_MSG("Total frames: " << NumFramesInTrajectory);

	// Reserve space for 10000 hydrogen bonds to minimize reallocations.
	hb.reserve(10000);
	// Now  determine the hydrogen bonds
	t_start = time(NULL);
	AtomNeighbors( &hb, &atom, Cell, match, rCutoff, angleCutoff );
	times["finding pairs"] = difftime(t_start, time(NULL));

	std::vector< struct HydrogenBondIterator_s >TrjIdx_iter;
	TrjIdx_iter = TrajectoryIndexIterator( &hb, NumFramesInTrajectory );

	VERBOSE_MSG("Looking for smallest hydrogen-acceptor bond lengths in all frames...");

	t_start = time(NULL);
	RemoveDuplicates ( &hb, &TrjIdx_iter );
	times["finding hydrogen bonds"] = difftime(t_start, time(NULL));
	VERBOSE_MSG("Hydrogen bonds:          " << hb.size() << ".");

	// Update TrjIdx_iter after removing elements.
	TrjIdx_iter = TrajectoryIndexIterator( &hb, NumFramesInTrajectory );

	unsigned int TrjIdx;

	// Each element of the vector points to a string of hbonds.
	// ListOfHBonds is a strings of hbonds.
	std::vector<ListOfHBonds *>HBStrings;
	HBStrings.reserve(5000);

	//Find all the strings.
	VERBOSE_MSG("Tracing HB strings.");

	t_start = time(NULL);
	for( unsigned int i=0; i < hb.size(); i++ )
	{
		ListOfHBonds *HBonds = new ListOfHBonds();
		if ( Trace( &HBonds, &TrjIdx_iter, hb.begin()+i) )
			HBStrings.push_back(HBonds);
		else
			delete HBonds;
	}
	times["finding hydrogen bond strings"] = difftime(t_start, time(NULL));

	// Make histograms.
	VERBOSE_MSG("Generating size histograms.");
	std::vector<struct Histograms_s> Histograms;
	for( TrjIdx = 0 ; TrjIdx < NumFramesInTrajectory; ++TrjIdx ) {
		Histograms.push_back( makeHistograms(HBStrings, TrjIdx) ); }

	VERBOSE_MSG("Generating neighbor histograms.");
	for( TrjIdx = 0 ; TrjIdx < NumFramesInTrajectory; ++TrjIdx ) {
		getNeighbors( &(Histograms.at(TrjIdx)), HBStrings, Cell );}

	VERBOSE_MSG("Lifetime of a hydrogen bond.");
	std::vector< std::vector<bool> >correlationData;
	correlationData = Lifetime(&TrjIdx_iter, hb.begin(), NumFramesInTrajectory);

	if ( 1 )
	{
		std::ofstream out;
		out.open("Correlations.txt",std::ios::out);
		if ( out.is_open() ) {
			Correlations(&out, &correlationData); }

		out.close();
	}

	// Save Hydrogen bond length H...Acceptor, and angle.
	if ( 1 )
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


	VERBOSE_MSG("Saving neighbor histograms.");
	if ( 1 )
	{
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
	DeleteVectorPointers(atom);
	atom.clear();

	DeleteVectorPointers(hb);
	hb.clear();

	for(unsigned int i=0; i < HBStrings.size(); ++i)
		delete HBStrings[i];
	HBStrings.clear();

	return(0);
}

// Save Iterators which point to just past the end of a Trajectory Index.
// TrjIdx_iter.at(1).begin points to first element of TrjIdx 1.
// TrjIdx_iter.at(1).end points to just past the last element of TrjIdx 1 The
// hbs are grouped by trajectory index number, however the order of the groups
// may not be in sequence.
std::vector< struct HydrogenBondIterator_s >
TrajectoryIndexIterator( std::vector<struct HydrogenBond *> *hb,
                         unsigned int N)
{
	struct HydrogenBondIterator_s HBit;
	HBit.begin = hb->begin();
	std::vector<struct HydrogenBondIterator_s>TrjIdx_iter(N,HBit);

	unsigned int curIdx=(*hb->begin())->TrajIdx;

	std::vector<struct HydrogenBond *>::iterator it_hb;
	for(it_hb = hb->begin(); it_hb < hb->end(); ++it_hb)
	{
		if ( (*it_hb)->TrajIdx != curIdx )
		{
			unsigned int newIdx = (*it_hb)->TrajIdx;
			TrjIdx_iter.at(newIdx).begin = it_hb;
			TrjIdx_iter.at(curIdx).end = it_hb;
			curIdx = (*it_hb)->TrajIdx;
		}
	}
	TrjIdx_iter.at(curIdx).end = hb->end();

	return(TrjIdx_iter);
}

template<class T> void DeleteVectorPointers( T v )
{
	for(unsigned int i =0; i < v.size(); ++i)
		delete v[i];
}

bool marked( struct HydrogenBond *hb )
{
	if ( hb->markedDuplicate )
		return true;

	return false;
}

void removeMarked( std::vector<struct HydrogenBond *> *hb )
{
	std::vector<struct HydrogenBond *>::iterator pos;

	pos = remove_if(hb->begin(), hb->end(), marked);

	std::vector<struct HydrogenBond *>::iterator toDelete;
	for( toDelete=pos; toDelete != hb->end(); toDelete++ ) {
		delete *toDelete; }

	hb->erase(pos, hb->end() );
}

/*
 * TODO: Make this work with threads.
 * Possible to have more than one:
 *
 *   - acceptor Oxygen (aO) with a single Hydrogen (H)
 *   - Hydrogen (H) with a single acceptor Oxygen (aO)
 *
 * Find the all duplicates and keep the shortest Hydrogen Bond length.
 */

void RemoveDuplicates( std::vector<struct HydrogenBond *> *hb,
                       std::vector<struct HydrogenBondIterator_s> *TrjIdx_iter)
{

	for(unsigned int i=0; i < TrjIdx_iter->size()-1; ++i) {
		RemoveDuplicatesThread(TrjIdx_iter->at(i));
	}

	VERBOSE_MSG("Removing Marked");
	removeMarked(hb);
	VERBOSE_MSG("Removed");
	return;

}

void RemoveDuplicatesThread( struct HydrogenBondIterator_s HBit )
{
	double MinLength;

	std::vector<struct HydrogenBond *>::iterator iter_hbmain;
	std::vector<struct HydrogenBond *>::iterator iter_hb;
	std::vector<struct HydrogenBond *>::iterator iter_hbmin;

	/*
	 * Look for acceptor duplicates
	 */

	for( iter_hbmain = HBit.begin; iter_hbmain != HBit.end-1; ++iter_hbmain )
	{
		// If this is already marked as a duplicate, skip it.
		if ( (*iter_hbmain)->markedDuplicate )
			continue;

		// Save the iterators for the minimum distance HBond atoms.
		// Set them to iter_?main initially.
		iter_hbmin = iter_hbmain;
		MinLength = (*iter_hbmin)->length;

		/*
		 * Go through entire vector looking for duplicate of
		 * iter_hbmain acceptor, and find the one with the shortest length
		 */
		for( iter_hb = iter_hbmain+1; iter_hb != HBit.end; ++iter_hb )
		{
			// If this is already marked as a duplicate, skip it.
			if ( (*iter_hb)->markedDuplicate )
				continue;

			if ( SameAtom( (*iter_hbmain)->acceptor, (*iter_hb)->acceptor) )
			{
				// duplicate = true;
				if ( (*iter_hb)->length < MinLength )
				{
					MinLength  = (*iter_hb)->length;
					// Mark the old iter_?min as duplicated;
					(*iter_hbmin)->markedDuplicate = true;
					iter_hbmin = iter_hb;
				}
				else
				{
					// Mark this one as duplicated;
					(*iter_hb)->markedDuplicate = true;
				}
			}
		}
	}

	//
	// Look for H duplicates
	//

	for(iter_hbmain = HBit.begin; iter_hbmain != HBit.end-1; ++iter_hbmain )
	{
		// If this is already marked as a duplicate, skip it.
		if ( (*iter_hbmain)->markedDuplicate )
			continue;

		// Save the iterator for the minimum distance HBond atoms.
		// Set them to iter_?main initially.
		iter_hbmin = iter_hbmain;
		MinLength = (*iter_hbmin)->length;

		// Go through entire vector looking for duplicate of
		// iter_Hmain, and find the one with the shortest length

		for( iter_hb = iter_hbmain+1 ; iter_hb != HBit.end; ++iter_hb )
		{
			// If this is already marked as a duplicate, skip it.
			if ( (*iter_hb)->markedDuplicate )
				continue;

			if ( SameAtom( (*iter_hbmain)->hydrogen, (*iter_hb)->hydrogen) )
			{
				if ( (*iter_hb)->length < MinLength )
				{
					MinLength = (*iter_hb)->length;
					// Mark the old iter_?min as duplicated;
					(*iter_hbmin)->markedDuplicate = true;
					iter_hbmin = iter_hb;
				}
				else
				{
					// This one isn't the minimum length.
					// Mark this one as duplicated;
					(*iter_hb)->markedDuplicate = true;
				}
			}
		}
	}

	return;
}

bool SameAtom( struct thbAtom *A,
               struct thbAtom *B)
{

	if ( A == B )
		return(true);


	return(false);
}

std::vector< std::vector<bool> >
Lifetime( std::vector< struct HydrogenBondIterator_s > *TrjIdx_iter,
          std::vector<struct HydrogenBond *>::iterator iter_hbmain,
          unsigned int NumFrames)
{
	// Look to see if this hydrogen bond exists in other frames.

	std::vector<struct HydrogenBond *>::iterator iter_hb;
	std::vector<struct HydrogenBond *>::iterator iter_begin;
	std::vector<struct HydrogenBond *>::iterator iter_end;

	iter_begin = TrjIdx_iter->at( 0 ).begin;
	iter_end   = TrjIdx_iter->at( 0 ).end;

	unsigned int NumHBsInFrameZero=0;

	for(iter_hb = iter_begin ; iter_hb < iter_end; ++iter_hb ) {
		NumHBsInFrameZero++; }
	std::cout << " Number of hydrogen bonds in initial frame: " << NumHBsInFrameZero << "\n";

	std::vector< std::vector<bool> >b(NumHBsInFrameZero, std::vector<bool>(NumFrames,false));

	if (NumHBsInFrameZero == 0)
		return b;

	for( unsigned int h=0; h < NumHBsInFrameZero; ++h)
	{
		b.at(h).at(0) = true;
		for( unsigned int f=1; f < NumFrames; ++f )
		{
			// The range of hydrogen bonds in frame f.
			iter_begin = TrjIdx_iter->at( f ).begin;
			iter_end   = TrjIdx_iter->at( f ).end;

			bool found = false;
			for(iter_hb = iter_begin ; iter_hb < iter_end; ++iter_hb )
			{
				if ( ((*iter_hb)->hydrogen == (*(iter_hbmain+h))->hydrogen) &&
				     ((*iter_hb)->acceptor == (*(iter_hbmain+h))->acceptor) )
				{
					// This matched the hydrogen bond in the first frame.
					found = true;
					break;
				}
			}
			if (found == true) {
				b.at(h).at(f) = true; }
		}
	}
	for( std::vector< std::vector<bool> >::iterator vbit=b.begin(); vbit < b.begin()+1; ++vbit)
	{
		for( std::vector<bool>::iterator bit=vbit->begin(); bit < vbit->end(); ++bit)
		{
			std::cout << (*bit==true?"|":" ");
		}
	std::cout << "\n";
	}

	return b;
}
/*
 * Return values:
 *
 * true   : A new string is found.
 * false  : Not a new string.
 *
 */
bool Trace( ListOfHBonds **HBonds,
            std::vector< struct HydrogenBondIterator_s > *TrjIdx_iter,
            std::vector<struct HydrogenBond *>::iterator iter_hbmain)
{
	// DonorO --- Hydrogen ... AcceptorO
	// ... Denotes the Hydrogen bond.

	std::vector<struct HydrogenBond *>::iterator iter_hb;

	std::vector<struct HydrogenBond *>::iterator iter_begin;
	std::vector<struct HydrogenBond *>::iterator iter_end;


	// If this hydrogen bond has already been assigned to a chain, skip it
	if ( ( (*iter_hbmain)->Next != NULL) || ( (*iter_hbmain)->Previous != NULL) )
		return(false);

	// The Range of HydrogenBonds to search.
	iter_begin = TrjIdx_iter->at( (*iter_hbmain)->TrajIdx ).begin;
	iter_end   = TrjIdx_iter->at( (*iter_hbmain)->TrajIdx ).end;

	// Starting a new chain.
	(*HBonds)->AddAtBegin(*iter_hbmain);

	bool StillLooking = true;
	while ( StillLooking )
	{
		bool FoundOne = false;
		for(iter_hb = iter_begin ; iter_hb < iter_end; ++iter_hb )
		{
			// If this hydrogen bond has already been assigned to a
			// chain, skip it
			if ( ( (*iter_hb)->Next != NULL) ||
			     ( (*iter_hb)->Previous != NULL) )
				continue;

			// Check that this hydrogen bond is not in the chain already.
			// Checking for only the H is sufficient.
			if ( !(*HBonds)->Find(*iter_hb)  )
			{
				if ( (*HBonds)->linksAtBegin(*iter_hb) )
				{
					// Found a new link at the beginning of the chain.
					(*HBonds)->AddAtBegin(*iter_hb);
					FoundOne = true;
				}
				else if ( (*HBonds)->linksAtEnd(*iter_hb) )
				{
					// Found a new link at the end of the chain.
					(*HBonds)->AddAtEnd(*iter_hb);
					FoundOne = true;
				}
			}
		}
		StillLooking = FoundOne;
	}

	return(true);
}

