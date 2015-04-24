#include "Print.h"
#include "OutputFormat.h"
#include "TraceHBonds.h"
#include "WorkerThreads.h"

extern bool THB_VERBOSE;

#ifdef PTHREADS
extern pthread_mutex_t inQueueLock;
extern pthread_mutex_t outQueueLock;
extern std::list< unsigned int > inQueue;
extern std::list< unsigned int > outQueue;
extern std::vector< struct worker_data_s > worker_data;
extern std::vector< std::vector<struct HydrogenBond *> *> worker_hb;
#endif

const long double PI = 3.14159265358979323846;

inline
double Round (double r,double f=1.0)
{
	return (r > 0.0) ? floor(r*f + 0.5)/f : ceil(r*f - 0.5)/f;
}

// Vector pointing from atom A to atom B.
inline
std::vector<double> getAtomSeparationVector( std::vector<double> *A,
                                             std::vector<double> *B)
{
	std::vector<double> sep;

	sep.push_back( B->at(0) - A->at(0) );
	sep.push_back( B->at(1) - A->at(1) );
	sep.push_back( B->at(2) - A->at(2) );

	return(sep);
}

// Calculate the angle A-B-C in degrees
// For Hydrogren bond, this would be donor-hydrogen-acceptor
double getBondAngle( std::vector<double> A,
                     std::vector<double> B,
                     std::vector<double> C)
{
	std::vector<double> BA = getAtomSeparationVector(&B,&A);
	std::vector<double> BC = getAtomSeparationVector(&B,&C);

	double BAdotBC = BA[0]*BC[0] +
	                 BA[1]*BC[1] +
	                 BA[2]*BC[2];

	
	double magBA = sqrt( BA[0]*BA[0] + BA[1]*BA[1] + BA[2]*BA[2] );
	double magBC = sqrt( BC[0]*BC[0] + BC[1]*BC[1] + BC[2]*BC[2] );

	return ( acos(BAdotBC/magBA/magBC)*180.0/PI );
}

// Minimum Image Vector pointing from atom A to atom B.
// returns x,y,z as a vector of doubles.
inline
std::vector<double> getMinimumImageVector( std::vector<double> *A,
                                           std::vector<double> *B,
                                           std::vector<double> *Cell)
{
	std::vector<double> d;
	
	d = getAtomSeparationVector(A,B);

	double Lx = Cell->at(0);
	double Ly = Cell->at(1);
	double Lz = Cell->at(2);

	d[0] -= Round(d[0]/Lx)*Lx;
	d[1] -= Round(d[1]/Ly)*Ly;
	d[2] -= Round(d[2]/Lz)*Lz;

	return(d);
}

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
		  std::vector<double>cell,
		  std::vector<struct thbAtom *>*hydrogens,
		  std::vector<struct thbAtom *>*acceptors,
		  double TrjIdx, double rCutoff, double angleCutoff,
		  unsigned int ThreadID, unsigned int Threads)
{

	std::vector<struct thbAtom *>::iterator it_h;
	std::vector<struct thbAtom *>::iterator it_a;

	std::vector<double> r;

	double rCutoff2 = pow(rCutoff,2.0);
	double r2;

	for( it_h = hydrogens->begin() + ThreadID;
	     it_h < hydrogens->end();
	     it_h += Threads)
	{
		std::vector<double>a;
		a.push_back( (*it_h)->x.at(TrjIdx) );
		a.push_back( (*it_h)->y.at(TrjIdx) );
		a.push_back( (*it_h)->z.at(TrjIdx) );

		for( it_a = acceptors->begin(); it_a < acceptors->end(); ++it_a)
		{
			// Make sure this acceptor is not covalently bonded to the hydrogen
			// we are looking at.  The angle check would catch this, but this
			// will skip a few calculations.
			if ( (*it_a) == (*it_h)->ConnectedAtom.at(0) )
				continue;

			std::vector<double>b;
			b.push_back( (*it_a)->x.at(TrjIdx) );
			b.push_back( (*it_a)->y.at(TrjIdx) );
			b.push_back( (*it_a)->z.at(TrjIdx) );
			r = getMinimumImageVector( &a, &b, &cell );
			r2 = pow(r[0],2.0) + pow(r[1],2.0) + pow(r[2],2.0);

			if ( r2 < rCutoff2)
			{
				// Distance cutoff is good, now check the angle.
				std::vector<double>c;
				c.push_back( (*it_h)->ConnectedAtom.at(0)->x.at(TrjIdx));
				c.push_back( (*it_h)->ConnectedAtom.at(0)->y.at(TrjIdx));
				c.push_back( (*it_h)->ConnectedAtom.at(0)->z.at(TrjIdx));
				double angle = getBondAngle(c,a,b);
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

#ifdef PTHREADS
	// How many jobs to split this work into; The number of jobs to put
	// in the queue. Use the same as the number of threads created.
	unsigned int NUM_JOBS=NumWorkerThreads();

	// Initialize the variables where the worker threads put their
	// results.
	for (unsigned int j=0; j < NUM_JOBS; ++j)
	{
		std::vector<struct HydrogenBond *> *thread_hb;
		thread_hb = new std::vector<struct HydrogenBond *>;
		worker_hb.push_back(thread_hb);
	}
	// Tell the threads to start looking at the inQueue for work to do.
	ContinueWorkerThreads();
#endif

	for( unsigned int TrjIdx=0; TrjIdx < Cell->frames; ++TrjIdx)
	{
		VERBOSE_RMSG("Processing frame " << TrjIdx+1 <<"/"<< Cell->frames << ". Hydrogen-acceptor pairs found: " << hb->size() << ".");
		

		std::vector<double>cell;
		cell.push_back( Cell->x.at(TrjIdx) );
		cell.push_back( Cell->y.at(TrjIdx) );
		cell.push_back( Cell->z.at(TrjIdx) );

#ifdef PTHREADS
		for (unsigned int j=0; j < NUM_JOBS; ++j)
		{
			// setup the job data.
			pthread_mutex_lock(&inQueueLock);

			worker_data.at(j).jobtype = THREAD_JOB_HBS;
			worker_data.at(j).jobnum = j;
			worker_data.at(j).num_threads = NUM_JOBS;
			worker_data.at(j).hb = worker_hb[j];
			worker_data.at(j).cell = cell;
			worker_data.at(j).hydrogens = &hydrogens;
			worker_data.at(j).acceptors = &acceptors;
			worker_data.at(j).TrjIdx = TrjIdx;
			worker_data.at(j).rCutoff = rCutoff;
			worker_data.at(j).angleCutoff = angleCutoff;

			inQueue.push_back(j);

			pthread_mutex_unlock(&inQueueLock);
		}
		// Wait for out queue to fill with all jobs.
		unsigned int Done_count = 0;
		while (Done_count != NUM_JOBS )
		{
			pthread_mutex_lock(&outQueueLock);
			Done_count = outQueue.size();
			pthread_mutex_unlock(&outQueueLock);
		}
		outQueue.clear();
		// Concatenate the worker_hb to the real hb, then clear the worker_hb
		// vector.
		for(unsigned int j=0; j < NUM_JOBS; ++j)
		{
			hb->reserve( hb->size() + worker_data.at(j).hb->size() );
			hb->insert(hb->end(),
			           worker_data.at(j).hb->begin(),
			           worker_data.at(j).hb->end() );

			worker_data.at(j).hb->clear();
		}
		
#else
		HBs( hb, cell, &hydrogens, &acceptors, TrjIdx, rCutoff, angleCutoff);
#endif
	}

#ifdef PTHREADS
	// Do not need the threads for now, so tell them to pause.
	PauseWorkerThreads();
	// Cleanup
	DeleteVectorPointers(worker_hb);
#endif


	VERBOSE_MSG("");
}

int doArcFile(char *ifilename,
              char *ofPrefix, char *ofSuffix,
              struct HydrogenBondMatching *match,
              double rCutoff, double angleCutoff,
              int NumBins, bool POVRAY)
{
	std::vector<struct HydrogenBond *> hb;
	std::vector<struct thbAtom *> atom;

	struct PBC *Cell;
	Cell = new struct PBC;

	ReadCarMdf( ifilename, &atom, Cell );
	std::vector<double>A, B, C;

	unsigned int NumFramesInTrajectory = 0;
	NumFramesInTrajectory = Cell->frames;

	VERBOSE_MSG("Total frames: " << NumFramesInTrajectory);

	// Now  determine the hydrogen bonds
	AtomNeighbors( &hb, &atom, Cell, match, rCutoff, angleCutoff );


	std::vector< std::vector<struct HydrogenBond *>::iterator >TrjIdx_iter;
	TrjIdx_iter = TrajectoryIndexIterator( &hb );

	VERBOSE_MSG("Looking for smallest hydrogen-acceptor bond lengths in all frames...");

	RemoveDuplicates ( &hb, &TrjIdx_iter );

	VERBOSE_MSG("Hydrogen bonds:          " << hb.size() << ".");

	// Update TrjIdx_iter after removing elements.
	TrjIdx_iter = TrajectoryIndexIterator( &hb );

	unsigned int TrjIdx;

	// Each element of the vector points to a string of hbonds.
	// ListOfHBonds is a strings of hbonds.
	std::vector<ListOfHBonds *>HBStrings;

	//Find all the strings.
	VERBOSE_MSG("Tracing HB strings.");

	for( unsigned int i=0; i < hb.size(); i++ )
	{
		ListOfHBonds *HBonds = new ListOfHBonds();
		if ( Trace( &HBonds, &TrjIdx_iter, hb.begin()+i) )
			HBStrings.push_back(HBonds);
		else
			delete HBonds;
	}

	VERBOSE_MSG("Done tracing HB strings.");


	const char *CC1 = "#";
	const char *CC2 = "//";
	std::string CC;

	// Povray uses a different comment string.
	if ( POVRAY )
		CC = CC2;
	else
		CC = CC1;
	for( TrjIdx = 0 ; TrjIdx < NumFramesInTrajectory; ++TrjIdx )
	{
		if (  ((TrjIdx+1)%50==0) || ((TrjIdx+1)==NumFramesInTrajectory)  )
			VERBOSE_RMSG("Preparing histograms, frame " << TrjIdx+1 << "/" << NumFramesInTrajectory);

		std::stringstream ofilename;
		ofilename << ofPrefix << TrjIdx+1 << ofSuffix;

		std::ofstream out;
		out.open(ofilename.str().c_str(),std::ios::out);
		if ( out.is_open() )
		{
			out << CC
				<< " PBC "
				<< Cell->x.at(TrjIdx)     << " "
				<< Cell->y.at(TrjIdx)     << " "
				<< Cell->z.at(TrjIdx)     << " "
				<< Cell->alpha.at(TrjIdx) << " "
				<< Cell->beta.at(TrjIdx)  << " "
				<< Cell->gamma.at(TrjIdx)
				<< "\n";

			out <<CC<< " Donor Oxygen atoms    : " << hb.size() << "\n";
			out <<CC<< " Hydrogen atoms        : " << hb.size() << "\n";
			out <<CC<< " Acceptor Oxygen atoms : " << hb.size() << "\n";

			// Make histograms, and printout the results.
			makeHistograms( &out, HBStrings, CC, NumBins, Cell, TrjIdx, POVRAY);
			out.close();
		}
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

/*
 * Make sure v is of size nelem, if not, initialize the needed number of
 * elements to val. Make sure not to touch/change any values that are already
 * in v.
 */
template<class T> bool alloc_vector( std::vector<T> *v,
                                     T val,
                                     unsigned int nelem)
{
	for ( unsigned int n=v->size(); n < nelem; n++)
	{
		try
		{
			v->push_back(val);
		}
		catch( std::exception const &e)
		{
			std::cout << "exception: " << e.what();
			std::cout << ". ran out of memory?" << "\n";
			return false;
		}
	}
	return true;
}

/*
 * Make sure v is of size nelem*melem, if not, initialize the needed number of
 * elements to val. Make sure we don't adjust the values already stored in v.
 */
template<class T> bool alloc_vector( std::vector< std::vector<T> > *v,
                                     T val,
                                     unsigned int nelem,
                                     unsigned int melem)
{
	for ( unsigned int n=0; n < nelem; n++)
	{
		if ( n == v->size() )
		{
			try
			{
				std::vector<T> Zero( nelem, val);
				v->push_back(Zero);
			}
			catch( std::exception const &e)
			{
				std::cout << "exception: " << e.what();
				std::cout << ". ran out of memory?" << "\n";
				return false;
			}
		}
		else
		{
			bool ret = alloc_vector(&(v->at(n)), val, melem);

			if ( ret == false )
				return false;
		}
	}
	return true;
}

/*
 * Add to the counts in h[bin], making sure enough space is allocated.
 * Also record the maximum bin. max may already be assigned a value.
*/
bool Bin(vui *h, unsigned int *max, unsigned int bin)
{
	if ( alloc_vector(h, 0U, bin+1) )
		h->at(bin)++;
	else
		return 1;

	if ( bin > *max)
		*max = bin;

	return true;
}

/*
 * Add to the counts in the h[bini_i][bin_j] h, making sure enough space is
 * allocated. Also record the maximum bin_j in hmax, making sure enough space
 * is allocated. hmax may already be assigned a value.
*/
bool Bin(vvui *h, vui *hmax, unsigned int hb, unsigned int c)
{

	if( alloc_vector(h, 0U, hb+1, c+1))
		(h->at(hb)).at(c)++;
	else
		return false;

	if( !alloc_vector(hmax, 0U, hb+1) )
		return 1;

	if (c > hmax->at(hb))
		hmax->at(hb) = c;
	return true;
}

// Save Iterators which point to just past the end of a Trajectory Index.
// TrjIdx_iter.at(1) points to first element of TrjIdx 1. TrjIdx_iter.at(2)
// points to just past the last element of TrjIdx 1, or the first element
// of TrjIdx 2.
std::vector< std::vector<struct HydrogenBond *>::iterator >
TrajectoryIndexIterator( std::vector<struct HydrogenBond *> *hb)
{
	std::vector< std::vector<struct HydrogenBond *>::iterator >TrjIdx_iter;

	TrjIdx_iter.push_back( hb->begin() );

	unsigned int counter=1;
	std::vector<struct HydrogenBond *>::iterator it_hb;
	for(it_hb = hb->begin(); it_hb < hb->end(); ++it_hb)
	{
		if ( (*it_hb)->TrajIdx == counter )
		{
			TrjIdx_iter.push_back( it_hb );
			counter++;
		}
	}
	TrjIdx_iter.push_back( hb->end() );

	return(TrjIdx_iter);
}

template<class T> void DeleteVectorPointers( T v )
{
	for(unsigned int i =0; i < v.size(); ++i)
		delete v[i];
}

void removeMarked( std::vector<struct HydrogenBond *> *hb )
{
	std::vector<struct HydrogenBond *>::iterator iter_hb = hb->begin();

	for( ; iter_hb < hb->end(); )
	{
		if ( (*iter_hb)->markedDuplicate )
		{
			delete *iter_hb;
			iter_hb  = hb->erase(iter_hb);
		}
		else
		{
			++iter_hb;
		}
	}
}

/*
 * Possible to have more than one:
 *
 *   - acceptor Oxygen (aO) with a single Hydrogen (H)
 *   - Hydrogen (H) with a single acceptor Oxygen (aO)
 *
 * Find the all duplicates and keep the shortest Hydrogen Bond length.
 */

void RemoveDuplicates( std::vector<struct HydrogenBond *> *hb,
                       std::vector< std::vector<struct HydrogenBond *>::iterator > *TrjIdx_iter)
{
	double MinLength;

	std::vector<struct HydrogenBond *>::iterator iter_hbmain;
	std::vector<struct HydrogenBond *>::iterator iter_hb;
	std::vector<struct HydrogenBond *>::iterator iter_hbmin;

	std::vector<struct HydrogenBond *>::iterator iter_begin;
	std::vector<struct HydrogenBond *>::iterator iter_end;

	if ( hb->size() == 0 )
		return;

	/*
	 * Look for acceptor duplicates
	 */

	for( iter_hbmain = hb->begin(); iter_hbmain < hb->end()-1; ++iter_hbmain )
	{
		// If this is already marked as a duplicate, skip it.
		if ( (*iter_hbmain)->markedDuplicate )
			continue;

		// Save the iterators for the minimum distance HBond atoms.
		// Set them to iter_?main initially.
		iter_hbmin = iter_hbmain;
		MinLength = (*iter_hbmin)->length;

		// The Range of HydrogenBonds to search.
		iter_end   = TrjIdx_iter->at( (*iter_hbmain)->TrajIdx+1 );

		/*
		 * Go through entire vector looking for duplicate of
		 * iter_hbmain acceptor, and find the one with the shortest length
		 */
		for( iter_hb = iter_hbmain+1; iter_hb < iter_end; ++iter_hb )
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

	for(iter_hbmain = hb->begin(); iter_hbmain < hb->end()-1; ++iter_hbmain )
	{
		// If this is already marked as a duplicate, skip it.
		if ( (*iter_hbmain)->markedDuplicate )
			continue;

		// Save the iterator for the minimum distance HBond atoms.
		// Set them to iter_?main initially.
		iter_hbmin = iter_hbmain;
		MinLength = (*iter_hbmin)->length;

		// The Range of HydrogenBonds to search.
		iter_end   = TrjIdx_iter->at( (*iter_hbmain)->TrajIdx+1 );

		// Go through entire vector looking for duplicate of
		// iter_Hmain, and find the one with the shortest length

		for( iter_hb = iter_hbmain+1 ; iter_hb < iter_end; ++iter_hb )
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

	removeMarked(hb);
	return;
}

bool SameAtom( struct thbAtom *A,
               struct thbAtom *B)
{

	if ( A == B )
		return(true);


	return(false);
}

/*
 * Return values:
 *
 * true   : A new string is found.
 * false  : Not a new string.
 *
 */
bool Trace( ListOfHBonds **HBonds,
            std::vector< std::vector<struct HydrogenBond *>::iterator > *TrjIdx_iter,
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
	iter_begin = TrjIdx_iter->at( (*iter_hbmain)->TrajIdx );
	iter_end   = TrjIdx_iter->at( (*iter_hbmain)->TrajIdx+1 );

	// Starting a new chain.
	(*HBonds)->AddAtBegin(*iter_hbmain);

	bool StillLooking = true;
	while ( StillLooking )
	{
		bool FoundOne = false;
		for(iter_hb = iter_begin ; iter_hb < iter_end; ++iter_hb )
		{
			// If this hydrogen bond has already been fully assigned to a
			// chain, skip it
			if ( ( (*iter_hb)->Next != NULL) &&
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

int makeHistograms( std::ostream *out,
                     std::vector<ListOfHBonds *> HBStrings,
                     std::string CC, unsigned int NumBins,
                     struct PBC *Cell, unsigned int TrjIdx,
                     bool POVRAY)
{
	unsigned int MaxChainLength = 0;
	unsigned int MaxLoopSize = 0;
	double EndToEndLength;

	/*
	 * Go through the vector of HBond strings and:
	 *  Bin Chain Lengths                             (1D)
	 *  Bin Chain Lengths of only Closed Loops        (1D)
	 *  Bin Molecule Switches for each Chain Length   (2D)
	 *  Bin Molecules in Chain, for each Chain Length (2D)
	 */

	// Zero all histogram bins. Set 20 elements initially.
	vui hChainLength(20,0);
	vui hClosedLoop(20,0);
	
	vvui hSwitchesInChain( 20, vui (20,0));
	vvui hMoleculesInChain( 20, vui (20,0));

	vui MaxSwitchesInChain(20,0);
	vui MaxMoleculesInChain(20,0);
	// All histograms zeroed.

	for( unsigned int i=0; i < HBStrings.size(); i++ )
	{
		if ( HBStrings[i]->TrajectoryIndex() != TrjIdx )
			continue;

		unsigned int HBCount        = HBStrings[i]->AtomCount();
		unsigned int SwitchingCount = HBStrings[i]->SwitchingCount();
		unsigned int MoleculeCount  = HBStrings[i]->MoleculeCount();

		// Bin the chain lengths.
		if( !Bin(&hChainLength, &MaxChainLength, HBCount) )
			return 1;

		// Bin the chain lengths for only closed loops.
		if ( HBStrings[i]->ClosedLoop() )
		{
			if( !Bin(&hClosedLoop, &MaxLoopSize, HBCount) )
				return 1;
		}

		// Bin the number of molecule switches for each chain length.
		if( !Bin(&hSwitchesInChain, &MaxSwitchesInChain,
		         HBCount, SwitchingCount) )
			return 1;

		// Tabulate the number of molecules in each chain length.
		if ( !Bin(&hMoleculesInChain,&MaxMoleculesInChain,
		          HBCount,MoleculeCount) )
			return 1;
	}

	OFmt colE2E(0,6);

	// Header for povray file.
	if (POVRAY)
	{
		*out << "#version 3.6;" << "\n";
		*out << "global_settings {  assumed_gamma 1.0 }" << "\n";
		*out << "Camera_LookAt( " << Cell->x.at(TrjIdx) << ", "
		                          << Cell->y.at(TrjIdx) << ", "
		                          << Cell->z.at(TrjIdx) << " )" << "\n";
		*out << "PBC( " << Cell->x.at(TrjIdx) << ", "
		                << Cell->y.at(TrjIdx) << ", "
		                << Cell->z.at(TrjIdx) << " )" << "\n";
	}

	// Printout information about each hbond string.
	for( unsigned int i=0; i < HBStrings.size(); i++ )
	{
		if ( HBStrings[i]->TrajectoryIndex() != TrjIdx )
			continue;

		*out << "\n\n";
		*out << CC << " Current Element : " << i << "\n";
		*out << CC << " Atoms in Chain : " << HBStrings[i]->AtomCount();
		*out << "\n";

		// Note if this is a closed loop.
		if ( HBStrings[i]->ClosedLoop() )
			*out << CC << " Closed Loop" << "\n";

		*out << CC << " Molecules : "
		          << HBStrings[i]->MoleculeCount() << "\n";

		*out << CC << " Unique forcefields : "
		          << HBStrings[i]->ForcefieldCount() << "\n";

		*out << CC
		          << " Times chain switched between Molecules (switching) : "
		          << HBStrings[i]->SwitchingCount() << "\n";

		*out << CC << " Periodic boundary conditions applied."
		          << "\n";
		// Show the Chain atoms, molecules and coordinates
		EndToEndLength = HBStrings[i]->PrintAll(out, *Cell, TrjIdx, POVRAY);
		*out << CC << " Chain end-to-end distance: ";
		*out << colE2E << EndToEndLength << "\n";
	}

	// Printout Histograms

	// Chain Length histogram table header
	unsigned int MaxBarLength = 62;
	*out << CC << "\n" << CC;
	*out << " Atoms/HBonds |Count| (For all Chains, including Closed Loops)";
	*out << "\n";
	// Minimum Chain length is 3 atoms (O-H...O). Chain length is always an odd 
	// number of atoms, so use a step length of 2.
	PrintHistogramChain( out,
	                     hChainLength,
	                     MaxChainLength, 3,
	                     2, MaxBarLength,
	                     NumBins,CC);

	// Closed Loop histogram table header
	*out << CC <<"\n" << CC 
	          << " Atoms/HBonds |Count| (For Closed Loops)" << "\n";
	// Minimum Chain length is 3 atoms (O-H...O). Chain length is always an odd 
	// number of atoms, so use a step length of 2.
	PrintHistogramChain( out,
	                     hClosedLoop,
	                     MaxLoopSize, 3,
	                     2, MaxBarLength,
	                     NumBins,CC);

	// Molecules in each Chain histogram

	// Minimum Chain length is 3 atoms (O-H..O). Chain length is always an odd 
	// number of atoms, so increase counter by 2 for each step.
	unsigned int MaxChainL;
	if ( NumBins > 0 )
		MaxChainL = NumBins;
	else
		MaxChainL = MaxChainLength;

	for(unsigned int chainL=3; chainL <= MaxChainL; chainL += 2)
	{
		// Number of Switches histogram table header
		*out << CC << "\n" << CC 
		          << " Switching |Count| (For Chain length of " << chainL << ")"
		          << "\n";

		if ( chainL > MaxChainLength )
		{
			vui dummy;
			PrintHistogramMolecules( out,
			                         dummy,
			                         0, 0,
			                         1, MaxBarLength,
			                         NumBins, CC);
		}
		else
			PrintHistogramMolecules( out,
									 hSwitchesInChain[chainL],
									 MaxSwitchesInChain[chainL],0,
									 1, MaxBarLength,
									 NumBins, CC);
	}

	// Minimum Chain length is 3 atoms (O-H..O). Chain length is always an odd 
	// number of atoms, so increase counter by 2 for each step.

	for(unsigned int chainL=3; chainL <= MaxChainL; chainL += 2)
	{
		// Number of Molecules histogram table header
		*out << CC << "\n" << CC 
		          << " Molecules |Count| (For Chain length of " << chainL << ")"
		          << "\n";

		if ( chainL > MaxChainLength)
		{
			vui dummy;
			PrintHistogramMolecules( out,
									 dummy,
									 0, 1,
									 1, MaxBarLength,
									 NumBins, CC);
		}
		else
			PrintHistogramMolecules( out,
									 hMoleculesInChain[chainL],
									 MaxMoleculesInChain[chainL], 1,
									 1, MaxBarLength,
									 NumBins, CC);
	}
	return(0);
}
