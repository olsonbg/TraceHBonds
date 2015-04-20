#include "Print.h"
#include "OutputFormat.h"
#include "TraceHBonds.h"

#ifdef DEBUG
#define DEBUG_MSG(str) do { std::cout << "DEBUG: " << str << std::endl; } while ( false )
#else
#define DEBUG_MSG(str) do { } while ( false )
#endif

extern bool THB_VERBOSE;

// General parsing/formatting

std::string toString(int n)
{
	std::stringstream out(std::ios_base::out);
	out << n;
	return out.str();
}

// Real code.

int doAllFiles(char *progname,
               char *fPrefix , char *fSuffix, int first, int last,
               char *ofPrefix, char *ofSuffix,
               int NumBins, bool POVRAY)
{
	std::vector<struct HydrogenBond *> hb;
	std::vector<struct thbAtom *> atom;
	struct PBC *Cell;
	Cell = new struct PBC;
	ReadCarMdf( fPrefix, &atom, Cell );
	return(0);

	////////////////////////
	unsigned int filecounter=1;
	for (int fidx=first; fidx <= last; fidx++)
	{
		std::stringstream ifile;
		std::stringstream ofile;
		
		ifile << fPrefix  << fidx << fSuffix;
		// Send to stdout, '-', is ofPrefix and ofSuffix are not
		// specified.
		if ( (ofPrefix == NULL) && (ofSuffix == NULL) )
			ofile << "-";
		else
			ofile << ofPrefix << fidx << ofSuffix;

		if (ifile == NULL)
		{
			Help(progname);
			return(1);
		}

		/*
		 * If VERBOSE was requested, show number of file processed, updated
		 * every 50.
		 */
		if ( THB_VERBOSE &&
		     ((filecounter%50 == 0) || (filecounter == 1)) )
			std::cout << "Processing file " << filecounter << "/" 
			          << last-first+1 << ".\r" << std::flush;

		doFrame(ifile.str().c_str(), ofile.str().c_str(), NumBins, POVRAY);
		filecounter++;
	}

	if ( THB_VERBOSE)
			std::cout << "Processing file "
			          << last-first+1 << "/"
			          << last-first+1 << "." << std::endl;

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
			std::cout << ". ran out of memory?" << std::endl;
			return false;
		}
	}
	return true;
}

/*
 * Make sure x,y, and z coordinates of atom are each of size nelem, if not,
 * initialize the needed number of elements to val. Make sure not to
 * touch/change any values that are already in atom.
 */
bool alloc_vector(struct thbAtom *atom,
                  double val,
                  unsigned int nelem)
{
	if ( !alloc_vector( &(atom->x), val, nelem) ) return(false);
	if ( !alloc_vector( &(atom->y), val, nelem) ) return(false);
	if ( !alloc_vector( &(atom->z), val, nelem) ) return(false);

	return(true);
}

/*
 * Make sure the PBC parameters in cell are each of size nelem, if not,
 * initialize the needed number of elements to val. Make sure not to
 * touch/change any values that are already in cell.
 */
bool alloc_vector(struct PBC *cell,
                  double val,
                  unsigned int nelem)
{
	if ( !alloc_vector( &(cell->x), val, nelem) ) return(false);
	if ( !alloc_vector( &(cell->y), val, nelem) ) return(false);
	if ( !alloc_vector( &(cell->z), val, nelem) ) return(false);

	if ( !alloc_vector( &(cell->alpha), val, nelem) ) return(false);
	if ( !alloc_vector( &(cell->beta) , val, nelem) ) return(false);
	if ( !alloc_vector( &(cell->gamma), val, nelem) ) return(false);

	return(true);
}

/*
 * Make sure v is of size nelem*melem, if not, initialize the needed number of
 * elements to val. Make sure we don't adjust the values already stored in v.
 */
bool alloc_vector(vvui *v,
                  unsigned int val,
                  unsigned int nelem,
                  unsigned int melem)
{
	for ( unsigned int n=0; n < nelem; n++)
	{
		if ( n == v->size() )
		{
			try
			{
				vui Zero( nelem, 0);
				v->push_back(Zero);
			}
			catch( std::exception const &e)
			{
				std::cout << "exception: " << e.what();
				std::cout << ". ran out of memory?" << std::endl;
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

int doFrame(const char *ifile, const char *ofile,
            unsigned int NumBins, bool POVRAY)
{
	// Redirect to either a file, or std::cout.
	std::streambuf *buf;
	std::ofstream ofs;

	if ( !strncmp(ofile,"-",1) )
		buf = std::cout.rdbuf();
	else
	{
		ofs.open(ofile,std::ios::out);
		if ( !ofs.is_open() )
		{
			std::perror(ofile);
			return(1);
		}
		buf = ofs.rdbuf();
	}
	std::ostream out(buf);

	const char *CC1 = "#";
	const char *CC2 = "//";
	std::string CC;

	// Povray uses a different comment string.
	if ( POVRAY )
		CC = CC2;
	else
		CC = CC1;

	std::vector<struct HydrogenBond *> hb;
	std::vector<struct thbAtom *> atom;
	struct PBC *Cell;

	// Reserve space to prevent reallocation. If more than
	// 5000 hydrogen bonds, it will start reallocation.
	hb.reserve(40000);
	// Reserve space to prevent reallocation. If more than
	// 15000 atoms, it will start reallocation.
	atom.reserve(5000);

	Cell = new struct PBC;

	unsigned int NumFramesInTrajectory = 0;
	NumFramesInTrajectory = ReadData( ifile, &hb, &atom, Cell );

	std::vector< std::vector<struct HydrogenBond *>::iterator >TrjIdx_iter;
	TrjIdx_iter = TrajectoryIndexIterator( &hb );

	/*
	 * Show some initial information
	 */
	out << CC << "--- Before removing duplicates." << std::endl;
	out << CC << " Donor Oxygen atoms    : " << hb.size() << std::endl;
	out << CC << " Hydrogen atoms        : " << hb.size() << std::endl;
	out << CC << " Acceptor Oxygen atoms : " << hb.size() << std::endl;
	std::cout << CC << " Unique atoms          : " << atom.size() << std::endl;

	DEBUG_MSG("Removing duplicates: " << hb.size() << "Initially" );
	if ( THB_VERBOSE )
		std::cout << "Removing duplicates: "
		          <<hb.size()
		          << " Initially."
		          << std::endl;
	RemoveDuplicates ( &hb, &TrjIdx_iter );
	if ( THB_VERBOSE )
		std::cout << "Duplicates Removed: "
		          << hb.size()
		          << " Remaining."
		          << std::endl;

	// Update TrjIdx_iter after removing elements.
	TrjIdx_iter = TrajectoryIndexIterator( &hb );


	/*
	 * Show some more information
	 */
	out << CC << "--- Removed duplicates." << std::endl;
	out << CC << " Filename : " << ifile << std::endl;
	if ( NumBins != 0 )
		out <<CC<< " Minimum # of bins set to : " <<NumBins << std::endl;

	unsigned int TrjIdx = 0;

	out << CC
	    << " PBC "
	    << Cell->x.at(TrjIdx)     << " "
	    << Cell->y.at(TrjIdx)     << " "
	    << Cell->z.at(TrjIdx)     << " "
	    << Cell->alpha.at(TrjIdx) << " "
	    << Cell->beta.at(TrjIdx)  << " "
	    << Cell->gamma.at(TrjIdx)
	    << std::endl;

	out <<CC<< " Donor Oxygen atoms    : " << hb.size() << std::endl;
	out <<CC<< " Hydrogen atoms        : " << hb.size() << std::endl;
	out <<CC<< " Acceptor Oxygen atoms : " << hb.size() << std::endl;

	// Each element of the vector points to a string of hbonds.
	// ListOfHBonds is a strings of hbonds.
	std::vector<ListOfHBonds *>HBStrings;

	//Find all the strings.
	if ( THB_VERBOSE ) std::cout << "Tracing HB strings." << std::endl;

	for( unsigned int i=0; i < hb.size(); i++ )
	{
		ListOfHBonds *HBonds = new ListOfHBonds();
		if ( Trace( &HBonds, &TrjIdx_iter, hb.begin()+i) )
			HBStrings.push_back(HBonds);
		else
			delete HBonds;
	}

	if (THB_VERBOSE) std::cout << "Done tracing HB strings." << std::endl;

	if (THB_VERBOSE) std::cout << "Preparing histograms..." << std::endl;
	for( ; TrjIdx < NumFramesInTrajectory; ++TrjIdx )
	{
		if (THB_VERBOSE && ((TrjIdx+1)%50==0) )
			std::cout << "\tframe " << TrjIdx+1 << "\r" << std::flush;

		std::ofstream odata;
		std::string OutFile = "/dev/shm/t/__U_";
		OutFile += toString(TrjIdx+1);
		OutFile += ".zzz";

		odata.open(OutFile.c_str(),std::ios::out);
		if ( odata.is_open() )
		{
			//Make histograms, and printout the results.
			makeHistograms( &odata, HBStrings, CC, NumBins, Cell, TrjIdx, POVRAY);
			out << CC << "NEXT" <<std::endl;
			odata.close();
		}
	}
	if (THB_VERBOSE) std::cout << "Done with histograms." << std::endl;
	// Cleanup.
	DeleteVectorPointers(hb);
	DeleteVectorPointers(atom);
	hb.clear();
	atom.clear();

	delete Cell;

	for(unsigned int i=0; i < HBStrings.size();++i)
		delete HBStrings[i];
	HBStrings.clear();
	
	if ( ofs.is_open() )
		ofs.close();
	// Done with cleanup.

	return(0);
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
		*out << "#version 3.6;" << std::endl;
		*out << "global_settings {  assumed_gamma 1.0 }" << std::endl;
		*out << "Camera_LookAt( " << Cell->x.at(TrjIdx) << ", "
		                          << Cell->y.at(TrjIdx) << ", "
		                          << Cell->z.at(TrjIdx) << " )" << std::endl;
		*out << "PBC( " << Cell->x.at(TrjIdx) << ", "
		                << Cell->y.at(TrjIdx) << ", "
		                << Cell->z.at(TrjIdx) << " )" << std::endl;
	}

	// Printout information about each hbond string.
	for( unsigned int i=0; i < HBStrings.size(); i++ )
	{
		if ( HBStrings[i]->TrajectoryIndex() != TrjIdx )
			continue;

		*out << std::endl << std::endl;
		*out << CC << " Current Element : " << i << std::endl;
		*out << CC << " Atoms in Chain : " << HBStrings[i]->AtomCount();
		*out << std::endl;

		// Note if this is a closed loop.
		if ( HBStrings[i]->ClosedLoop() )
			*out << CC << " Closed Loop" << std::endl;

		*out << CC << " Molecules : "
		          << HBStrings[i]->MoleculeCount() << std::endl;

		*out << CC << " Unique forcefields : "
		          << HBStrings[i]->ForcefieldCount() << std::endl;

		*out << CC
		          << " Times chain switched between Molecules (switching) : "
		          << HBStrings[i]->SwitchingCount() << std::endl;

		*out << CC << " Periodic boundary conditions applied."
		          << std::endl;
		// Show the Chain atoms, molecules and coordinates
		EndToEndLength = HBStrings[i]->PrintAll(out, *Cell, TrjIdx, POVRAY);
		*out << CC << " Chain end-to-end distance: ";
		*out << colE2E << EndToEndLength << std::endl;
	}

	// Printout Histograms

	// Chain Length histogram table header
	unsigned int MaxBarLength = 62;
	*out << CC << std::endl << CC;
	*out << " Atoms/HBonds |Count| (For all Chains, including Closed Loops)";
	*out << std::endl;
	// Minimum Chain length is 3 atoms (O-H...O). Chain length is always an odd 
	// number of atoms, so use a step length of 2.
	PrintHistogramChain( out,
	                     hChainLength,
	                     MaxChainLength, 3,
	                     2, MaxBarLength,
	                     NumBins,CC);

	// Closed Loop histogram table header
	*out << CC <<std::endl << CC 
	          << " Atoms/HBonds |Count| (For Closed Loops)" << std::endl;
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
		*out << CC << std::endl << CC 
		          << " Switching |Count| (For Chain length of " << chainL << ")"
		          << std::endl;

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
		*out << CC << std::endl << CC 
		          << " Molecules |Count| (For Chain length of " << chainL << ")"
		          << std::endl;

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
