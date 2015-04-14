#include "Print.h"
#include "OutputFormat.h"

extern bool thbVERBOSE;

int doAllFiles(char *progname,
               char *fPrefix , char *fSuffix, int first, int last,
               char *ofPrefix, char *ofSuffix,
               int NumBins, bool POVRAY)
{
	unsigned int filecounter=1;
	for (int fidx=first; fidx <= last; fidx++)
	{
		std::stringstream ifile;
		std::stringstream ofile;
		
		ifile << fPrefix  << fidx << fSuffix;
		if ( (ofPrefix == NULL) && (ofSuffix == NULL) )
			ofile << "-";
		else
			ofile << ofPrefix << fidx << ofSuffix;

		if (ifile == NULL)
		{
			Help(progname);
			return(1);
		}

		if ( thbVERBOSE &&
		     ((filecounter%50 == 0) || (filecounter == 1)) )
			std::cout << "Processing file " << filecounter << "/" << last-first+1 << ".\r" << std::flush;

		doFrame(ifile.str().c_str(), ofile.str().c_str(), NumBins, POVRAY);
		filecounter++;
	}
	if ( thbVERBOSE)
			std::cout << "Processing file " << last-first+1 << "/" << last-first+1 << "." << std::endl;;

	return(0);
}
// Make sure v is of size nelem, if not, initialize the needed number of 
// elements to val.
bool alloc_vector(std::vector<unsigned int> *v,
                  unsigned int val,
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

// Make sure v is of size nelem*melen, if not, initialize the needed number of 
// elements to val.
bool alloc_vector(std::vector< std::vector<unsigned int> >*v,
                  unsigned int val,
                  unsigned int nelem,
                  unsigned int melem)
{
	for ( unsigned int n=0; n < nelem; n++)
	{
		if ( n == v->size() )
		{
			std::vector<unsigned int>Zero( nelem, 0);
			v->push_back(Zero);
		}
		else
		{
			for ( unsigned int m=v[n].size(); m < melem; m++)
			{
				try
				{
					v->at(n).push_back(val);
				}
				catch( std::exception const &e)
				{
					std::cout << "exception: " << e.what();
					std::cout << ". ran out of memory?" << std::endl;
					return false;
				}
			}
		}
	}
	return true;
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

	std::vector<struct HBondAtom *> vec_H;
	std::vector<struct HBondAtom *> vec_aO;
	std::vector<struct HBondAtom *> vec_dO;
	struct PBC *Cell;

	// Reserve space to prevent reallocation. If more than
	// 5000 atoms, it will start reallocation.
	vec_H.reserve(5000);
	vec_aO.reserve(5000);
	vec_dO.reserve(5000);

	Cell = new struct PBC;

	ReadData( ifile, &vec_H, &vec_aO, &vec_dO, Cell );

	/*
	 * Show some initial information
	 */
	out <<CC<< "--- Before removing duplicates." << std::endl;
	out <<CC<< " Donor Oxygen atoms    : " << vec_dO.size() << std::endl;
	out <<CC<< " Hydrogen atoms        : " << vec_H.size() << std::endl;
	out <<CC<< " Acceptor Oxygen atoms : " << vec_aO.size() << std::endl;

	RemoveDuplicates ( &vec_H, &vec_aO, &vec_dO );

	/*
	 * Show some more information
	 */
	out <<CC<< "--- Removed duplicates." << std::endl;
	out <<CC<< " Filename : " << ifile << std::endl;
	if ( NumBins != 0 )
		out <<CC<< " Minimum # of bins set to : " <<NumBins << std::endl;
	out <<CC<< " PBC " << Cell->x << " " << Cell->y << " " << Cell->z 
	          << " " << Cell->alpha << " " << Cell->beta << " " << Cell->gamma
	          << std::endl;
	out <<CC<< " Donor Oxygen atoms    : " << vec_dO.size() << std::endl;
	out <<CC<< " Hydrogen atoms        : " << vec_H.size() << std::endl;
	out <<CC<< " Acceptor Oxygen atoms : " << vec_aO.size() << std::endl;

	if (POVRAY)
	{
		out << "#version 3.6;" << std::endl;
		out << "global_settings {  assumed_gamma 1.0 }" << std::endl;
		out << "Camera_LookAt( " << Cell->x << ", "
		                               << Cell->y << ", "
		                               << Cell->z << " )" << std::endl;
		out << "PBC( " << Cell->x << ", "
		                     << Cell->y << ", "
		                     << Cell->z << " )" << std::endl;
	}

	// Each element of the vector points to a string of hbonds.
	// ListOfHBonds is a strings of hbonds.
	std::vector<ListOfHBonds *>HBStrings;

	//Find all the strings.
	for( unsigned int i=0; i < vec_H.size(); i++ )
	{
		ListOfHBonds *HBonds = new ListOfHBonds();
		if ( Trace( &HBonds, vec_H, vec_aO, vec_dO, i) )
			HBStrings.push_back(HBonds);
		else
			delete HBonds;
	}

	//Make histograms, and printout the results.
	makeHistograms( &out, HBStrings, CC, NumBins, POVRAY, Cell);

	// Cleanup.
	DeleteVectorPointers(vec_H);
	DeleteVectorPointers(vec_aO);
	DeleteVectorPointers(vec_dO);
	vec_H.clear();
	vec_aO.clear();
	vec_dO.clear();

	delete Cell;

	for(unsigned int i=0; i < HBStrings.size();++i)
		delete HBStrings[i];
	HBStrings.clear();
	
	if ( ofs.is_open() )
		ofs.close();
	// Done with cleanup.

	return(0);
}

void DeleteVectorPointers( std::vector<struct HBondAtom *> v)
{
	for(unsigned int i =0; i < v.size(); ++i)
		delete v[i];
}

void removeMarked( std::vector<struct HBondAtom *> *H,
                   std::vector<struct HBondAtom *> *aO,
                   std::vector<struct HBondAtom *> *dO )
{
	std::vector<struct HBondAtom *>::iterator iter_H  = H->begin();
	std::vector<struct HBondAtom *>::iterator iter_aO = aO->begin();
	std::vector<struct HBondAtom *>::iterator iter_dO = dO->begin();

	for( ; iter_H < H->end(); )
	{
		if ( (*iter_H)->markedDuplicate )
		{
			delete *iter_H;
			delete *iter_aO;
			delete *iter_dO;
			iter_H  = H->erase(iter_H);
			iter_aO = aO->erase(iter_aO);
			iter_dO = dO->erase(iter_dO);
		}
		else
		{
			++iter_H;
			++iter_aO;
			++iter_dO;
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
 *
 * aO->length is the distance between H...aO (HBond length)
 */

void RemoveDuplicates( std::vector<struct HBondAtom *> *H,
                       std::vector<struct HBondAtom *> *aO,
                       std::vector<struct HBondAtom *> *dO )
{
	/*
	 * H  : Hydrogen
	 * aO : acceptor Oxygen
	 * dO : donor Oxygen
	 *
	 * DonorO --- Hydrogen ... AcceptorO
	 * ... Denotes the Hydrogen bond.
	 * 
	 * aO->length : Hydrogen - Acceptor Length (HBond Length)
	 */

	double MinLength;
	// bool duplicate;
	std::vector<struct HBondAtom *>::iterator iter_Hmain;
	std::vector<struct HBondAtom *>::iterator iter_aOmain;
	std::vector<struct HBondAtom *>::iterator iter_dOmain;

	std::vector<struct HBondAtom *>::iterator iter_H;
	std::vector<struct HBondAtom *>::iterator iter_aO;
	std::vector<struct HBondAtom *>::iterator iter_dO;

	std::vector<struct HBondAtom *>::iterator iter_Hmin;
	std::vector<struct HBondAtom *>::iterator iter_aOmin;
	std::vector<struct HBondAtom *>::iterator iter_dOmin;
	/*
	 * Look for aO duplicates
	 */
	iter_Hmain = H->begin();
	iter_aOmain = aO->begin();
	iter_dOmain = dO->begin();

	for( ; iter_Hmain < H->end()-1; ++iter_Hmain, ++iter_aOmain, ++iter_dOmain )
	{
		// If this is already marked as a duplicate, skip it.
		if ( (*iter_Hmain)->markedDuplicate )
			continue;

		// Save the iterators for the minimum distance HBond atoms.
		// Set them to iter_?main initially.
		iter_Hmin  = iter_Hmain;
		iter_aOmin = iter_aOmain;
		iter_dOmin = iter_dOmain;
		MinLength = (*iter_aOmin)->length;

		/*
		 * Go through entire vector looking for duplicate of
		 * iter_aOmain, and find the one with the shortest length
		 */
		// duplicate = false;
		iter_H = iter_Hmain+1;
		iter_aO = iter_aOmain+1;
		iter_dO = iter_dOmain+1;
		for( ; iter_H < H->end(); ++iter_H, ++iter_aO, ++iter_dO )
		{
			// If this is already marked as a duplicate, skip it.
			if ( (*iter_H)->markedDuplicate )
				continue;

			if ( SameAtom( *iter_aOmain, *iter_aO) )
			{
				// duplicate = true;
				if ( (*iter_aO)->length < MinLength )
				{
					MinLength  = (*iter_aO)->length;
					// Mark the old iter_?min as duplicated;
					(*iter_Hmin)->markedDuplicate = true;
					(*iter_aOmin)->markedDuplicate = true;
					(*iter_dOmin)->markedDuplicate = true;
					// Update the iterators for the minimum distance
					// HBond atoms.
					iter_Hmin  = iter_H;
					iter_aOmin = iter_aO;
					iter_dOmin = iter_dO;
				}
				else
				{
					// Mark this one as duplicated;
					(*iter_H)->markedDuplicate = true;
					(*iter_aO)->markedDuplicate = true;
					(*iter_dO)->markedDuplicate = true;
				}
			}
		}
	}

	//
	// Look for H duplicates
	//
	iter_Hmain = H->begin();
	iter_aOmain = aO->begin();
	iter_dOmain = dO->begin();

	for( ; iter_Hmain < H->end()-1; ++iter_Hmain, ++iter_aOmain, ++iter_dOmain )
	{
		// If this is already marked as a duplicate, skip it.
		if ( (*iter_Hmain)->markedDuplicate )
			continue;

		// Save the iterators for the minimum distance HBond atoms.
		// Set them to iter_?main initially.
		iter_Hmin  = iter_Hmain;
		iter_aOmin = iter_aOmain;
		iter_dOmin = iter_dOmain;
		MinLength = (*iter_aOmin)->length;

		/*
		 * Go through entire vector looking for duplicate of
		 * iter_Hmain, and find the one with the shortest length
		 */
		// duplicate = false;
		iter_H = iter_Hmain+1;
		iter_aO = iter_aOmain+1;
		iter_dO = iter_dOmain+1;
		for( ; iter_H < H->end(); ++iter_H, ++iter_aO, ++iter_dO )
		{
			// If this is already marked as a duplicate, skip it.
			if ( (*iter_H)->markedDuplicate )
				continue;

			if ( SameAtom( *iter_Hmain, *iter_H) )
			{
				// duplicate = true;
				if ( (*iter_aO)->length < MinLength )
				{
					MinLength = (*iter_aO)->length;
					// Mark the old iter_?min as duplicated;
					(*iter_Hmin)->markedDuplicate = true;
					(*iter_aOmin)->markedDuplicate = true;
					(*iter_dOmin)->markedDuplicate = true;
					// Update the iterators for the minimum distance
					// HBond atoms.
					iter_Hmin  = iter_H;
					iter_aOmin = iter_aO;
					iter_dOmin = iter_dO;
				}
				else
				{
					// This one isn't the minimum length.
					// Mark this one as duplicated;
					(*iter_H)->markedDuplicate = true;
					(*iter_aO)->markedDuplicate = true;
					(*iter_dO)->markedDuplicate = true;
				}
			}
		}
	}

	removeMarked(H, aO, dO);
	return;
}

bool SameAtom( struct HBondAtom *A,
               struct HBondAtom *B)
{
	// if ( (A->x == B->x) &&
	//      (A->y == B->y) &&
	//      (A->z == B->z) )
	//     return true;

	// Safer to compare strings, rather than floats/doubles.
	if ( !strcmp(A->Name, B->Name) && !strcmp(A->dendrimer,B->dendrimer) )
		return true;

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
            std::vector<struct HBondAtom *> H,
            std::vector<struct HBondAtom *> aO,
            std::vector<struct HBondAtom *> dO,
            unsigned int current)
{
	// H : Hydrogen
	// aO : acceptor Oxygen
	// dO : donor Oxygen
	// DonorO --- Hydrogen ... AcceptorO
	// ... Denotes the Hydrogen bond.
	//
	std::vector<struct HBondAtom *>::iterator iter_Hmain;
	std::vector<struct HBondAtom *>::iterator iter_aOmain;
	std::vector<struct HBondAtom *>::iterator iter_dOmain;

	std::vector<struct HBondAtom *>::iterator iter_H;
	std::vector<struct HBondAtom *>::iterator iter_aO;
	std::vector<struct HBondAtom *>::iterator iter_dO;

	// Check that requested element is not beyond the size
	// of the vector.
	if ( current >= aO.size() )
		return(false);

	iter_aOmain = aO.begin()+current;
	iter_dOmain = dO.begin()+current;
	iter_Hmain = H.begin()+current;

	// If this atom triplet has already been assigned to a chain, skip it
	// Checking for only the H is sufficient.
	if ( ( (*iter_Hmain)->Next != NULL) && ( (*iter_Hmain)->Previous != NULL) )
		return(false);

	// Starting a new chain.
	(*HBonds)->AddAtStart(*iter_aOmain);
	(*HBonds)->AddAtStart(*iter_Hmain);
	(*HBonds)->AddAtStart(*iter_dOmain);

	bool StillLooking = true;
	while ( StillLooking )
	{
		iter_H = H.begin();
		iter_aO = aO.begin();
		iter_dO = dO.begin();
		bool FoundOne = false;
		for( ; iter_H < H.end(); ++iter_H, ++iter_aO, ++iter_dO )
		{
			// Check that this triplet is not in the chain, already.
			// Checking for only the H is sufficient.
			if ( !(*HBonds)->Find(*iter_H)  )
			{
				if ( (*HBonds)->IsSameAsFirst(*iter_aO) )
				{
					// Found a new link at the beginning of the chain.
					(*HBonds)->AddAtStart(*iter_H);
					(*HBonds)->AddAtStart(*iter_dO);
					(*iter_aO)->Next = (*iter_aO)->Previous = NULL;
					FoundOne = true;
				} 
				else if ( (*HBonds)->IsSameAsLast(*iter_dO) )
				{
					// Found a new link at the end of the chain.
					(*HBonds)->AddAtEnd(*iter_H);
					(*HBonds)->AddAtEnd(*iter_aO);
					(*iter_dO)->Next = (*iter_dO)->Previous = NULL;
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
                     bool POVRAY, struct PBC *Cell)
{
	unsigned int MaxChainLength = 0;
	unsigned int MaxLoopSize = 0;
	double EndToEndLength;
	//Go through the vector of HBond strings and tabulate:
	//  Chain Lengths
	//  Closed Loops
	//  Molecule Switches
	//  Molecules
	//
	// Zero all histogram bins
	// Set 20 elements initially.
	std::vector<unsigned int>ChainLength_hist(20,0);
	std::vector<unsigned int>ClosedLoop_hist(20,0);
	
	std::vector< std::vector<unsigned int> >SwitchesInChain_hist(
	                                      20,
	                                      std::vector<unsigned int>(20,0)
	                                      );
	std::vector< std::vector<unsigned int> >MoleculesInChain_hist(
	                                      20,
	                                      std::vector<unsigned int>(20,0)
	                                      );

	std::vector<unsigned int>MaxSwitchesInChain(20,0);
	std::vector<unsigned int>MaxMoleculesInChain(20,0);
	// All histograms zeroed.

	for( unsigned int i=0; i < HBStrings.size(); i++ )
	{
		unsigned int HBCount = HBStrings[i]->Count();
		unsigned int SwitchingCount = HBStrings[i]->SwitchingCount();
		unsigned int MoleculeCount = HBStrings[i]->MoleculeCount();

		// Tabulate the chain lengths.
		if ( alloc_vector(&ChainLength_hist,0, HBCount+1) )
			ChainLength_hist.at(HBCount)++;
		else
			return 1;

		if ( HBCount > MaxChainLength)
			MaxChainLength = HBCount;

		// If this is a closed loop, tabulate it.
		if ( HBStrings[i]->ClosedLoop() )
		{
			if ( alloc_vector(&ClosedLoop_hist, 0, HBCount+1) )
				ClosedLoop_hist.at(HBCount)++;
			else
				return 1;

			if ( HBCount > MaxLoopSize )
				MaxLoopSize = HBCount;
		}

		// Tabulate the number of molecule switches for each chain length.
		if( alloc_vector(&SwitchesInChain_hist, 0, HBCount+1, SwitchingCount+1))
		{
			SwitchesInChain_hist[HBCount][SwitchingCount]++;
		}
		else
			return 1;

		if( !alloc_vector(&MaxSwitchesInChain, 0, HBCount+1) )
			return 1;

		if (SwitchingCount > MaxSwitchesInChain[HBCount])
			MaxSwitchesInChain[HBCount] = SwitchingCount;

		// Tabulate the number of molecules in each chain length.
		if( alloc_vector(&MoleculesInChain_hist, 0, HBCount+1, MoleculeCount+1))
			MoleculesInChain_hist[HBCount][MoleculeCount]++;
		else
			return 1;

		if( !alloc_vector(&MaxMoleculesInChain, 0, HBCount+1) )
			return 1;

		if (MoleculeCount > MaxMoleculesInChain[HBCount])
			MaxMoleculesInChain[HBCount] = MoleculeCount;
	}
	OFmt colE2E(0,6);
	// Printout information about each hbond string.
	for( unsigned int i=0; i < HBStrings.size(); i++ )
	{
		*out << std::endl << std::endl;
		*out << CC << " Current Element : " << i << std::endl;
		*out << CC << " Atoms in Chain : " << HBStrings[i]->Count();
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
		EndToEndLength = HBStrings[i]->PrintAll(out, *Cell, POVRAY);
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
	                     ChainLength_hist,
	                     MaxChainLength, 3,
	                     2, MaxBarLength,
	                     NumBins,CC);

	// Closed Loop histogram table header
	*out << CC <<std::endl << CC 
	          << " Atoms/HBonds |Count| (For Closed Loops)" << std::endl;
	// Minimum Chain length is 3 atoms (O-H...O). Chain length is always an odd 
	// number of atoms, so use a step length of 2.
	PrintHistogramChain( out,
	                     ClosedLoop_hist,
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
	// if ( MaxChainLength < NumBins )
	//     MaxChainL = NumBins;
	// else
	//     MaxChainL = MaxChainLength;

	for(unsigned int chainL=3; chainL <= MaxChainL; chainL += 2)
	{
		// Number of Switches histogram table header
		*out << CC << std::endl << CC 
		          << " Switching |Count| (For Chain length of " << chainL << ")"
		          << std::endl;

		if ( chainL > MaxChainLength )
		{
			std::vector<unsigned int>dummy;
			PrintHistogramMolecules( out,
			                         dummy,
			                         0, (unsigned int)0,
			                         (unsigned int)1, MaxBarLength,
			                         NumBins, CC);
		}
		else
			PrintHistogramMolecules( out,
									 SwitchesInChain_hist[chainL],
									 MaxSwitchesInChain[chainL],(unsigned int)0,
									 (unsigned int)1, MaxBarLength,
									 NumBins, CC);
	}

	// Minimum Chain length is 3 atoms (O-H..O). Chain length is always an odd 
	// number of atoms, so increase counter by 2 for each step.
	// if ( MaxChainLength < NumBins )
	//     MaxChainL = NumBins;
	// else
	//     MaxChainL = MaxChainLength;

	for(unsigned int chainL=3; chainL <= MaxChainL; chainL += 2)
	{
		// Number of Molecules histogram table header
		*out << CC << std::endl << CC 
		          << " Molecules |Count| (For Chain length of " << chainL << ")"
		          << std::endl;

		if ( chainL > MaxChainLength)
		{
			std::vector<unsigned int>dummy;
			PrintHistogramMolecules( out,
									 dummy,
									 0, 1,
									 1, MaxBarLength,
									 NumBins, CC);
		}
		else
			PrintHistogramMolecules( out,
									 MoleculesInChain_hist[chainL],
									 MaxMoleculesInChain[chainL], 1,
									 1, MaxBarLength,
									 NumBins, CC);
	}
	return(0);
}
