#include "Histograms.h"
#include "sizehistPrint.h"

#ifdef PTHREADS
extern Queue<struct worker_data_s> inQueue;
extern Queue<struct worker_data_s> outQueue;
#endif

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

/** \todo Check why this makes nelem*melem matrix, instead of melem's of
 * varying nelem length. In other words, why doesn't each v[melem] have only
 * needed number of elements, instead of nelem for all.
 */

/*
 * Make sure v is of size nelem*melem, if not, initialize the needed number of
 * elements to val. Make sure we don't adjust the values already stored in v.
 */
template<class T> bool alloc_vector( std::vector< std::vector<T> > *v,
                                     T val,
                                     unsigned int melem,
                                     unsigned int nelem)
{
	for ( unsigned int n=0; n < melem; n++)
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
			bool ret = alloc_vector(&(v->at(n)), val, nelem);

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
 * Add to the counts in the h[bin_i][bin_j] h, making sure enough space is
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
		return false;

	if (c > hmax->at(hb))
		hmax->at(hb) = c;
	return true;
}

void getNeighbors( struct Histograms_s *Histograms,
                   std::vector<ListOfHBonds *> *HBStrings,
                   struct PBC *Cell )
{
	unsigned int TrjIdx = Histograms->TrjIdx;
	Point cellp = Cell->p.at(TrjIdx);
	// Starting index of NearestNeighbor index for a specific chain length.
	unsigned int offset=0;

	for( unsigned int i = 0; i < HBStrings->size(); ++i )
	{
		if ( HBStrings->at(i)->TrajectoryIndex() != TrjIdx )
			continue;

		unsigned int N=HBStrings->at(i)->HydrogenBondCount();
		if (N==0) return;

		// Hydrogen Bond count to index
		offset = N*(N-1)/2;

		std::vector<Point *>pCoord = HBStrings->at(i)->nonHydrogenCoordinates();

		// Initialize elements to -1.0.
		alloc_vector(&(Histograms->NearestNeighbors), -1.0,
		             offset+N,1);

		// Make sure the coordinates are the minimum Image of their neighbors.
		std::vector< Point >p;
		p.push_back( *(pCoord.at(0)) );

		for( unsigned int j = 1; j < pCoord.size(); ++j) {
			Point image = p.at(j-1) + pCoord.at(j-1)->minimumImage( *(pCoord.at(j)), cellp );
			p.push_back(image);
		}

		for( unsigned int j = 0; j < p.size()-1; ++j) {
			for (unsigned int k = j+1; k < p.size(); ++k)
			{
				double d = p.at(j).distance( p.at(k) );
				if( Histograms->NearestNeighbors.at(offset+k-j-1).at(0) == -1.0)
					Histograms->NearestNeighbors.at(offset+k-j-1).at(0) = d;
				else
					Histograms->NearestNeighbors.at(offset+k-j-1).push_back(d);
			}
		}
	}
}


/*
 * Go through the vector of HBond strings and:
 *  Bin Chain Lengths                             (1D)
 *  Bin Chain Lengths of only Closed Loops        (1D)
 *  Bin Molecule Switches for each Chain Length   (2D)
 *  Bin Molecules in Chain, for each Chain Length (2D)
 */
void
makeHistograms( struct Histograms_s *Histogram,
                std::vector<ListOfHBonds *> *HBStrings,
                unsigned int TrjIdx)
{
	for( unsigned int i=0; i < HBStrings->size(); i++ )
	{
		if ( HBStrings->at(i)->TrajectoryIndex() != TrjIdx )
			continue;

		unsigned int HBCount        = HBStrings->at(i)->AtomCount();
		unsigned int SwitchingCount = HBStrings->at(i)->SwitchingCount();
		unsigned int MoleculeCount  = HBStrings->at(i)->MoleculeCount();

		// Bin the chain lengths.
		if( !Bin(&Histogram->ChainLength, &Histogram->MaxChainLength, HBCount) )
			return;

		// Bin the chain lengths for only closed loops.
		if ( HBStrings->at(i)->ClosedLoop() )
		{
			if( !Bin(&Histogram->ClosedLoop, &Histogram->MaxClosedLoop, HBCount) )
				return;
		}

		// Bin the number of molecule switches for each chain length.
		if( !Bin(&Histogram->SwitchesInChain, &Histogram->MaxSwitchesInChain,
		         HBCount, SwitchingCount) )
			return;

		// Tabulate the number of molecules in each chain length.
		if ( !Bin(&Histogram->MoleculesInChain,&Histogram->MaxMoleculesInChain,
		          HBCount,MoleculeCount) )
			return;
	}
}

void
prntHistograms( std::ostream *out,
                std::vector<ListOfHBonds *> *HBStrings,
                struct Histograms_s *Histogram,
                std::string CC, unsigned int NumBins,
                struct PBC *Cell, unsigned int TrjIdx,
                unsigned int flags)
{
	double EndToEndLength;

	OFmt colE2E(0,6);

	// Header for povray file.
	if ( flags & Flags::POVRAY)
	{
		*out << "#version 3.6;" << "\n";
		*out << "global_settings {  assumed_gamma 1.0 }" << "\n";
		*out << "Camera_LookAt( " << Cell->p.at(TrjIdx).x() << ", "
		                          << Cell->p.at(TrjIdx).y() << ", "
		                          << Cell->p.at(TrjIdx).z() << " )" << "\n";
		*out << "PBC( " << Cell->p.at(TrjIdx).x() << ", "
		                << Cell->p.at(TrjIdx).y() << ", "
		                << Cell->p.at(TrjIdx).z() << " )" << "\n";
	}

	// If JSON, only print data.
	if ( flags & Flags::JSON )
	{
		*out << "[";
		for( unsigned int i=0; i < HBStrings->size(); i++ )
		{
			if ( HBStrings->at(i)->TrajectoryIndex() != TrjIdx )
				continue;
			HBStrings->at(i)->PrintAll(out, *Cell, TrjIdx, flags);
		}
		return;
	}

	// Printout information about each hbond string.
	for( unsigned int i=0; i < HBStrings->size(); i++ )
	{
		if ( HBStrings->at(i)->TrajectoryIndex() != TrjIdx )
			continue;

		*out << "\n\n";
		*out << CC << " Current Element : " << i << "\n";
		*out << CC << " Atoms in Chain : " << HBStrings->at(i)->AtomCount();
		*out << "\n";

		// Note if this is a closed loop.
		if ( HBStrings->at(i)->ClosedLoop() )
			*out << CC << " Closed Loop" << "\n";

		*out << CC << " Molecules : "
		          << HBStrings->at(i)->MoleculeCount() << "\n";

		*out << CC << " Unique forcefields : "
		          << HBStrings->at(i)->ForcefieldCount() << "\n";

		*out << CC
		          << " Times chain switched between Molecules (switching) : "
		          << HBStrings->at(i)->SwitchingCount() << "\n";

		*out << CC << " Periodic boundary conditions applied."
		          << "\n";
		// Show the Chain atoms, molecules and coordinates
		EndToEndLength = HBStrings->at(i)->PrintAll(out, *Cell, TrjIdx, flags);
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
	                     Histogram->ChainLength,
	                     Histogram->MaxChainLength, 3,
	                     2, MaxBarLength,
	                     NumBins,CC);

	// Closed Loop histogram table header
	*out << CC <<"\n" << CC
	          << " Atoms/HBonds |Count| (For Closed Loops)" << "\n";
	// Minimum Chain length is 3 atoms (O-H...O). Chain length is always an odd
	// number of atoms, so use a step length of 2.
	PrintHistogramChain( out,
	                     Histogram->ClosedLoop,
	                     Histogram->MaxClosedLoop, 3,
	                     2, MaxBarLength,
	                     NumBins,CC);

	// Molecules in each Chain histogram

	// Minimum Chain length is 3 atoms (O-H..O). Chain length is always an odd
	// number of atoms, so increase counter by 2 for each step.
	unsigned int MaxChainL;
	if ( NumBins > 0 )
		MaxChainL = NumBins;
	else
		MaxChainL = Histogram->MaxChainLength;

	for(unsigned int chainL=3; chainL <= MaxChainL; chainL += 2)
	{
		// Number of Switches histogram table header
		*out << CC << "\n" << CC
		          << " Switching |Count| (For Chain length of " << chainL << ")"
		          << "\n";

		if ( chainL > Histogram->MaxChainLength )
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
			                         Histogram->SwitchesInChain[chainL],
			                         Histogram->MaxSwitchesInChain[chainL],0,
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

		if ( chainL > Histogram->MaxChainLength)
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
			                         Histogram->MoleculesInChain[chainL],
			                         Histogram->MaxMoleculesInChain[chainL], 1,
			                         1, MaxBarLength,
			                         NumBins, CC);
	}
}
