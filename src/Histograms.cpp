#include "Histograms.h"
#include "Print.h"

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
				std::vector<T> Zero( melem, val);
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

void getNeighbors( struct Histograms_s *Histograms,
                   std::vector<ListOfHBonds *> HBStrings,
                   struct PBC *Cell )
{
	unsigned int TrjIdx = Histograms->TrjIdx;
	Point cellp = Cell->p.at(TrjIdx);
	// Starting index of NearestNeighbor index for a specific chain length.
	unsigned int offset=0;

	for( unsigned int i = 0; i < HBStrings.size(); ++i )
	{
		if ( HBStrings.at(i)->TrajectoryIndex() != TrjIdx )
			continue;


		unsigned int N=HBStrings.at(i)->HydrogenBondCount();
		if (N==0) return;

		offset = N*(N-1)/2;

		std::vector<Point *>pCoord = HBStrings.at(i)->nonHydrogenCoordinates();

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

void Correlations( std::ostream *out,
                   std::vector< std::vector<bool> > *v )
{
	unsigned int numHBs = v->size();
	unsigned int numFrames = v->at(0).size();
	unsigned int fcutoff=numFrames/2;

	// Initialize the histograms to zero (0).
	vvui continuous  ( numHBs, vui(fcutoff, 0) );
	vvui intermittent( numHBs, vui(fcutoff, 0) );

	for( unsigned int h=0; h < numHBs; ++h)
	{
		for( unsigned int f1=0; f1 < fcutoff; ++f1)
		{
			unsigned int f2=f1+1;
			if ( v->at(h).at(f1) ) {
				continuous.at(h).at(0)++;
				intermittent.at(h).at(0)++;
				for( ; f2 < f1+fcutoff; ++f2)
				{
					if( ! v->at(h).at(f2) ) break;

					continuous.at(h).at(f2-f1)++;
					intermittent.at(h).at(f2-f1)++;
				}
				// Continue until fcutoff to calculate the intermittent
				// hydrogen bond autocorrelation.
				for(f2 = f2+1; f2 < f1+fcutoff; ++f2)
				{
					if( v->at(h).at(f2) )
						intermittent.at(h).at(f2-f1)++;
				}

			}
		}
	}

	// for( unsigned int h=0; h < c.size(); ++h)
	// {
		// for( unsigned int i=0; i < c.at(h).size(); ++i)
			// *out << i << "\t" << c.at(h).at(i) << "\n";
		// *out << "\n";
	// }

	// Save the continuous hydrogen bond autocorrelation data.
	double max=0.0;
	for( unsigned int i=0; i < fcutoff; ++i)
	{
		double d=0.0;
		for( unsigned int h=0; h < continuous.size(); ++h) {
			d += continuous.at(h).at(i); }

		d = d/continuous.size();

		if ( i == 0 ) max = d;

		*out << i << "\t" << d/max << "\n";
	}

	*out << "\n";

	// Save the intermittent hydrogen bond autocorrelation data.
	max=0.0;
	for( unsigned int i=0; i < fcutoff; ++i)
	{
		double d=0.0;
		for( unsigned int h=0; h < intermittent.size(); ++h) {
			d += intermittent.at(h).at(i); }

		d = d/intermittent.size();

		if ( i == 0 ) max = d;

		*out << i << "\t" << d/max << "\n";
	}
}

/*
 * Go through the vector of HBond strings and:
 *  Bin Chain Lengths                             (1D)
 *  Bin Chain Lengths of only Closed Loops        (1D)
 *  Bin Molecule Switches for each Chain Length   (2D)
 *  Bin Molecules in Chain, for each Chain Length (2D)
 */
struct Histograms_s
makeHistograms( std::vector<ListOfHBonds *> HBStrings,
                unsigned int TrjIdx)
{
	// Zero all histogram bins. Set 20 elements initially.
	struct Histograms_s Histogram = { TrjIdx,
	                                  vui(20,0), vui(20,0), 0, 0,
	                                  vvui(20,vui(20,0)),
	                                  vvui(20,vui(20,0)),
	                                  vui(20,0), vui(20,0),
	                                  vvd(1,vd(1,-1.0)) };

	for( unsigned int i=0; i < HBStrings.size(); i++ )
	{
		if ( HBStrings[i]->TrajectoryIndex() != TrjIdx )
			continue;

		unsigned int HBCount        = HBStrings[i]->AtomCount();
		unsigned int SwitchingCount = HBStrings[i]->SwitchingCount();
		unsigned int MoleculeCount  = HBStrings[i]->MoleculeCount();

		// Bin the chain lengths.
		if( !Bin(&Histogram.ChainLength, &Histogram.MaxChainLength, HBCount) )
			return Histogram;

		// Bin the chain lengths for only closed loops.
		if ( HBStrings[i]->ClosedLoop() )
		{
			if( !Bin(&Histogram.ClosedLoop, &Histogram.MaxClosedLoop, HBCount) )
				return Histogram;
		}

		// Bin the number of molecule switches for each chain length.
		if( !Bin(&Histogram.SwitchesInChain, &Histogram.MaxSwitchesInChain,
		         HBCount, SwitchingCount) )
			return Histogram;

		// Tabulate the number of molecules in each chain length.
		if ( !Bin(&Histogram.MoleculesInChain,&Histogram.MaxMoleculesInChain,
		          HBCount,MoleculeCount) )
			return Histogram;
	}

	return(Histogram);
}

void
prntHistograms( std::ostream *out,
                std::vector<ListOfHBonds *> HBStrings,
                struct Histograms_s *Histogram,
                std::string CC, unsigned int NumBins,
                struct PBC *Cell, unsigned int TrjIdx,
                bool POVRAY)
{
	double EndToEndLength;

	OFmt colE2E(0,6);

	// Header for povray file.
	if (POVRAY)
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
