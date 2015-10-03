#include "ListOfHBonds.h"
#include "MessageDefines.h"
#include "Trace.h"
#include "Histograms.h"
#include "sizehist.h"

void sizehist(unsigned int NumFramesInTrajectory,
              char *ofPrefix, char *ofSuffix,
              HBVec *hb,
              struct PBC *Cell,
              int NumBins,
              bool povray,
              std::vector<ListOfHBonds *>* HBStrings)
{
	time_t timer = time(NULL);
	unsigned int TrjIdx;


	// Make histograms.

	VERBOSE_MSG("Generating size histograms.");
	std::vector<struct Histograms_s> Histograms;
	Histograms.reserve(NumFramesInTrajectory);

	for( TrjIdx = 0 ; TrjIdx != NumFramesInTrajectory; ++TrjIdx ) {
		// Zero all histogram bins. Set 20 elements initially.
		struct Histograms_s Histogram = { TrjIdx,
		                                  vui(20,0), vui(20,0), 0, 0,
		                                  vvui(20,vui(20,0)),
		                                  vvui(20,vui(20,0)),
		                                  vui(20,0), vui(20,0),
		                                  vvd(1,vd(1,-1.0)) };

		makeHistograms(&Histogram, HBStrings, TrjIdx );
		// Histograms are not guaranteed to be in the same order as the frames;
		// check Histograms.at(i).TrjIdx to know which frame Histograms.at(i)
		// corresponds to.
		Histograms.push_back( Histogram );
	}

	// Save the histograms.
	const char *CC1 = "#";
	const char *CC2 = "//";
	std::string CC;

	// Povray uses a different comment string.
	if ( povray )
		CC = CC2;
	else
		CC = CC1;

	for( TrjIdx = 0 ; TrjIdx != NumFramesInTrajectory; ++TrjIdx )
	{
		std::stringstream ofilename;
		ofilename << ofPrefix << TrjIdx+1 << ofSuffix;


		std::ofstream out;
		out.open(ofilename.str().c_str(),std::ios::out);
		if ( out.is_open() )
		{
			if ( (difftime(time(NULL), timer) > 1.0) ||
				 (TrjIdx == NumFramesInTrajectory-1) )
			{
				BRIEF_RMSG("Saving size histograms, frame " << TrjIdx+1 << "/" << NumFramesInTrajectory << " (" << ofilename.str() << ")");
				timer = time(NULL);
			}
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

			out <<CC<< " Donor Oxygen atoms    : " << hb->size() << "\n";
			out <<CC<< " Hydrogen atoms        : " << hb->size() << "\n";
			out <<CC<< " Acceptor Oxygen atoms : " << hb->size() << "\n";

			// print the histograms and chains, need to find the histogram that corresponds to
			// TrjIdx first.
			for( unsigned int i = 0 ; i != NumFramesInTrajectory; ++i )
			{
				if ( Histograms.at(i).TrjIdx == TrjIdx )
				{
					prntHistograms( &out, HBStrings, &Histograms.at(i),
					                CC, NumBins, Cell, TrjIdx, povray );
					break;
				}
			}

			out.close();
		} else {
			BRIEF_MSG("ERROR: Can not save " <<ofilename.str() << "!" );
		}
	}
	BRIEF_MSG("");
}
