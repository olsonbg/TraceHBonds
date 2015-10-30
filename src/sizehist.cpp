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
              unsigned int flags,
              std::vector<ListOfHBonds *>* HBStrings)
{
	time_t timer = time(NULL);
	unsigned int TrjIdx;


	// Make histograms.

	// Zero all histogram bins. Set 20 elements initially.
	struct Histograms_s Histogram =  { 0,
	                                  vui(20,0), vui(20,0), 0, 0,
	                                  vvui(20,vui(20,0)),
	                                  vvui(20,vui(20,0)),
	                                  vui(20,0), vui(20,0),
	                                  vvd(1,vd(1,-1.0)) };
	// Fill Histograms
	std::vector<struct Histograms_s> Histograms(NumFramesInTrajectory,
	                                            Histogram);

	VERBOSE_MSG("Generating size histograms.");
	for( TrjIdx = 0 ; TrjIdx != NumFramesInTrajectory; ++TrjIdx ) {
		Histograms.at(TrjIdx).TrjIdx = TrjIdx;

#ifdef PTHREADS
		struct worker_data_s wd;
		wd.jobtype     = THREAD_JOB_SIZEHIST;
		wd.jobnum      = TrjIdx;
		wd.HBStrings   = HBStrings;
		wd.TrjIdx      = TrjIdx;
		wd.Histogram   = &Histograms.at(TrjIdx);

		inQueue.push(wd);
#else
		makeHistograms(&Histograms.at(TrjIdx), HBStrings, TrjIdx );
#endif
	}

#ifdef PTHREADS
	// Wait for all the worker threads to finish.
	for( TrjIdx = 0 ; TrjIdx != NumFramesInTrajectory; ++TrjIdx ) {
		outQueue.pop();

		if (  (difftime(time(NULL),timer) > 1.0) ||
		      ((TrjIdx+1)==NumFramesInTrajectory)  )
		{
			VERBOSE_RMSG("\tFrame "
			             << TrjIdx+1
			             <<"/"
			             << NumFramesInTrajectory
			             << ".");
			timer = time(NULL);
		}
	}
#endif
	VERBOSE_MSG("");

	// Save the histograms.
	const char *CC1 = "#";
	const char *CC2 = "//";
	std::string CC;

	// Povray uses a different comment string.
	if ( flags & Flags::POVRAY )
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
			if ( ! (flags & Flags::JSON) ) 
			{
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
			}

			// print the histograms and chains
			prntHistograms( &out, HBStrings, &Histograms.at(TrjIdx),
			                CC, NumBins, Cell, TrjIdx, flags );

			if ( flags & Flags::JSON )
			{
				out << "{\n    \"PBC\": [\n"
				    << "        { \"xyz\": ["
				    << Cell->p.at(TrjIdx).x()      << ", "
				    << Cell->p.at(TrjIdx).y()      << ", "
				    << Cell->p.at(TrjIdx).z()      << " "
				    << "] },\n"
				    << "        { \"angles\": ["
				    << Cell->angles.at(TrjIdx).x() << ", "
				    << Cell->angles.at(TrjIdx).y() << ", "
				    << Cell->angles.at(TrjIdx).z()
				    << "] }\n"
				    << "    ]\n}]\n";
			}
			out.close();
		} else {
			BRIEF_MSG("ERROR: Can not save " <<ofilename.str() << "!" );
		}
	}
	BRIEF_MSG("");
}
