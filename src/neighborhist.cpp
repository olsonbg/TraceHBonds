#include "ListOfHBonds.h"
#include "MessageDefines.h"
#include "Trace.h"
#include "neighborhistPrint.h"
#include "neighborhist.h"

// Use the same stream but write to different file
template <typename Stream>
void reopen(Stream& pStream,
            const char * pFile, const char *mFile, const char *sFile,
            std::ios_base::openmode pMode = std::ios::out)
{
	std::stringstream ofilename;
	ofilename << pFile << mFile << sFile;

	pStream.close();
	pStream.clear();
	pStream.open(ofilename.str().c_str(), pMode);
	BRIEF_MSG("\tSaving " << ofilename.str() << ".");
}

void neighborhist(unsigned int NumFramesInTrajectory,
                  char *ofPrefix, char *ofSuffix,
                  struct PBC *Cell,
                  std::vector<ListOfHBonds *>*HBStrings)
{
//	time_t timer = time(NULL);
	unsigned int TrjIdx;


	// Make histograms.

	VERBOSE_MSG("Generating size histograms.");
	std::vector<struct Histograms_s> Histograms(NumFramesInTrajectory);

	VERBOSE_MSG("Generating neighbor histograms.");
	for( TrjIdx = 0 ; TrjIdx < NumFramesInTrajectory; ++TrjIdx ) {
		Histograms.at(TrjIdx).TrjIdx = TrjIdx;

#ifdef PTHREADS
		struct worker_data_s wd;
		wd.jobtype     = THREAD_JOB_NEIGHBORHIST;
		wd.jobnum      = TrjIdx;
		wd.HBStrings   = HBStrings;
		wd.Cell        = Cell;
		wd.Histogram   = &( Histograms.at(TrjIdx) );

		inQueue.push(wd);
#else
		getNeighbors( &(Histograms.at(TrjIdx)), HBStrings, Cell );
#endif
	}

#ifdef PTHREADS
	// Wait for all the worker threads to finish.
	for( TrjIdx = 0 ; TrjIdx != NumFramesInTrajectory; ++TrjIdx ) {
		outQueue.pop();
	}
#endif

	std::stringstream ofilename;
	ofilename << ofPrefix << "-NN-AllFrames" << ofSuffix;

	std::ofstream out(ofilename.str().c_str(),std::ios::out);
	if ( out.is_open() )
	{
		BRIEF_MSG("\tSaving " << ofilename.str() << ".");
		Print_AllFrames(&out, &Histograms);
	} else {
		BRIEF_MSG("ERROR: Can not save file!"); }

	reopen(out, ofPrefix, "-NN-Combined", ofSuffix, std::ios::out);
	if ( out.is_open() ) {
		Print_CombineFrames(&out, &Histograms);
	} else {
		BRIEF_MSG("ERROR: Can not save file!"); }

	reopen(out, ofPrefix, "-NN-only", ofSuffix, std::ios::out);
	if ( out.is_open() ) {
		Print_CombineChains(&out, &Histograms);
	} else {
		BRIEF_MSG("ERROR: Can not save file!"); }
}
