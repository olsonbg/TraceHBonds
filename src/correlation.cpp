#include "correlation.h"

#ifdef PTHREADS
extern Queue<struct worker_data_s> inQueue;
extern Queue<struct worker_data_s> outQueue;
#endif

void CorrelationsThread(vd *C, vd *I,
                        vvui *continuous, vvui *intermittent,
                        unsigned int NumThreads=1, unsigned int ThreadID=0 )
{
	unsigned int fcutoff = continuous->at(0).size();

	for( unsigned int i=ThreadID; i < fcutoff; i += NumThreads)
	{
		for( unsigned int h=0; h < continuous->size(); ++h) {
			C->at(i) += continuous->at(h).at(i);
			I->at(i) += intermittent->at(h).at(i);
		}

		// Normalize
		C->at(i) = C->at(i)/continuous->size();
		I->at(i) = I->at(i)/intermittent->size();
	}

}

/**
 * Continuous and Intermittent lifetimes for each hydrogen bond
 */
void CorrelationsTableThread( std::vector< std::vector<bool> > *v,
                              vvui *continuous, vvui *intermittent,
                              unsigned int numHBs,
                              unsigned int fcutoff,
                              unsigned int NumThreads=1,
                              unsigned int ThreadID=0)
{
	for( unsigned int h=ThreadID; h < numHBs; h += NumThreads)
	{
		for( unsigned int f1=0; f1 < fcutoff; ++f1)
		{
			unsigned int f2=f1+1;
			if ( v->at(h).at(f1) ) {
				continuous->at(h).at(0)++;
				intermittent->at(h).at(0)++;
				for( ; f2 < f1+fcutoff; ++f2)
				{
					if( ! v->at(h).at(f2) ) break;

					continuous->at(h).at(f2-f1)++;
					intermittent->at(h).at(f2-f1)++;
				}
				// Continue until fcutoff to calculate the intermittent
				// hydrogen bond autocorrelation.
				for(f2 = f2+1; f2 < f1+fcutoff; ++f2)
				{
					if( v->at(h).at(f2) )
						intermittent->at(h).at(f2-f1)++;
				}
			}
		}
	}
}

void Correlations( std::ostream *out,
                   std::vector< std::vector<bool> > *v )
{
	unsigned int numHBs = v->size();
	unsigned int numFrames = v->at(0).size();
	// Using a sliding window, therefore cutoff should be half
	// the total.
	unsigned int fcutoff=numFrames/2;

	VERBOSE_MSG("\tCalculating autocorrelation for all hydrogen bonds.");
	VERBOSE_MSG("\t\tUsing moving window of " << fcutoff << " frames.");
	// Initialize the histograms to zero (0).
	vvui continuous  ( numHBs, vui(fcutoff, 0) );
	vvui intermittent( numHBs, vui(fcutoff, 0) );

#ifdef PTHREADS
	for( unsigned int jobnum=0; jobnum < NumberOfCPUs(); ++jobnum)
	{
		struct worker_data_s wd;
		wd.jobtype     = THREAD_JOB_CORR_TABLE;
		wd.jobnum      = jobnum;
		wd.num_threads = NumberOfCPUs();
		wd.b           = v;
		wd.numHBs      = numHBs;
		wd.fcutoff     = fcutoff;
		wd.vvuiC       = new vvui(numHBs, vui(fcutoff, 0) );
		wd.vvuiI       = new vvui(numHBs, vui(fcutoff, 0) );

		inQueue.push(wd);
	}

	// Get the results back from the worker threads.
	for( unsigned int jobnum=0; jobnum < NumberOfCPUs(); ++jobnum)
	{
		struct worker_data_s wd = outQueue.pop();

		for( unsigned int h=0; h < numHBs; ++h)
		{
			for( unsigned int f = 0; f < fcutoff; ++f)
			{
				continuous.at(h).at(f)   += wd.vvuiC->at(h).at(f);
				intermittent.at(h).at(f) += wd.vvuiI->at(h).at(f);
			}
		}
		delete wd.vvuiC;
		delete wd.vvuiI;
	}
#else
	CorrelationsTableThread( v, &continuous, &intermittent,
	                         numHBs, fcutoff );
#endif // PTHREADS

	VERBOSE_MSG("\tAveraging autocorrelations over all hydrogen bonds.");
	// Combine all hydrogen bonds.
	vd C(fcutoff,0.0);
	vd I(fcutoff,0.0);

#ifdef PTHREADS
	for( unsigned int jobnum=0; jobnum < NumberOfCPUs(); ++jobnum)
	{
		struct worker_data_s wd;
		wd.jobtype = THREAD_JOB_CORR;
		wd.jobnum = jobnum;
		wd.num_threads = NumberOfCPUs();
		wd.vvuiC = &continuous;
		wd.vvuiI = &intermittent;

		wd.vdC = new vd(fcutoff,0.0); // Cthread(fcutoff,0.0);
		wd.vdI = new vd(fcutoff,0.0); // Ithread(fcutoff,0.0);
		inQueue.push(wd);
	}

	// Get the results back from the worker threads.
	for( unsigned int jobnum=0; jobnum < NumberOfCPUs(); ++jobnum)
	{
		struct worker_data_s wd = outQueue.pop();
		for ( unsigned int i=0; i < wd.vdC->size(); ++i )
		{
			C.at(i) += wd.vdC->at(i);
			I.at(i) += wd.vdI->at(i);
		}
		delete wd.vdC;
		delete wd.vdI;
	}
#else
	CorrelationsThread( &C, &I, &continuous, &intermittent );
#endif // PTHREADS

	BRIEF_MSG("\tSaving autocorrelation data.");

	// Save the continuous and intermittent hydrogen bond autocorrelation data.
	*out << "# Continuous Intermittent\n";

	for( unsigned int i=0; i < fcutoff; ++i) {
		*out << i << "\t" << C.at(i)/C.at(0) << "\t" << I.at(i)/I.at(0) << "\n";
	}
}

