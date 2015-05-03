#include <cstdlib>
#include <getopt.h>
#include "TraceHBonds.h"
#include "Print.h"
#include "Thread.h"
#include "queue.h"
#include "WorkerThreads.h"
#include "cpu.h"

bool THB_VERBOSE = false;

#ifdef PTHREADS
extern Queue<struct worker_data_s> inQueue;
#endif // PTHREADS

int main(int argc, char *argv[])
{
	unsigned int NumBins = 0;
	char *fArc        =  NULL; // arc filename (without the .arc extension).
	char *ofPrefix    =  NULL; // Prefix for output filename.
	char *ofSuffix    =  NULL; // Suffix for output filename.
	double rCutoff    =   0.0;  // Distance Cutoff for Hydrogen bonds.
	double angleCutoff= 180.0; // Angle Cutoff for Hydrogen bonds.
	// Matching keys for hydrogens and acceptors.
	struct HydrogenBondMatching match;
	int   flag_verbose= 0;
	int   POVRAY = 0;

	// Read command line arguments.
	int c;
	while(1)
	{
		static struct option long_options[] =
		{
			/* These options set a flag. */
			{"verbose", no_argument,       &flag_verbose, 1},
			{"brief",   no_argument,       &flag_verbose, 0},
			{"povray",  no_argument,       &POVRAY,       1},
			/* These options don’t set a flag.
			   We distinguish them by their indices. */
			{"input",       required_argument, 0, 'i'},
			{"outprefix",   required_argument, 0, 'p'},
			{"outsuffix",   required_argument, 0, 's'},
			{"bins",        required_argument, 0, 'b'},
			{"rcutoff",     required_argument, 0, 'r'},
			{"anglecutoff", required_argument, 0, 'a'},
			{"hydrogen",    required_argument, 0, 'H'},
			{"acceptor",    required_argument, 0, 'A'},
			{"help",        no_argument,       0, 'h'},
			{0, 0, 0, 0}
		};

		/* getopt_long stores the option index here. */
		int option_index = 0;

		c = getopt_long (argc, argv, "i:p:s:f:l:b:r:a:H:A:h",
		                 long_options, &option_index);

		/* Detect the end of the options. */
		if (c == -1)
			break;

		switch (c)
		{
			case 0:
				/* If this option set a flag, do nothing else now. */
				if (long_options[option_index].flag != 0)
					break;
				std::cout << "Option: " << long_options[option_index].name;
				if (optarg)
					std::cout << " with arg " <<  optarg <<"\n";
				break;
			case 'i':
				fArc = optarg;
			case 'p':
				ofPrefix = optarg;
				break;
			case 's':
				ofSuffix = optarg;
				break;
			case 'b':
				NumBins = atoi(optarg);
				break;
			case 'r':
				rCutoff = atof(optarg);
				break;
			case 'a':
				angleCutoff = atof(optarg);
				break;
			case 'H':
				match.Hydrogens.insert(optarg);
				break;
			case 'A':
				match.Acceptors.insert(optarg);
				break;
			case 'h':
			case '?':
			default:
				Help(argv[0]);
				return(1);
		}
	}

	if ( flag_verbose ) THB_VERBOSE=true;

	// Done reading command line arguments
	// 
	// 

	if ( (fArc == NULL) || (ofPrefix == NULL) )
	{
		Help(argv[0]);
		return(1);
	}

#ifdef PTHREADS
	VERBOSE_MSG("Starting " << NumberOfCPUs() << " threads.");

	std::vector<MyThread *>MyThreads;
	for (unsigned int j=0; j < NumberOfCPUs(); ++j)
	{
		// Setup and start a thread.
		MyThreads.push_back( new MyThread() );
		MyThreads.at(j)->start();
	}
#endif

	doArcFile(fArc, ofPrefix, ofSuffix,
	          &match,
	          rCutoff, angleCutoff,
	          NumBins, POVRAY);

#ifdef PTHREADS
	// Tell the threads to exit.
	// Use a dummy worker_data_s. All that matters is the jobtype. The other
	// parameters are set just to avoid 'uninitialize parameter' warning by
	// the compiler.
	struct worker_data_s wd;
	wd.jobtype = THREAD_JOB_EXIT;
	wd.jobnum  = 0;
	wd.hydrogens = NULL;
	wd.acceptors = NULL;
	wd.hb = NULL;
	wd.TrjIdx = 0;
	wd.num_threads = NumberOfCPUs();
	wd.rCutoff = 0.0;
	wd.angleCutoff = 0.0;

	for (unsigned int j=0; j < NumberOfCPUs(); ++j) {
		inQueue.push(wd); }

	// Wait for the threads to return.
	for (unsigned int j=0; j < NumberOfCPUs(); ++j) {
		MyThreads.at(j)->join(); }

	std::vector<MyThread *>::iterator it=MyThreads.begin();
	for (; it < MyThreads.end(); ++it) {
		delete *it; }
#endif

	return(0);
}

