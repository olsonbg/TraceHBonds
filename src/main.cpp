#include <cstdlib>
#include <getopt.h>
#include "TraceHBonds.h"
#include "Print.h"
#include "WorkerThreads.h"

bool THB_VERBOSE = false;

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

	unsigned int numCPUs = sysconf( _SC_NPROCESSORS_ONLN );

	VERBOSE_MSG("Using " << numCPUs+1 << " threads.");

	if ( !StartWorkerThreads(numCPUs+1) )
	{
		std::cout << "Error; Could not start worker threads" << std::endl;
		return(1);
	}
	// Pause the threads for now.
	PauseWorkerThreads();
#endif
	doArcFile(fArc, ofPrefix, ofSuffix,
	          &match,
	          rCutoff, angleCutoff,
	          NumBins, POVRAY);
#ifdef PTHREADS
	StopWorkerThreads();
#endif

	return(0);
}

