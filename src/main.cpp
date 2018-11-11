/**
 * \file
 * \brief main()
 */
#include <cstdlib>
#include <getopt.h>
#include "TraceHBonds.h"
#include "help.h"
#include "Thread.h"
#include "queue.h"
#include "WorkerThreads.h"
#include "cpu.h"
#include "flags.h"

/**
 * Flag for being verbose
 */
bool THB_VERBOSE = false;

#ifdef PTHREADS
extern Queue<struct worker_data_s> inQueue;
#endif // PTHREADS

/**
 * main function
 */
int main(int argc, char *argv[])
{
	unsigned int NumBins = 0;
	// Connections filename without the .arc extension for discover, with
	// extension for LAMMPS.
	char *fileData    =  NULL;
	char *fileTrj     =  NULL; // LAMMPS trajectory file.
	char *fileMols    =  NULL; // Molecular definition file for LAMMPS.
	char *ofPrefix    =  NULL; // Prefix for output filename.
	char *ofSuffix    =  NULL; // Suffix for output filename.
	double rCutoff    =   0.0;  // Distance Cutoff for Hydrogen bonds.
	double angleCutoff= 180.0; // Angle Cutoff for Hydrogen bonds.
	// Matching keys for hydrogens and acceptors.
	struct HydrogenBondMatching match;
	int   flag[(unsigned int)Flags::COUNT] = {}; // Initialize all elements of flag to zero.
	unsigned int flags = 0;

	// Read command line arguments.
	int c;
	while(1)
	{
		static struct option long_options[] =
		{
			/* These options set a flag. */
			{"verbose"     , no_argument, &flag[0],  (unsigned int)Flags::VERBOSE      },
			{"brief"       , no_argument, &flag[0],  0            },
			{"povray"      , no_argument, &flag[1],  (unsigned int)Flags::POVRAY       },
			{"lifetime"    , no_argument, &flag[2],  (unsigned int)Flags::LIFETIME     },
			{"lengths"     , no_argument, &flag[3],  (unsigned int)Flags::LENGTHS      },
			{"angles"      , no_argument, &flag[4],  (unsigned int)Flags::ANGLES       },
			{"sizehist"    , no_argument, &flag[5],  (unsigned int)Flags::SIZE_HIST    },
			{"neighborhist", no_argument, &flag[6],  (unsigned int)Flags::NEIGHBOR_HIST},
			{"json"        , no_argument, &flag[7],  (unsigned int)Flags::JSON         },
			{"incell"      , no_argument, &flag[8],  (unsigned int)Flags::INCELL       },
			{"jsonall"     , no_argument, &flag[9],  (unsigned int)Flags::JSONALL      },
			{"all"         , no_argument, &flag[10], (unsigned int)Flags::ALL          },
			{"list"        , no_argument, &flag[11],  (unsigned int)Flags::LIST       },
			/* These options don’t set a flag.
			   We distinguish them by their indices. */
			{"input",       required_argument, 0, 'i'},
			{"trajectory",  required_argument, 0, 't'},
			{"molecules",   required_argument, 0, 'm'},
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

		c = getopt_long (argc, argv, "i:t:m:p:s:f:l:b:r:a:H:A:h",
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
				fileData = optarg;
				break;
			case 't':
				fileTrj = optarg;
				break;
			case 'm':
				fileMols = optarg;
				break;
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
				match.Hydrogens.push_back(optarg);
				break;
			case 'A':
				match.Acceptors.push_back(optarg);
				break;
			case 'h':
			case '?':
			default:
				Help(argv[0]);
				return(1);
		}
	}
	// Done reading command line arguments

	if ( match.Hydrogens.size() == 0 ) {
		Help(argv[0]);
		BRIEF_MSG("\nError: Must specify at least one type of Hydrogen.");
		return EXIT_FAILURE;
	}
	if ( match.Acceptors.size() == 0 ) {
		Help(argv[0]);
		BRIEF_MSG("\nError: Must specify at least one type of Acceptor.");
		return EXIT_FAILURE;
	}
	// Set the flags;
	for ( unsigned int i=0; i < (unsigned int)Flags::COUNT; i++ ) {
		flags |= flag[i];
	}

	if ( flags & Flags::VERBOSE ) THB_VERBOSE=true;

	if ( (fileData == NULL) )
	{
		Help(argv[0]);
		BRIEF_MSG("\nError: Must specify an input file.");
		return EXIT_FAILURE;
	}

	if ( (fileTrj != NULL) && (fileMols == NULL) )
	{
		Help(argv[0]);
		BRIEF_MSG("\nError: When specifying a trajectory file, must also specify a molecule file.");
		return EXIT_FAILURE;
	}

	if ( (fileMols != NULL) && (fileTrj == NULL) )
	{
		Help(argv[0]);
		BRIEF_MSG("\nError: When specifying a molecule file, must also specify an trajectory file.");
		return EXIT_FAILURE;
	}

	if ( (flags & Flags::SIZE_HIST) && (ofPrefix == NULL) )
	{
		Help(argv[0]);
		BRIEF_MSG("\nError: Must specify a prefix for the output file.");
		return EXIT_FAILURE;
	}

	VERBOSE_MSG("\t--- Calculations ---");
	if ( flags & Flags::LIFETIME )     VERBOSE_MSG("\tHydrogen bond lifetime correlations");
	if ( flags & Flags::SIZE_HIST)     VERBOSE_MSG("\tChain lengths in each frame");
	if ( flags & Flags::NEIGHBOR_HIST) VERBOSE_MSG("\t- Consolidated chain lengths");
	if ( flags & Flags::LENGTHS)       VERBOSE_MSG("\tHydrogen - Acceptor distances");
	if ( flags & Flags::ANGLES)        VERBOSE_MSG("\tHydrogen bond angles");
	if ( flags & Flags::LIST)          VERBOSE_MSG("\tHydrogen bond list");
	VERBOSE_MSG("\t--------------------\n");

#ifdef PTHREADS
	VERBOSE_MSG("Starting " << NumberOfCPUs() << " threads.");

	std::vector<MyThread *>MyThreads;
	for (unsigned int j=0; j != NumberOfCPUs(); ++j)
	{
		// Setup and start a thread.
		MyThreads.push_back( new MyThread() );
		MyThreads.at(j)->start();
	}
#endif

	doArcFile(fileData, fileTrj, fileMols,
	          ofPrefix, ofSuffix,
	          &match,
	          rCutoff, angleCutoff,
	          NumBins, flags);

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

	for (unsigned int j=0; j != NumberOfCPUs(); ++j) {
		inQueue.push(wd); }

	// Wait for the threads to return.
	for (unsigned int j=0; j != NumberOfCPUs(); ++j) {
		MyThreads.at(j)->join(); }

	while ( !MyThreads.empty() ) {
		delete MyThreads.back();
		MyThreads.pop_back();
	}
#endif

	return(0);
}

