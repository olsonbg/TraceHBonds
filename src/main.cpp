#include <cstdlib>
#include <getopt.h>
#include "TraceHBonds.h"
#include "Print.h"

bool THB_VERBOSE = false;

int main(int argc, char *argv[])
{
	unsigned int NumBins = 0;
	char *fArc        = NULL; // arc filename (without the .arc extension).
	char *fPrefix     = NULL; // Prefix for filename.
	char *fSuffix     = NULL; // Suffix for filename.
	int   fidx_f      = 0;    // First number for file index.
	int   fidx_l      = 0;    // Last number for file index.
	bool  flag_fidx_f = false;// Flag indicating whether fidx_f is defined.
	bool  flag_fidx_l = false;// Flag indicating whether fidx_l is defined.
	char *ofPrefix    = NULL; // Prefix for output filename.
	char *ofSuffix    = NULL; // Suffix for output filename.
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
			{"arc",       required_argument, 0, 'a'},
			{"inprefix",  required_argument, 0, 'p'},
			{"insuffix",  required_argument, 0, 's'},
			{"outprefix", required_argument, 0, 'P'},
			{"outsuffix", required_argument, 0, 'S'},
			{"first",     required_argument, 0, 'f'},
			{"last",      required_argument, 0, 'l'},
			{"bins",      required_argument, 0, 'b'},
			{"help",      no_argument,       0, 'h'},
			{0, 0, 0, 0}
		};

		/* getopt_long stores the option index here. */
		int option_index = 0;

		c = getopt_long (argc, argv, "a:p:s:P:S:f:l:b:h",
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
					std::cout << " with arg " <<  optarg <<std::endl;
				break;
			case 'a':
				fArc = optarg;
			case 'p':
				fPrefix = optarg;
				break;
			case 's':
				fSuffix = optarg;
				break;
			case 'P':
				ofPrefix = optarg;
				break;
			case 'S':
				ofSuffix = optarg;
				break;
			case 'f':
				fidx_f = atoi(optarg);
				flag_fidx_f = true;
				break;
			case 'l':
				fidx_l = atoi(optarg);
				flag_fidx_l = true;
				break;
			case 'b':
				NumBins = atoi(optarg);
				break;
			case 'h':
			case '?':
			default:
				Help(argv[0]);
				return(1);
		}
	}

	if ( (fArc == NULL) &&
	    ( (flag_fidx_f == false) ||
	     (flag_fidx_l == false) ||
	     (fidx_l < fidx_f)      ||
	     (fidx_l < 0)           ||
	     (fidx_f < 0)          ) )
	{
		Help(argv[0]);
		return(1);
	}

	if ( flag_verbose ) THB_VERBOSE=true;
	// Done reading command line arguments
	doAllFiles(argv[0],
	           fPrefix , fSuffix, fidx_f, fidx_l,
	           ofPrefix,ofSuffix,
	           NumBins, POVRAY);


	return(0);
}

