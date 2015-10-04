#include "help.h"

char * word_wrap (char* string, int line_width) {
	int i = 0;
	int k, counter;

	while( (i+line_width-1) < (int)strlen( string ) )
	{
		counter = i+line_width-1;

		// check for nearest whitespace back in string
		for ( k = counter; k > i; k--)
		{
			if ( isspace( string[k] ) )
			{
				string[k] = '\n';
				i = k + 1;
				break;
			}
		}
		if ( k == i ) // No white space found, then search forward.
		{
			for ( k = counter; k < (int)strlen(string); k++)
			{
				if ( isspace( string[k] ) )
				{
					string[k] = '\n';
					i = k+1;
					break;
				}
			}
			if ( k >= (int)strlen(string) )
				return(string);
		}
	}

	return(string);
}

void HelpIndented(const char *text, unsigned int indent, unsigned int wrap)
{
	char wrapped_text[strlen(text)+1];
	strcpy(wrapped_text,text);

	word_wrap(wrapped_text,wrap-indent);

	for( unsigned int i = 0; i < indent; i++)
		std::cerr << " ";

	for(unsigned int i=0; i < strlen(wrapped_text); i++)
	{
		std::cerr << wrapped_text[i];
		if ( wrapped_text[i] == '\n')
		{
			for( unsigned int ii = 0; ii < indent; ii++)
				std::cerr << " ";
		}
	}
	std::cerr << "\n";
}

void HelpOption(const char *v1, const char *v2, const char *text)
{
	int INDENT=8;
	int WRAP=76;

	if (v1 != NULL) HelpIndented(v1, INDENT, WRAP);
	if (v2 != NULL) HelpIndented(v2, INDENT, WRAP);
	if (text != NULL) HelpIndented(text, 2*INDENT, WRAP);

	std::cerr << "\n";
}

void Help(char *name)
{
	std::cerr << "\n" << name << " Version ";
	std::cerr << TraceHBonds_VERSION_MAJOR << ".";
	std::cerr << TraceHBonds_VERSION_MINOR;
	std::cerr << "\n";
	std::cerr << "\n";

	std::cerr << "USAGE: " << "\n";
	const char *usage = "-i <arc file> -p <prefix> -s <suffix> -r <distance cutoff> -a <angle cutoff> -H <hydrogen forcefield> -A <acceptor forcefield> [-b <number>] [--verbose] [--brief] [--povray] [--lifetime] [--lengths] [--angles] [--sizehist] [--neighborhist] [--all]";
	std::string USAGE = name;
	USAGE += " ";
	USAGE += usage;
	HelpOption(USAGE.c_str(),
	           NULL,
	           NULL);
	std::cerr << "OPTIONS:" << "\n";
	HelpOption("--input <arc file>",
	           "-i <arc file>",
	           "The archive file generated from Discover. ");
	HelpOption("--outprefix <prefix>",
	           "-p <prefix>",
	           "All output will have this string as a prefix to the filenames. For example, to save data as `HBonds1.dat`, use `-p HBonds` as the prefix");
	HelpOption("--outsuffix <suffix>",
	           "-s <suffix>",
	           "All output will have this string as a suffix to the filenames. For example, to save data as 'HBonds1.dat', use `-s .dat` as the suffix");
	HelpOption("--rcutoff <Rc>",
	           "-r <Rc>",
	           "Set the cutoff length, in angstroms, for the determination of a hydrogen bond. ");
	HelpOption("--anglecutoff <Ac>",
	           "-a <Ac>",
	           "Set the cutoff angle, in degrees, for the determination of a hydrogen bond.");
	HelpOption("--hydrogen <force field>",
	           "-H <force field>",
	           "Set the force field of donor hydrogens for hydrogen bonding (e.g. -H h1o). More than one force field may be used by specifying additional -H force field parameters.  **NOTE** the short option is a capital 'H.'");
	HelpOption("--acceptor <force field>",
	           "-A <force field>",
	           "Set the force field of acceptor atoms for hydrogen bonding. More than one force field may be used by specifying additional -A force field parameters (e.g. -A o2h -A o1=). **NOTE** the short option is a capital 'A.'");
	HelpOption("--bins <number>",
	           "-b <number>",
	           "Minimum number of bins to show in histograms.");
	HelpOption("--every <number>",
	           "-e <number",
	           "Load every <number> frames from trajectory (e.g. -e 10)");
	HelpOption("--povray",
	           NULL,
	           "Output in povray format, relevant for --sizehist only.");
	HelpOption("--verbose",
	           NULL,
	           "Show verbose messages while running.");
	HelpOption("--brief",
	           NULL,
	           "Show brief messages while running.");
	HelpOption("--lifetime",
	           NULL,
	           "Calculate hydrogen bond lifetime correlations.");
	HelpOption("--lengths",
	           NULL,
	           "Save length of all hydrogen bonds.");
	HelpOption("--angles",
	           NULL,
	           "Save angle of all hydrogen bonds.");
	HelpOption("--sizehist",
	           NULL,
	           "Save hydrogen bond strings and histograms.");
	HelpOption("--neighborhist",
	           NULL,
	           "Save neighbor length lists.");
	HelpOption("--all",
	           NULL,
	           "Do all calculations and save all data.");
	HelpOption("--help",
	           "-h",
	           "This help screen");
	std::cerr << "\n";
	std::cerr << "Compiled on " << __DATE__;
	std::cerr << " at " << __TIME__ << "." <<"\n";
	std::cerr << "Author: Brian G. Olson (olsonbg@gmail.com)" <<"\n";
}

