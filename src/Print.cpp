#include "Print.h"

inline int Round( double v)
{
	double r;
	r = (v > 0.0) ? floor(v + 0.5) : ceil(v - 0.5);

	return(int(r));
}


void PrintHistogramChain( std::ostream *out,
                          std::vector<unsigned int>Histogram, 
                          unsigned int Stop, unsigned int Start,
                          unsigned int Step, int BarLength,
                          unsigned int NumBins, std::string CC )
{
	double Scale;
	unsigned int Max = 0;

	if (Step == 0) // That does not make sense, abort!
	{
		std::cerr << "Step should be greater than zero (0)!" << "\n";
		return;
	}

	// Find Max of Histogram
	for (unsigned int i = Start; i <= Stop; i += Step)
		if (Histogram[i] > Max) { Max = Histogram[i]; }

	// Define 'Scale' so that the maximum of the histogram will have a length 
	// of 'BarLength'
	Scale = double(Max)/double(BarLength);

	// Using uniformStop to make it easy to aggregate
	// many TraceHBond output files
	unsigned int uniformStop;
	if ( NumBins > 0 )
		uniformStop = Step*NumBins + Start- (unsigned int)1;
	else
		uniformStop = Stop;

	for (unsigned int i = Start;
	     i <= uniformStop;
	     i += Step)
	{
//		fprintf(*out,"%s %4d /%4d   | %4d|",CC.c_str(),i,(i-1)/2,Histogram[i]);
		*out << CC.c_str() << " ";
		*out << std::setw(4) << i << " /";
		*out << std::setw(4) << (i-1)/2 << "   | ";
		if ( i >= Histogram.size() )
		{
			// If uniformStop brings us past the size of Histogram,
			// show zeros.
			*out << std::setw(4) << 0 << "|";
		}
		else
		{
			*out << std::setw(4) << Histogram[i] << "|";
			for (int j = 0; j < Round(Histogram[i]/Scale); j++) *out << "*";
		}
		*out << "\n";
	}
}

void PrintHistogramMolecules( std::ostream *out,
                              std::vector<unsigned int>Histogram,
                              unsigned int Stop, unsigned int Start,
                              unsigned int Step, int BarLength,
                              unsigned int NumBins,
                              std::string CC )
{
	double Scale;
	unsigned int Max = 0;

	if (Step == 0) // That does not make sense, abort!
	{
		std::cerr << "Step should be greater than zero (0)!" << "\n";
		return;
	}

	// Find Max of Histogram
	if (Histogram.size() > 0 )
	{
		for (unsigned int i = Start; i <= Stop; i += Step)
			if (Histogram[i] > Max) { Max = Histogram[i]; }
	}
	else
	{
		Max = 0;
	}

	// Define 'Scale' so the the maximum of the histogram will have a length 
	// of 'BarLength'
	Scale = double(Max)/double(BarLength);

	// Using uniformStop to make it easy to aggregate
	// many TraceHBond output files
	unsigned int uniformStop;
	// if ( Stop < Step*NumBins + Start - (unsigned int)1)
	if ( NumBins > 0 )
		uniformStop = Step*NumBins + Start - (unsigned int)1;
	else
		uniformStop = Stop;

	for (unsigned int i = Start; i <= uniformStop; i += Step )
	{
//		printf("%s %4d      | %4d|",CC.c_str(),i,Histogram[i]);
		*out << CC.c_str() << " ";
		*out << std::setw(4) << i << "      | ";
		if ( i >= Histogram.size() )
		{
			// If uniformStop brings us past the size of Histogram,
			// show zeros.
			*out << std::setw(4) << 0 << "|";
		}
		else
		{
			*out << std::setw(4) << Histogram[i] << "|";
			for (int j = 0; j < Round(Histogram[i]/Scale); j++) *out << "*";
		}
		*out << "\n";
	}
}

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
				// std::cerr << "Newline" <<"\n";
				string[k] = '\n';
				// std::cerr << "[" << string << "]" <<"\n";
				// set string index to character after this one
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
			for( unsigned int i = 0; i < indent; i++)
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
	// fprintf(stdout,"\n%s Version %d.%d\n\n",
	std::cerr << "\n" << name << " Version ";
	std::cerr << TraceHBonds_VERSION_MAJOR << ".";
	std::cerr << TraceHBonds_VERSION_MINOR;
	std::cerr << "\n";
	std::cerr << "\n";

	std::cerr << "USAGE: " << "\n";
	const char *usage = "(-a <arc file> | -p <prefix> -s <suffix> -f <number> -l <number>) -P <prefix> -S <suffix> [-b <number>] [--povray]";
	std::string USAGE = name;
	USAGE += " ";
	USAGE += usage;
	HelpOption(USAGE.c_str(),
	           NULL,
	           NULL);
	std::cerr << "OPTIONS:" << "\n";
	HelpOption("--arc <arc file>",
	           "-a <arc file>",
	           "<Arc file> is the archive file generated from Discover.");
	HelpOption("--inprefix <prefix>",
	           "-p <prefix>",
	           "<Prefix> of the input filename. The text before the integer in the filename. For a filename of 'HBonds1.dat' the <prefix> would be 'HBonds'");
	HelpOption("--insuffix <suffix>",
	           "-s <suffix>",
	           "<Suffix> of the input filename. The text after the integer in the filename. For a filename of 'HBonds1.dat' the <suffix> would be '.dat'");
	HelpOption("--outprefix <prefix>",
	           "-P <prefix>",
	           "<Prefix> of the output filename. The text before the integer in the filename. For a filename of 'HBonds1.txt' the <prefix> would be HBonds'");
	HelpOption("--outsuffix suffix",
	           "-S suffix",
	           "<Suffix> of the output filename. The text after the integer in the filename. For a filename of 'HBonds1.txt' the <suffix> would be '.dat'");
	HelpOption("--first <number>",
	           "-f <number>",
	           "The first <number> to start the processing at. <number> is the integer in the filename. For filenames of 'HBonds1.dat,' 'HBonds2.dat,' ..., 'HBonds1000.dat,' the first <nmber> would be 1.");
	HelpOption("--last <number>",
	           "-l <number>",
	           "The last <number> to start the processing at. <number> is the integer in the filename. For a filenames of 'HBonds1.dat,' 'HBonds2.dat,' ..., 'HBodns1000.dat,' the last <number> would be 1000.");
	HelpOption("--bins <number>",
	           "-b <number>",
	           "Minimum <number> of bins to show in histograms.");
	HelpOption("--povray",
	           NULL,
	           "Output in povray format.");
	HelpOption("--verbose",
	           NULL,
	           "Show some extra information while running.");
	HelpOption("-h",
	           NULL,
	           "This help screen");
	std::cerr << "\n";
	std::cerr << "Compiled on " << __DATE__;
	std::cerr << " at " << __TIME__ << "." <<"\n";
	std::cerr << "Author: Brian G. Olson (olsonbg@gmail.com)" <<"\n";
}

