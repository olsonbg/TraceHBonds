#include "sizehistPrint.h"

inline int Round( double v)
{
	double r;
	r = (v > 0.0) ? floor(v + 0.5) : ceil(v - 0.5);

	return(int(r));
}


void PrintHistogramChain( std::ostream *out,
                          vui Histogram,
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
                              vui Histogram,
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
	if ( NumBins > 0 )
		uniformStop = Step*NumBins + Start - (unsigned int)1;
	else
		uniformStop = Stop;

	for (unsigned int i = Start; i <= uniformStop; i += Step )
	{
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


