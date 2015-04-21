#ifndef _Print_h
#define _Print_h
#include <iomanip>
#include "TraceHBondsConfig.h"
#include "TraceHBonds.h"

void PrintHistogramChain( std::ostream *out,
                          std::vector<unsigned int>Histogram,
                          unsigned int Max, unsigned int Min,
                          unsigned int Step, int BarLength,
                          unsigned int NumBins,
                          std::string CC );

void PrintHistogramMolecules( std::ostream *out,
                              std::vector<unsigned int>Histogram,
                              unsigned int Max, unsigned int Min,
                              unsigned int Step, int BarLength,
                              unsigned int NumBins,
                              std::string CC );

void Help(char *callname);

#endif
