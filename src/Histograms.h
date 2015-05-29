#ifndef _Histograms_h
#define _Histograms_h

#include "ListOfHBonds.h"
#include "queue.h"
#include "WorkerThreads.h"
#include "cpu.h"
// Macros

typedef std::vector<unsigned int> vui;
typedef std::vector< vui > vvui;
typedef std::vector< double > vd;
typedef std::vector< vd > vvd;


struct Histograms_s
{
	unsigned int TrjIdx;
	/*
	 * Chain Lengths                             (1D)
	 * Chain Lengths of only Closed Loops        (1D)
	 * Molecule Switches for each Chain Length   (2D)
	 * Molecules in Chain, for each Chain Length (2D)
	 */
	vui ChainLength;
	vui ClosedLoop;
	unsigned int MaxChainLength;
	unsigned int MaxClosedLoop;
	
	vvui SwitchesInChain;
	vvui MoleculesInChain;
	vui MaxSwitchesInChain;
	vui MaxMoleculesInChain;

	// vui ChainLengths;
	/*
	 * NearestNeighbors (NN):
	 * 1st element, NN[0]:vector of nearest neighbor distances.
	 * 2nd element, NN[1]:vector of next nearest neighbor distances.
	 * 3rd element, NN[2]:vector of nextnext nearest neighbor distances.
	 * ...
	 * nth element, NN[n]:vector of nth nearest neighbor distances.
	 */
	std::vector< std::vector<double> >NearestNeighbors;
};

struct Histograms_s
makeHistograms( std::vector<ListOfHBonds *> HBStrings,
                unsigned int TrjIdx);

void Correlations( std::ostream *out,
                   std::vector< std::vector<bool> > *v );

void CorrelationsThread(vd *C, vd *I, 
                        vvui *continuous, vvui *intermittent,
                        unsigned int NumThreads, unsigned int ThreadID );

void CorrelationsTableThread( std::vector< std::vector<bool> > *v,
                              vvui *continuous, vvui *intermittent,
                              unsigned int numHBs,
                              unsigned int fcutoff,
                              unsigned int NumThreads,
                              unsigned int ThreadID);

void getNeighbors( struct Histograms_s *Histograms,
                   std::vector<ListOfHBonds *> HBStrings,
                   struct PBC *Cell);
void
prntHistograms( std::ostream *out,
                std::vector<ListOfHBonds *> HBStrings,
                struct Histograms_s *Histogram,
                std::string CC, unsigned int NumBins,
                struct PBC *Cell, unsigned int TrjIdx,
                bool POVRAY);

template<class T> bool alloc_vector(std::vector< std::vector<T> > *v,
                                    T val,
                                    unsigned int nelem,
                                    unsigned int melem);

template<class T> bool alloc_vector( std::vector<T> *v,
                                     T val,
                                     unsigned int nelem);

#endif // _Histograms_h
