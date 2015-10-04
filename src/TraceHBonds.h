/**
 * \file
 * \brief   Trace hydrogen bonded atom to find all connections
 * \author  Brian G. Olson (olsonbg@gmail.com)
 * \date    30 July 2009
 */
#ifndef _TraceHBonds_h
#define _TraceHBonds_h
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <unistd.h>
#include <math.h>
#include <string.h>
#include <vector>
#include <iostream>
#include <set>
#include "MessageDefines.h"
#include "ReadCarMdf.h"
#include "ListOfHBonds.h"
#include "Lifetime.h"
#include "Point.h"

extern bool THB_VERBOSE;

/**
 * Flags for command line options
 */
enum {
	VERBOSE       = 0x01, /**< Be verbose                         */
	LIFETIME      = 0x02, /**< Calculate lifetime correlations    */
	SIZE_HIST     = 0x04, /**< Generate size histograms           */
	LENGTHS       = 0x08, /**< Generate hydrogen bond length list */
	POVRAY        = 0x10, /**< Ouput in PovRay format             */
	NEIGHBOR_HIST = 0x20, /**< Generge neighbor histograms        */
	ANGLES        = 0x40, /**< Generate hydrogen bond angle list  */
	/** Do all possible calculations */
	ALL           = (SIZE_HIST|LENGTHS|NEIGHBOR_HIST|LIFETIME|ANGLES),
};

/**
 * Atoms which contribute to hydrogen bonding, specified by their forcefields,
 * thbAtom::ForceField.
 */
struct HydrogenBondMatching
{
	/**
	 * Forcefield (thbAtom::ForceField) of hydrogens that may contribute to
	 * hydrogen bonding
	 */
	std::vector<std::string>Hydrogens;
	/**
	 * Forcefield  (thbAtom::ForceField) of atoms which may act as hydrogen
	 * bond acceptors
	 */
	std::vector<std::string>Acceptors;
};

/**
 * Iterators which point to the first (begin), and just past the last (end),
 * elements in a series of HydrogenBond structures. Each begin and end pair
 * correspond to a single frame in a trajectory.
 */
struct HydrogenBondIterator_s
{
	/** Iterator which points to the first element */
	std::vector<struct HydrogenBond *>::iterator begin;
	/** Iterator which point to just past the last element */
	std::vector<struct HydrogenBond *>::iterator end;
};

/**
 *
 * Main function to drive calculations.
 *
 * \param[in] ifilename    Name of file to open
 * \param[in] ofPrefix     Prefix of filename to save data (e.g. "HBonds")
 * \param[in] ofSuffix     Suffix of filename to save data (e.g. ".dat")
 * \param[in] match        Atoms which contribute to hydrogen bonding, specified by their forcefields,
 * \param[in] rCutoff      Maximum distance to consider for hydrogen bonding
 * \param[in] angleCutoff  Minimum angle to consider for hydrogen bonding
 * \param[in] NumBins      Number of bins to output
 * \param[in] EveryNth     Load 'EveryN' frame from trajectory
 * \param[in] flags        User specified flags from command line
 *
 * \return integer \c 0 on success, \c 1 on error encountered
 */
int doArcFile(char *ifilename,
              char *ofPrefix, char *ofSuffix,
              struct HydrogenBondMatching *match,
              double rCutoff, double angleCutoff,
              unsigned int NumBins, unsigned int EveryNth,
              unsigned char flags);



/**
 * Used with std::remove_if() from DeleteVectorPointers()
 *
 * delete \p v
 *
 * \param v  Pointer to an element
 *
 * \return \c TRUE
 */
template <typename T> bool deleteVectorPointers( T* v );

/**
 * Delete and then remove each element in vector of pointers. Uses
 * std::remove_if() with deleteVectorPointers().
 *
 * \param v  Vector of pointers
 */
template<class T> void DeleteVectorPointers( std::vector<T*> v );


/**
 * Generate Iterators which point to hydrogen bonds in specific frames.
 * TrjIdx_iter.at(5).begin points to first hydrogen bond (\p hb) in frame 5
 * (frame numbers start at 0).  TrjIdx_iter.at(5).end points to just past the
 * last hydrogen bond in frame 5.  The hydrogen bonds in \p hb are grouped by
 * trajectory index number, however the order of the groups may not be in
 * sequence because they are obtained from threads.
 *
 * \param[in,out] TrjIdx_iter  Iterators for begin and end of \p hb in each
 *                             frame
 * \param[in]     hb           All hydrogem bonds
 */
void
TrajectoryIndexIterator( std::vector< struct HydrogenBondIterator_s > *TrjIdx_iter,
                         std::vector< struct HydrogenBond *> *hb);
#endif
