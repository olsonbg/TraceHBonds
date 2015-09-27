/**
 * \file
 * \author Brian G. Olson
 * \date   23 May 2015
 * \brief Get hydrogen bonds.
 */
#ifndef _HydrogenBonds_h
#define _HydrogenBonds_h

#include "TraceHBonds.h"
#include "queue.h"
#include "WorkerThreads.h"
#include "cpu.h"

extern bool THB_VERBOSE;
/*
#ifdef PTHREADS
extern Queue<struct worker_data_s> inQueue;
extern Queue<struct worker_data_s> outQueue;
#endif
*/
/** Vector of pointers to HydrogenBond struct */
typedef std::vector<struct HydrogenBond *> HBVec;
/** Vector of HydrogenBond struct Iterators. */
typedef std::vector<struct HydrogenBondIterator_s> HBVecIter;

/**
 * Get all hydrogens and acceptors that can participate in hydrogen bonding
 *
 * The forcefield of matching atoms may be specified once for each number of times
 * they may participate in hydrogen bonding.
 *
 * \param[in]  atom           Vector containing all atoms in the system
 * \param[out] hydrogendonors Vector containing all hydrogens that may hydrogen
 *                            bond
 * \param[out] acceptors      Vector containing all atoms that may act as
 *                            hydrogen bond acceptors.
 * \param[in]  match          Forcefields of hydrogens and acceptors for hydrogen
 *                            bonding.
 *
 */
void getHydrogenBondElements( std::vector<struct thbAtom *> *atom,
                              std::vector<struct thbAtom *> *hydrogendonors,
                              std::vector<struct thbAtom *> *acceptors,
                              struct HydrogenBondMatching *match);
/**
 * \anchor HBfunctions
 * \name Find all hydrogen bonds
 *
 * One standard version, and another which can save memory for certain user
 * selected options.
 */
/**@{*/
/**
 * \anchor SaveMemory
 * SaveMemory version of routine to find all hydrogen bonds
 *
 * If neither NEIGHBOR_HIST or SIZE_HIST are specified, we can save some memory
 * by storing the atom coordinates only as long as we need them, otherwise, use
 * \ref standard HBs function.
 *
 * This function is called as a thread, put in the Queue by PositionsCAR().
 *
 * \todo Add links to NEIGHBOR_HIST and SIZE_HIST
 * \todo Why is TrjIdx a double?
 *
 * \param[out] hb         Vector of all hydrogen bonds
 * \param[in]  cell       Dimension of periodic cell
 * \param[in]  hydrogens  Vector of all hydrogens that may participate in
 *                        hydrogen bonding
 * \param[in]  acceptors  Vector of all atoms that may act as hydrogen bond
 *                        acceptors
 * \param[in]  Coordinates  X
 * \param[in]  TrjIdx       Trajectory index to operate on (Frame number)
 * \param[in]  rCutoff      Maximum distance to consider for hydrogen bonding
 * \param[in]  angleCutoff  Minimum angle to consider for hydrogen bonding
 *
 */
void HBs( std::vector<struct HydrogenBond *> *hb,
          Point cell,
          std::vector<struct thbAtom *>*hydrogens,
          std::vector<struct thbAtom *>*acceptors,
          std::vector<Point> *Coordinates,
          double TrjIdx,
          double rCutoff, double angleCutoff);
/**
 * \anchor standard
 *
 * Standard version of routine to find all hydrogen bonds
 *
 * This function is called as a thread, put in the Queue by PositionsCAR().
 *
 * \todo Add links to NEIGHBOR_HIST and SIZE_HIST
 * \todo Why is TrjIdx a double?
 *
 * \param[out] hb           Vector of all hydrogen bonds
 * \param[in]  cell         Dimension of periodic cell
 * \param[in]  hydrogens    Vector of all hydrogens that may participate in
 *                          hydrogen bonding
 * \param[in]  acceptors    Vector of all atoms that may act as hydrogen bond
 *                          acceptors
 * \param[in]  TrjIdx       Trajectory index to operate on (Frame number)
 * \param[in]  rCutoff      Maximum distance to consider for hydrogen bonding
 * \param[in]  angleCutoff  Minimum angle to consider for hydrogen bonding
 *
 */
void HBs( std::vector<struct HydrogenBond *> *hb,
          Point cell,
          std::vector<struct thbAtom *>*hydrogens,
          std::vector<struct thbAtom *>*acceptors,
          double TrjIdx, double rCutoff, double angleCutoff);
/**@}*/

#endif // _HydrogenBonds_h
