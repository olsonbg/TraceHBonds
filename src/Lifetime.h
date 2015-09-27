#ifndef _Lifetime_h
#define _Lifetime_h
/**
 * \file
 * \author Brian G. Olson
 * \date   September 24th, 2015
 * \brief  Track bonded and non-bonded states of each hydrogen bond
 */
#include "TraceHBonds.h"
#include "queue.h"
#include "WorkerThreads.h"
#include "cpu.h"

extern bool THB_VERBOSE;

/** Vector of pointers to HydrogenBond struct */
typedef std::vector<struct HydrogenBond *> HBVec;
/** Vector of HydrogenBond struct Iterators. */
typedef std::vector<struct HydrogenBondIterator_s> HBVecIter;

/**
 * Track bonded and non-bonded states of each hydrogen bond
 *
 * Puts LifetimeThread() in the Queue, and combines results from each job.
 *
 * \param[out] b            Boolean indicating wheter hydrogen bonds are formed,
 *                          and in which frame
 * \param[in]  TrjIdx_iter  Structure containing iterator for beginning and end
 *                          of each frame.
 */
void
Lifetime( std::vector< std::vector<bool> >*b,
          std::vector<struct HydrogenBondIterator_s > *TrjIdx_iter );

/**
 * Determine which hydrogen bonds are formed in each frame, a value of \c TRUE
 * means the hydrogen bond is bonded in that frame. A value of \c FALSE means
 * that the atoms comprising a potential hydrogen bond are NOT bonded in that
 * frame.
 *
 * This function is called as a thread, put in the Queue by Lifetime().
 *
 * \param[out] b            Boolean indicating whether hydrogen bonds are formed,
 *                          and in which frame
 * \param[in]  TrjIdx_iter  Structure containing iterator for beginning and end
 *                          of a frame.
 * \param[in]  NumThreads   Number of jobs this calculation has been split into
 *                          (usually NumberOfCPUs())
 * \param[in]  ThreadID     Job ID of this call <tt>(0 <= ThreadID < NumThreads)</tt>
 *
 */
void
LifetimeThread(std::vector< std::vector<bool> >*b,  HBVecIter *TrjIdx_iter,
               unsigned int NumThreads=1, unsigned int ThreadID=0);
#endif // _Lifetime_h

