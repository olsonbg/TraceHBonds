/**
 * \file
 * \author Brian G. Olson
 * \date   30 April 2015
 * \brief  lifetime correlations
 *
 **/
#ifndef _lifetime_h
#define _lifetime_h

#include "ListOfHBonds.h"
#include "VectorTypes.h"
#include "cpu.h"
#include "queue.h"
#include "WorkerThreads.h"

/**
 * Calculate autocorrelation.
 *
 * Puts jobs for both CorrelationsTableThread() and CorrelationsThread() in the
 * Queue, and combines results from all jobs.
 *
 * \param[in] out   Stream to send results to.
 * \param[in] v     Boolean indicating which hydrogen bonds are formed, and
 *                  in which frame
 */
void Correlations( std::ostream *out,
                   std::vector< std::vector<bool> > *v );

/**
 * Combine autocorrelation results for all hydrogen bonds.
 *
 * For both continuous and intermittent hydrogen bond autocorrelations, sum
 * results for all hydrogen bonds in a frame (time slice), then divide by
 * number total number.
 *
 * \param[out] C            Average continuous hydrogen bond autocorrelation
 * \param[out] I            Average intermittent hydrogen bond autocorrelation
 * \param[in]  continuous   Continuous autocorrelation of each hydrogen bond,
 *                          from CorrelationsTableThread()
 * \param[in]  intermittent Intemittent autocorrelation of each hydrogen bond,
 *                          from CorrelationsTableThread()
 * \param[in]  NumThreads   Number of jobs this calculation has been split
 *                          into, determined by NumberOfCPUs().
 * \param[in]  ThreadID     Job ID of this call 0 <=\p ThreadID <\p NumThreads.
 *
 */
void CorrelationsThread(vd *C, vd *I,
                        vvui *continuous, vvui *intermittent,
                        unsigned int NumThreads, unsigned int ThreadID );

/**
 * Calculate autocorrelation of hydrogen bonds.
 *
 * Calculate both the continuous and intermittent hydrogen bond
 * autocorrelations for each hydrogen bond in the system. Uses a sliding
 * window, therefore \p fcutoff should be half the total number of frames.
 *
 * The average over all hydrogen bonds in the system can be obtained by
 * subsequently calling CorrelationsThread().
 *
 * \param[in]  v            Boolean indicating which hydrogen bonds are formed,
 *                          and in which frame
 * \param[out] continuous   Continuous hydrogen bond autocorrelation
 * \param[out] intermittent Intermittent hydrogen bond autocorrelation
 * \param[in]  numHBs       Number of hydrogen bonds.
 * \param[in]  fcutoff      Time cutoff (in number of frame units).
 * \param[in]  NumThreads   Number of jobs this calculation has been split
 *                          into, determined by NumberOfCPUs().
 * \param[in]  ThreadID     Job ID of this call 0 <=\p ThreadID <\p NumThreads.
 *
 */
void CorrelationsTableThread( std::vector< std::vector<bool> > *v,
                              vvui *continuous, vvui *intermittent,
                              unsigned int numHBs,
                              unsigned int fcutoff,
                              unsigned int NumThreads,
                              unsigned int ThreadID);

#endif // _lifetime_h
