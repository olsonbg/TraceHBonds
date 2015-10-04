/**
 * \file
 * \date  24 April 2015
 * \brief Calls worker thread for jobs in queue
 */
#ifndef _WorkerThreads_h
#define _WorkerThreads_h

#ifdef PTHREADS
#include <vector>
#include <list>
#include <queue>
#include <pthread.h>
#include "VectorTypes.h"
#include "Point.h"
#include "TraceHBonds.h"
#include "ListOfHBonds.h"
#include "Trace.h"
#include "HydrogenBonds.h"
#include "RemoveDuplicates.h"
#include "Histograms.h"
#include "correlation.h"


/**
 * \anchor JobTypes
 * \name   JobTypes
 *
 * Indicate type of job this worker thread should perform
 */
//**@{*/
const unsigned int THREAD_JOB_HBS           =  1; /**< HBs().                     */
const unsigned int THREAD_JOB_RMDUPS        =  2; /**< RemoveDuplicatesThread().  */
const unsigned int THREAD_JOB_TRACE         =  3; /**< TraceThread().             */
const unsigned int THREAD_JOB_CORR          =  4; /**< CorrelationsThread().      */
const unsigned int THREAD_JOB_LIFETIME      =  5; /**< LifetimeThread().          */
const unsigned int THREAD_JOB_POSITIONS_CAR =  6; /**< PositionsCAR().            */
const unsigned int THREAD_JOB_HBS2          =  7; /**< HBs().                     */
const unsigned int THREAD_JOB_CORR_TABLE    =  8; /**< CorrelationsTableThread(). */
const unsigned int THREAD_JOB_SIZEHIST      =  9; /**< makeHistograms().          */
const unsigned int THREAD_JOB_NEIGHBORHIST  = 10; /**< getNeighbors().            */
const unsigned int THREAD_JOB_PAUSE         = 90; /**< Pause thread.              */
const unsigned int THREAD_JOB_EXIT          = 99; /**< Exit thread.               */
//**@}*/

/**
 *
 * Structure to hold data for worker threads to use. The current list of
 * functions for worker threads is:
 *   - HBs()
 *   - RemoveDuplicatesThread()
 *   - TraceThread()
 *   - CorrelationsThread()
 *   - LifetimeThread()
 *   - PositionsCAR()
 *   - HBs()
 *   - CorrelationsTableThread()
 *
 * This list may not be exhaustive.
 *
 * \todo Look into splitting this struct into separate ones for each thread
 * function.
 */
struct worker_data_s
{
	/**
	 * \name
	 * Common to all jobs.
	 */
	//**@{*/
	/** The type of job to run (e.g. THREAD_JOB_HBS)
	 * See also \ref JobTypes
	 */
	unsigned int jobtype;
	//**@}*/

	/**
	 * \name
	 * Common to all jobs. Used for loops which are split between more than one
	 * thread.  e.g.  \c "for(i=jobnum; i < end ; i += num_threads)"
	 */
	//**@{*/
	/** Job number, 0..\p num_threads - 1 */
	unsigned int jobnum;

	/** Total number of threads this jobs has been split into */
	unsigned int num_threads;
	//**@}*/

	/**
	 * \name
	 *
	 * See \p HBs functions at \ref HBfunctions for description of variables.
	 */
	//**@{*/
	std::vector<struct HydrogenBond *> *hb;
	Point cell;
	std::vector<struct thbAtom *>*hydrogens;
	std::vector<struct thbAtom *>*acceptors;
	unsigned int TrjIdx;
	double rCutoff;
	double angleCutoff;
	//**@}*/

	/** \name
	 *
	 * Used for calls to TraceThread(), makeHistograms(), and
	 * RemoveDuplicatesThread(). See their respective functions for a
	 * description of these variables.
	 */
	//**@{*/
	struct HydrogenBondIterator_s *HBit;
	std::vector<ListOfHBonds *> *HBStrings;
	struct Histograms_s *Histogram;

	//**@}*/

	/** \name
	 *
	 * Used for CorrelationsThread(), CorrelationsTableThread() and
	 * LifetimeThread(). See their respective functions for a description of
	 * these variables.
	 */
	//**@{*/
	vvui *vvuiC;
	vvui *vvuiI;
	vd *vdC;
	vd *vdI;
	std::vector< std::vector<bool> > *b;
	std::vector<struct HydrogenBondIterator_s> *TrjIdx_iter;
	unsigned int numHBs;
	unsigned int fcutoff;
	//**@}*/

	/** \name
	 * Use for calls to PositionsCAR().
	 */
	//**@{*/
	const char *filename;
	unsigned int everyNth;
	std::vector<struct thbAtom *> *atom;
	struct PBC *Cell;
	std::vector<Point> *coordinates;
	bool saveMemory;
	//**@}*/
};

#endif // PTHREADS
#endif // _WorkerThreads.h
