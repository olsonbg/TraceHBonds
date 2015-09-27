/**
 * \file
 * \brief Remove hydrogen bonds which have multiple bonding
 *
 * Possible to have more than one:
 *
 *   - acceptor  with a single hydrogen
 *   - Hydrogen with a single acceptor
 *
 * Find the all duplicates and keep the one with the shortest Hydrogen Bond
 * length.
 *
 * \todo Make this work for atoms which may hydrogen bond more than once.
 */
#ifndef _RemoveDuplicates_h
#define _RemoveDuplicates_h

#include "TraceHBonds.h"
#include "queue.h"
#include "WorkerThreads.h"
#include "cpu.h"

extern bool THB_VERBOSE;

#ifdef PTHREADS
extern Queue<struct worker_data_s> inQueue;
extern Queue<struct worker_data_s> outQueue;
#endif

/** Vector of pointers to HydrogenBond struct */
typedef std::vector<struct HydrogenBond *> HBVec;
/** Vector of HydrogenBond struct Iterators */
typedef std::vector<struct HydrogenBondIterator_s> HBVecIter;

/**
 * Structure holding function for compare bool operator of the  std::sort()
 * function.
 */
struct sort_Neighbors {
	/**
	 * Sort by length of hydrogen bond, the shorter comes first.
	 *
	 * \param[in] left   left argument of std::sort()
	 * \param[in] right  right argument of std::sort()
	 *
	 * \return \c TRUE if hydrogen bond length of \c left is less than \c
	 * right, \c FALSE otherwise
	 */
	bool operator()(const HBVec::iterator &left, const HBVec::iterator &right)
	{
		return (*left)->length < (*right)->length;
	}
};

/**
 * Unary predicate used in std::remove_if().
 *
 * \param[in] hb   Pointer to a HydrogenBond structure
 *
 * \return \c TRUE if \p hb has been marked as a duplicate, \c FALSE otherwise
 */
bool deleteMarked( struct HydrogenBond *hb );

/**
 * Remove all hydrogen bonds which have been marked as duplicate.
 *
 * \param[in] hb   List of hydrogen bonds
 */
void removeMarked( HBVec *hb );
/**
 * Find the all duplicates and keep the shortest Hydrogen Bond length.
 *
 * Put NumberOfCPUs() RemoveDupliateThread() jobs in Queue.
 *
 * \param[in] hb           List of hydrogen bonds
 * \param[in] TrjIdx_iter  Indicates first and last objects of \p hb for each
 *                         frame
 */
void RemoveDuplicates( std::vector<struct HydrogenBond *> *hb,
                       std::vector< struct HydrogenBondIterator_s > *TrjIdx_iter);

/**
 * Find the all duplicates in single frame, and keep the shortest Hydrogen Bond
 * length.
 *
 * Called as a thread from RemoveDupicates().
 *
 * \param[in] HBit  Indicates first and last objects in \p hb for this
 *                  frame
 */
void RemoveDuplicatesThread( struct HydrogenBondIterator_s HBit );

/**
 *
 * \param[in] A   Pointer to an atom
 * \param[in] B   Pointer to an atom
 *
 * \return \c TRUE if \p A and \p B point to the same object.
 */
bool SameAtom( struct thbAtom *A,
               struct thbAtom *B);

#endif // _RemoveDuplicates_h

