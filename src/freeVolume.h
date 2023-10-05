/**
 * \file
 * \author Brian G. Olson
 * \date   17 September 2023
 * \brief  free volume determination
 *
 **/
#ifndef _freevolume_h
#define _freevolume_h

#include "ListOfHBonds.h"
#include "VectorTypes.h"
#include "cpu.h"
#include "queue.h"
#include "WorkerThreads.h"

/**
 * Place holder for now
 */
void freeVolume( std::vector<thbAtom *> *atom,
                 Point cell,
                 unsigned int TrjIdx);

#endif // _freevolume_h
