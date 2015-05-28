#ifndef _HydrogenBonds_h
#define _HydrogenBonds_h

#include "TraceHBonds.h"
#include "queue.h"
#include "WorkerThreads.h"
#include "cpu.h"

void getHydrogenBondElements( std::vector<struct thbAtom *> *atom,
                              std::vector<struct thbAtom *> *hydrogendonors,
                              std::vector<struct thbAtom *> *acceptors,
                              struct HydrogenBondMatching *match);

void HBs( std::vector<struct HydrogenBond *> *hb,
          Point cell,
          std::vector<struct thbAtom *>*hydrogens,
          std::vector<struct thbAtom *>*acceptors,
          double TrjIdx, double rCutoff, double angleCutoff);

void AtomNeighbors( std::vector<struct HydrogenBond *> *hb,
                    std::vector<struct thbAtom *> *atom,
                    struct PBC *Cell, struct HydrogenBondMatching *match,
                    double rCutoff, double angleCutoff );
#endif // _HydrogenBonds_h