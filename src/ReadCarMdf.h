/**
 * \file
 * \date  20 April 2015
 * \brief Read Discover car/mdf file pairs
 */
#ifndef _ReadCarMdf_h
#define _ReadCarMdf_h
#include <boost/iostreams/filtering_stream.hpp>
// #include <boost/iostreams/device/file.hpp>
// #include <boost/iostreams/copy.hpp>
// #include <boost/iostreams/filter/gzip.hpp>
// #include <boost/iostreams/filter/bzip2.hpp>
#include <vector>
#include <math.h>
#include "OpenFile.h"
#include "queue.h"
#include "WorkerThreads.h"
#include "cpu.h"
#include "ListOfHBonds.h"
#include "TraceHBonds.h"

/**
 * Read Discover CAR or ARC file.
 *
 * Reads atoms, coordinates, forcefields, and periodic boundary conditions from
 * CAR or a single frame of ARC file. Consecutive calls to this function will
 * read next frame.
 *
 * \param[in]  in          boost filtering_stream previously setup for file
 * \param[out] atom        atom vector containing description of all atoms
 * \param[out] Cell        periodic boundary conditions read from file
 * \param[out] Coordinates Coordinates of atoms in a frame.
 * \param[in]  SaveMemory  Bool stating memory saving routines can be used.
 *
 * \return \c TRUE on success, \c FALSE otherwise
 */
bool ReadCar(boost::iostreams::filtering_stream<boost::iostreams::input> *in,
             std::vector<struct thbAtom *> *atom,
             struct PBC *Cell, std::vector<Point> *Coordinates,
             bool SaveMemory);

/**
 * Read Discover MDF file.
 *
 * Reads atoms connections, and molecule and residue information from MDF file.
 *
 * \param[in]     in      boost filtering_stream previously setup for file
 * \param[in,out] atom    atom vector containing description of all atoms
 * \param[in,out] atom    vector containing description of all molecule
 *
 * \return \c TRUE on success, \c FALSE otherwise
 */
bool ReadMdf(boost::iostreams::filtering_stream<boost::iostreams::input> *in,
             std::vector<struct thbAtom *> *atom,
             std::vector<struct thbMolecule *> *molecule);

/**
 * Read connections from MDF file
 *
 * Calls openfile(), then ReadMdf()
 *
 * \param[in]     filename   Name of file to open
 * \param[in,out] atom       vector containing description of all atoms
 * \param[in,out] atom       vector containing description of all molecule
 *
 * \return \c TRUE on success, \c FALSE othewise
 */
bool ConnectionsMDF(const char *filename,
                    std::vector<struct thbAtom *> *atom,
                    std::vector<struct thbMolecule *> *molecules,
                    std::vector<struct thbBond *> *bonds);
/**
 * CAR/ARC
 *   - calls openfile()
     - Loops over all frames
 *      - calls ReadCar() to get next frame
 *      - For this frame, makes puts the \ref SaveMemory HBs function in Queue,
 *        otherwise the \ref standard HBs function is queued.
 *
 * \param[in] filename        Name of CAR/ARC file
 * \param[out] atom           atom vector containing description of all atoms
 * \param[out] Cell           periodic boundary conditions read from file
 * \param[in]  hydrogens      Vector containing all hydrogens that may hydrogen
 *                            bond
 * \param[in]  acceptors      Vector containing all atoms that may act as
 *                            hydrogen bond acceptors
 * \param[in]  rCutoff        Maximum distance for hydrogen bonding
 * \param[in]  angleCutoff    Minimum angle for hydrogen bonding
 * \param[in]  SaveMemory     Bool stating memory saving routines can be used
 */
bool PositionsCAR(const char *filename,
                  std::vector<struct thbAtom *> *atom,
                  struct PBC *Cell,
                  std::vector<struct thbAtom *> *hydrogens,
                  std::vector<struct thbAtom *> *acceptors,
                  double rCutoff, double angleCutoff, bool SaveMemory);

#endif

