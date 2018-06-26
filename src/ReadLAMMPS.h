#ifndef _ReadLAMMPS_h
#define _ReadLAMMPS_h
#include <vector>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include "queue.h"
#include "WorkerThreads.h"
#include "cpu.h"
#include "ListOfHBonds.h"
#include "TraceHBonds.h"
#ifdef USE_LZMA
#include "lzma.h"
#endif

struct AtomTypes
{
	unsigned int id;
	double mass;
	char comment[1024];
};

/**
 * Molecule definitions which show the atoms contained in a molecule of id. The
 * atoms must be consecutively numbered. The atom numbers correspond to the
 * numbers in the LAMMPS data file.
 */
struct MoleculeDefs
{
	/**
	 * ID of molecule
	 */
	unsigned int id;
	/**
	 * First atom that belongs to this molecule
	 */
	unsigned int atomFirst;
	/**
	 * Last atom that belongs to this molecule
	 */
	unsigned int atomLast;
};

/**
 * Read Molecule definition file.
 *
 * Reads molecule ID, first atom, and last atom for this molecule.
 *
 * \param[in]  in          boost filtering_stream previously setup for file
 * \param[in,out] moldefs  MolecularDefs
 *
 * \return \c TRUE on success, \c FALSE otherwise
 */
bool ReadLAMMPSMolDefs(boost::iostreams::filtering_stream<boost::iostreams::input> *in,
                       std::vector<struct MoleculeDefs *> *moldefs);

/**
 * Read number of atoms and bonds, along with the bond definitions from LAMMPS
 * datafile.
 */
bool ReadLAMMPSData(boost::iostreams::filtering_stream<boost::iostreams::input> *in,
                    std::vector<struct thbAtom     *> *atom,
                    std::vector<struct thbBond     *> *bonds,
                    std::vector<struct MoleculeDefs *> moldefs );

/**
 * Read LAMMPS datafile and trajectory file, along with a molecule definition
 * file.
 */
bool ReadLAMMPSConnections(char *fileData,
                           char *fileMols,
                           std::vector<struct thbAtom     *> *atom,
                           std::vector<struct thbMolecule *> *molecules,
                           std::vector<struct thbBond     *> *bonds );

int doLAMMPSFile(char *fileData, char *fileTrj, char *fileMols,
                 char *ofPrefix, char *ofSuffix,
                 struct HydrogenBondMatching *match,
                 double rCutoff, double angleCutoff,
                 int NumBins, unsigned int flags);

bool ReadLAMMPSPositions(const char *fileTrj,
                         std::vector<struct thbAtom *> *atom,
                         struct PBC *Cell,
                         std::vector<struct thbAtom *> *hydrogens,
                         std::vector<struct thbAtom *> *acceptors,
                         double rCutoff, double angleCutoff, bool SaveMemory);
#endif
