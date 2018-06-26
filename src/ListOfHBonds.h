/**
 * \file
 * \author  Brian G. Olson
 * \date    30 July 2009
 * \brief   Class for hydrogen bonds
 */
#ifndef _ListOfHBonds_h
#define _ListOfHBonds_h
#include <map>
#include <vector>
#include <math.h>
#include <string.h>
#include "VectorTypes.h"
#include "OutputFormat.h"
#include "Point.h"
#include "flags.h"

/**
 *
 * Structure to hold parameters for the periodic boundry conditions of all
 * frames(snapshots) in a trajectory.
 *
 * Only periodic boundry conditions with all 90 degree angles are
 * supported.
 *
 **/
struct PBC
{
	std::vector<class Point> p;      /**< x,y, and z dimensions          */
	std::vector<class Point> angles; /**< alpha, beta, and gamma angles  */

	unsigned int frames;             /**< Number of frames in trajectory */

	// Assign default values
	public:
	PBC() { frames = 0; }
};

/**
 * Structure to hold bonds
 **/
struct thbBond
{
	struct thbAtom *A;      /**< Atom participating in this bond */
	struct thbAtom *B;      /**< Atom participating in this bond */
	int Order;              /**< Bond order */
	unsigned int A_ID;      /**< Atom ID participating in this bond (LAMMPS) */
	unsigned int B_ID;      /**< Atom ID participating in this bond (LAMMPS) */
};


/**
 * Structure to hold molecule information
 **/
struct thbMolecule
{
	std::string Name;
	/** All atoms associated in this molecule **/
	std::vector<struct thbAtom *> atoms;
	std::vector<struct thbBond *> bonds;
};

/**
 * Structure to hold atom information
 **/
struct thbAtom
{
	/**
	 * The coordinates for this atom in all frames of a trajectory.
	 */
	std::vector<class Point> p;

	// These should not change with different frames(snapshots).
	unsigned int ID;              /**< Unique ID of this atom                */
	std::string Type;             /**< Type of atom (E.g., "C", "H", "N")    */
	std::string Name;             /**< Unique name of atom                   */
	std::string Residue;          /**< Residue this atom belongs to          */
	unsigned int ResidueNum;      /**< Residue number which this atom belongs*/
	struct thbMolecule *Molecule; /**< Molecule to which this atom belongsi  */
	std::string ForceField;       /**< Forcefield used for this atom         */
	unsigned int AtomTypeID;      /**< Atom type, mass #, from LAMMPS */
	/**
	 * How many times this atom can contribute to hydrogen bonding
	 **/
	unsigned int HydrogenBondMax;
	/** Bonds this atom participates in */
	std::vector<struct thbBond *> Bonds;

	/**
	 * \name Atoms connected to this one.
	 **/
	/**@{*/

	/** All atoms connected to this one */
	std::vector<struct thbAtom *> ConnectedAtom;
	/** Name of all atoms connected to this one */
	std::vector<std::string>  ConnectedAtomName;
	/** Residue of all atom connected to this one */
	std::vector<std::string>  ConnectedAtomResidue;
	/** Residue number of all atom connected to this one */
	std::vector<unsigned int> ConnectedAtomResidueNum;
	/** Bond order of all atom connected to this one */
	std::vector<float>        ConnectedAtomBondOrder;
	/**@}*/

	// Assign default values
	public:
	thbAtom() { HydrogenBondMax = 0; }
};

/**
 * Struct containing hydrogen bond parameters
 **/
struct HydrogenBond
{
	/**
	 * Distance between Hydrogen and acceptor atom.
	 **/
	double length;
	/**
	 * Distance between acceptor and donor atoms
	 **/
	double acceptorDonorDistance;
	/**
	 * Angle formed between two vectors, the first pointing from the hydrogen
	 * to the acceptor, and the second pointing from the hydrogen to the atom
	 * covalently bonded to it. See also Point::angle.
	 */
	double angle;

	struct thbAtom *hydrogen; /**< Hydrogen atom */
	struct thbAtom *donor;    /**< Donor atom    */
	struct thbAtom *acceptor; /**< Acceptor atom */

	/**
	 * Frame of trajectory this hydrogen bond exists (Trajectory Index)
	 *
	 * \todo I use 'frame' and 'trajectory index' interchangeably, I should
	 * pick one and stick with it.
	 *
	 **/
	unsigned int TrajIdx;

	/**
	 * \name Linked list of hydrogen bonds
	 */
	/**@{*/
	struct HydrogenBond *Next;      /**< Next hydrogen bond in list     */
	struct HydrogenBond *Previous;  /**< Previous hydrogen bond in list */
	/**@}*/

	/**
	 * Used for housekeeping later.
	 * \todo Elaborate on this.
	 **/
	bool markedDuplicate; /**< Used to mark hydrogen bond as duplicate *  for housecleaning later */

	// Assign default values
	public:
	HydrogenBond()
	{
		Next = Previous = NULL;
		hydrogen = NULL;
		donor = NULL;
		acceptor = NULL;
		length = angle = acceptorDonorDistance = 0.0;
		markedDuplicate = false;
	}
};

/**
 * Maintains list of hydrogen bonds which are connected to each other.
 *
 * Builds strings of hydrogen bonds which are connected to each other.
 *
 **/
class ListOfHBonds
{
	public:
		// List Builder
		ListOfHBonds();
		/**
		 * Number of atoms in this list (string) of hydrogen bonds
		 *
		 * \returns Number of atoms
		 */
		unsigned int AtomCount();
		/**
		 * Number of hydrogen bonds in this list (string) of hydrogen bonds
		 *
		 * \return Number of hydrogen bonds
		 **/
		unsigned int HydrogenBondCount();
		// List Operations
		/**
		 * \return Pointer to first hydrogen bond in this list
		 **/
		struct HydrogenBond *Begin();
		/**
		 * \return Pointer to last hydrogen bond in this list
		 **/
		struct HydrogenBond *End();
		/**
		 *
		 * Adds \p Item to end of list of hydrogen bonds.
		 *
		 * \param[in] Item  Pointer to HydrogenBond
		 *
		 * \return Number of atoms in the new list
		 *
		 **/
		int AddAtEnd(struct HydrogenBond *Item);
		/**
		 *
		 * Adds \p Item to beginning of list of hydrogen bonds.
		 *
		 * \param[in] Item  Pointer to HydrogenBond
		 *
		 * \return Number of atoms in the new list
		 *
		 **/
		int AddAtBegin(struct HydrogenBond *Item);
		/**
		 *
		 * \param[in] pos  Index of hydrogen bond in list.
		 *
		 * \return Pointer to HydrogenBond
		 *
		 **/
		struct HydrogenBond *Retrive(int pos);
		/**
		 * Delete this list of hydrogen bonds
		 **/
		bool DeleteList();
		/**
		 *
		 * \param[in] Item Pointer to HydrogenBond
		 *
		 * \return \c TRUE if \p Item exists in this list of hydrogen bonds.
		 *
		 **/
		bool Find( struct HydrogenBond *Item);
		/**
		 * Frame number (Trajectory index) this list of hydrogen bonds exists
		 * in.
		 *
		 * \return Trajectory index number.
		 *
		 **/
		unsigned int TrajectoryIndex();
		/**
		 *
		 * \return Number of times this list of hydrogen bonds switches between
		 * different molecules.
		 *
		 **/
		unsigned int SwitchingCount();
		/**
		 *
		 * \return Number of unique atom forcefields in this list of hydrogen
		 * bonds
		 *
		 **/
		unsigned int ForcefieldCount();
		/**
		 *
		 * \return Number of unique molecules in this list of hydrogen bonds
		 *
		 **/
		unsigned int MoleculeCount();
		/**
		 *
		 * Counts number of unique strings. Used by MoleculeCount() and
		 * ForcefieldCount().
		 *
		 * \param[in] s  Vector of strings
		 * \return Number of unique strings in \p s
		 *
		 **/
		unsigned int CountUniqStr( std::vector< std::string >s );
		/**
		 *
		 * \return Vector of all donor-acceptor distances in this list of
		 * hydrogen bonds.
		 *
		 **/
		vd donorAcceptorDistances();
		/**
		 *
		 * \return Coordinates of all atoms in this list of hydrogen bonds,
		 * expect for the hydrogens.
		 *
		 **/
		std::vector< Point *>nonHydrogenCoordinates();
		/**
		 *
		 * If the acceptor of \p Item is the same atom as the donor of the
		 * first element in the list, then Item should come before Begin().
		 *
		 * \param[in] Item Pointer to a HydrogenBond.
		 *
		 * \return \c TRUE if \p Item should become first hydrogen bond in
		 * list. \c FALSE otherwise.
		 *
		 **/
		bool linksAtBegin( struct HydrogenBond *Item);
		/**
		 * If the donor of \p Item is the same atom as the acceptor of the last
		 * element in the list, then Item should come after End().
		 *
		 * \param[in] Item Pointer to a HydrogenBond.
		 *
		 * \return \c TRUE if \p Item should become last hydrogen bond in list.
		 * \c FALSE otherwise.
		 *
		 **/
		bool linksAtEnd  ( struct HydrogenBond *Item);
		/**
		 *
		 * \return \c TRUE if this list is a closed loop, \c FALSE otherwise
		 *
		 **/
		bool ClosedLoop(void);
		/**
		 *
		 * Calculates the minimum image vector of atom \p A, in frame \p
		 * TrjIdx, from coordinate \p r.
		 *
		 * \param[in] A       Pointer to atom
		 * \param[in] TrjIdx  Trajectory index (frame) of interest
		 * \param[in] r       Coordinate of point of interest
		 * \param[in] Cell    Parameters for periodic boundary conditions
		 *
		 **/
		Point MinimumImage( struct thbAtom *A,
		                    unsigned int TrjIdx,
		                    Point r,
		                    struct PBC Cell);

		/**
		 *
		 * Print out this list of hydrogen bonds.
		 *
		 * \param[in]  out     Stream to send output to
		 * \param[in]  Cell    Parameters for periodic boundary conditions
		 * \param[in]  TrjIdx  Trajectory index (frame) of interest
		 * \param[in]  POVRAY  Bool indicating whether to output in PovRay
		 * format, or plain text.
		 *
		 * \return End-to-end distance of the list of hydrogen bonds
		 *
		 **/
		double PrintAll(std::ostream *out,
		                struct PBC Cell,
		                unsigned int TrjIdx,
		                unsigned int flags);

	private:
		unsigned int size;
		unsigned int HBsize;
		struct HydrogenBond *Beginning;
		Point  Round (Point p, double f=1.0);
		double Round (double r, double f=1.0);

};
#endif // _ListOfHBonds_h
