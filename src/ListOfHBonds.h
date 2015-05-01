#ifndef _ListOfHBonds_h
#define _ListOfHBonds_h
#include <map>
#include <vector>
#include <math.h>
#include <string.h>
#include "OutputFormat.h"
#include "Point.h"

struct PBC
{
	// Store the PBC parameters for many frames(snapshots) in a trajectory.
	std::vector<class Point> p;
	std::vector<class Point> angles; // alpha, beta, gamma;
	unsigned int frames;
	// 
	// Assign default values
	public:
	PBC()
	{
		frames = 0;
	}
};

struct thbAtom
{
	// Store the coordinates for many frames(snapshots) in a trajectory.
	std::vector<class Point> p;

	// These should not change with different frames(snapshots).
	std::string Type; // E.g., C, H, N, Au.
	std::string Name;
	std::string Residue;
	unsigned int ResidueNum;
	std::string Molecule;
	std::string ForceField;
	// float BondOrder;
	// Atoms connected to this one.
	std::vector<struct thbAtom *> ConnectedAtom;
	std::vector<std::string>  ConnectedAtomName;
	std::vector<std::string>  ConnectedAtomResidue;
	std::vector<unsigned int> ConnectedAtomResidueNum;
	std::vector<float>        ConnectedAtomBondOrder;

	// Assign default values
	public:
	thbAtom()
	{
		// Hydrogen = false;
		// x = y = z = 0.0;
		// length = angle = 0.0;
		// markedDuplicate = false;
	}
};

//struct HBondAtom
struct HydrogenBond
{
	double length;
	double acceptorDonorDistance;
	double angle;
	struct thbAtom *hydrogen;
	struct thbAtom *donor;
	struct thbAtom *acceptor;
	// Index of the coordinates in thbAtom for the hydrogen, donor, and
	// acceptor atoms. All three are at the same index. The same index is used
	// for the coordinates and angles in the PBC struct too.
	unsigned int TrajIdx;

	// char dendrimer[80];
	struct HydrogenBond *Next;
	struct HydrogenBond *Previous;
	// bool Hydrogen;
	// char Name[80];
	bool markedDuplicate;

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

class ListOfHBonds
{
	private:
		unsigned int size;
		unsigned int HBsize;
		struct HydrogenBond *Beginning;
		double Round (double r, double f);

	public:
		// List Builder
		ListOfHBonds();
		unsigned int AtomCount();
		unsigned int HydrogenBondCount();
		// List Operations
		struct HydrogenBond *Begin();
		struct HydrogenBond *End();
		int AddAtEnd(struct HydrogenBond *Item);
		int AddAtBegin(struct HydrogenBond *Item);
		struct HydrogenBond *Retrive(int pos);
		// struct HydrogenBond *Last();
		// struct HydrogenBond *First();
		bool DeleteList();
		bool Find( struct HydrogenBond *Item);
		unsigned int TrajectoryIndex();
		unsigned int SwitchingCount();
		unsigned int ForcefieldCount();
		unsigned int MoleculeCount();
		unsigned int CountUniqStr( std::vector< std::string >s );
		std::vector< double >donorAcceptorDistances();
		std::vector< Point *>nonHydrogenCoordinates();
		bool linksAtBegin( struct HydrogenBond *Item);
		bool linksAtEnd  ( struct HydrogenBond *Item);
		bool ClosedLoop(void);
		Point MinimumImage( struct thbAtom *A,
		                    unsigned int TrjIdx,
		                    Point r,
		                    struct PBC Cell);
#ifdef DEBUG
		void PrintAll(unsigned int TrjIdx);
		void PrintAllPovRay(unsigned int TrjIdx);
#endif
		double PrintAll(std::ostream *out,
		                struct PBC Cell,
		                unsigned int TrjIdx,
		                bool POVRAY);
};
#endif // _ListOfHBonds_h
