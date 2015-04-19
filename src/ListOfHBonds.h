#ifndef _ListOfHBonds_h
#define _ListOfHBonds_h
#include <map>
#include <vector>
#include <math.h>
#include <string.h>
#include "OutputFormat.h"

struct PBC
{
	double x, y, z;
	double alpha, beta, gamma;
};

struct thbAtom
{
	// Store the coordinates for many frames(snapshots)
	// std::vector<double> x, y, z;
	double x, y, z;

	// These should not change with different frames(snapshots).
	std::string Type; // E.g., C, H, N, Au.
	std::string Name;
	std::string Residue;
	unsigned int ResidueNum;
	std::string Molecule;
	std::string ForceField;
	float BondOrder;
	// Atoms connected to this one.
	std::vector<struct thbAtom *> Connected;

	// Assign default values
	public:
	thbAtom()
	{
		// Hydrogen = false;
		x = y = z = 0.0;
		// length = angle = 0.0;
		// markedDuplicate = false;
	}
};

//struct HBondAtom
struct HydrogenBond
{
	double length;
	double angle;
	struct thbAtom *hydrogen;
	struct thbAtom *donor;
	struct thbAtom *acceptor;
	// char dendrimer[80];
	struct HydrogenBond *Next;
	struct HydrogenBond *Previous;
	// bool Hydrogen;
	// char Name[80];
	bool markedDuplicate;
//

	// Assign default values
	public:
	HydrogenBond()
	{
		Next = Previous = NULL;
		hydrogen = NULL;
		donor = NULL;
		acceptor = NULL;
		length = angle = 0.0;
		markedDuplicate = false;
	}
};

class ListOfHBonds
{
	private:
		int size;
		struct HydrogenBond *Start;
		double Round (double r, double f);

	public:
		// List Builder
		ListOfHBonds();
		unsigned int AtomCount();
		// List Operations
		int AddAtEnd(struct HydrogenBond *Item);
		int AddAtStart(struct HydrogenBond *Item);
		struct HydrogenBond *Retrive(int pos);
		struct HydrogenBond *Last();
		struct HydrogenBond *First();
		bool DeleteList();
		bool Find( struct HydrogenBond *Item);
		unsigned int SwitchingCount();
		unsigned int ForcefieldCount();
		unsigned int MoleculeCount();
		unsigned int CountUniqStr( std::vector< std::string >s );
		bool IsSameAsFirst( struct HydrogenBond *Item);
		bool IsSameAsLast( struct HydrogenBond *Item);
		bool ClosedLoop(void);
		std::vector<double> MinimumImage( struct thbAtom *A,
		                                  std::vector<double> r,
		                                  struct PBC Cell);
		void PrintAll(void);
		void PrintAllPovRay(void);
		double PrintAll(std::ostream *out, struct PBC Cell, bool POVRAY);
};
#endif
