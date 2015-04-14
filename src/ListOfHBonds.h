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

struct HBondAtom
{
	double x,y,z;
	double length;
	double angle;
	char dendrimer[80];
	struct HBondAtom *Next;
	struct HBondAtom *Previous;
	bool Hydrogen;
	char ffType[80];
	char Name[80];
	bool markedDuplicate;

	// Assign default values
	public:
	HBondAtom()
	{
		Hydrogen = false;
		Next = Previous = NULL;
		length = angle = 0.0;
		x = 0.0;
		y = 0.0;
		z = 0.0;
		markedDuplicate = false;
	}
};

class ListOfHBonds
{
	private:
		int size;
		struct HBondAtom *Start;
		double Round (double r);
		double Round (double r, double f);

	public:
		// List Builder
		ListOfHBonds();
		unsigned int Count();
		// List Operations
		int AddAtEnd(struct HBondAtom *Item);
		int AddAtStart(struct HBondAtom *Item);
		struct HBondAtom *Retrive(int pos);
		struct HBondAtom *Last();
		struct HBondAtom *First();
		bool DeleteList();
		bool Find( struct HBondAtom *Item);
		unsigned int SwitchingCount();
		unsigned int ForcefieldCount();
		unsigned int MoleculeCount();
		unsigned int CountUniqStr( std::vector<char *>s );
		bool IsSameAsFirst( struct HBondAtom *Item);
		bool IsSameAsLast( struct HBondAtom *Item);
		bool ClosedLoop(void);
		void PrintAll(void);
		void PrintAllPovRay(void);
		double PrintAll(std::ostream *out, struct PBC Cell, bool POVRAY);
};
#endif
