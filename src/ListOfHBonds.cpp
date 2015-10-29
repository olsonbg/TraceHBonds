#include "ListOfHBonds.h"

ListOfHBonds::ListOfHBonds()
  : size(0), HBsize(0), Beginning(NULL)
{
}

unsigned int ListOfHBonds::AtomCount() {
	// return (2*HBsize + 1);
	return size;
}

unsigned int ListOfHBonds::HydrogenBondCount() {
	return HBsize; }

struct HydrogenBond *ListOfHBonds::Begin() {
	return(Beginning); }

struct HydrogenBond *ListOfHBonds::End()
{
	struct HydrogenBond *current = Begin();

	if ( current == NULL )
		return(NULL);

	while ( current->Next != NULL )
		current = current->Next;

	return(current);
}

// All the HydrogenBonds in this list should be of the same
// Trajectory Index.
unsigned int ListOfHBonds::TrajectoryIndex()
{
	return (Begin()->TrajIdx);
}

// How many times chain switched to another molecule
// Note: Can switch many times between 2 molecules.
//       Switching is not equal to the number of molecules
//       in the chain.
unsigned int ListOfHBonds::SwitchingCount()
{
	int Switches = 0;
	struct HydrogenBond *current = Begin();

	// Donor and Hydrogen are always on the same molecule,
	// So checking which molecule the Hydrogen and Acceptor
	// are on is sufficient.

	while (current != NULL)
	{
		if ( current->hydrogen->Molecule != current->acceptor->Molecule )
			Switches++;

		current = current->Next;
	}

	return(Switches);
}

unsigned int ListOfHBonds::CountUniqStr( std::vector< std::string >s )
{
	std::map<std::string, int>word_count;

	for(unsigned int i=0; i < s.size(); ++i)
		word_count[ s.at(i) ]++;

	return(word_count.size());
}

vd ListOfHBonds::donorAcceptorDistances()
{
	vd distances;
	struct HydrogenBond *current = Begin();

	while (current != NULL)
	{
		distances.push_back(current->acceptorDonorDistance);
		current = current->Next;
	}
	return(distances);
}

// Coordinates of the donor and acceptor atoms in the chain.
std::vector< Point *>
ListOfHBonds::nonHydrogenCoordinates()
{
	std::vector< Point *>p;
	struct HydrogenBond *current = Begin();

	while ( current != NULL)
	{
		p.push_back( &(current->donor->p.at(TrajectoryIndex())) );

		if ( current->Next == NULL )
		{
			p.push_back( &(current->acceptor->p.at(TrajectoryIndex())) );
		}

		current = current->Next;
	}


	return(p);
}
// Count how many unique molecules are on this chain of
// hydrogen bonds.
unsigned int ListOfHBonds::MoleculeCount()
{
	std::vector< std::string > molecules;
	struct HydrogenBond *current = Begin();

	while (current != NULL)
	{
		// Donor and Hydrogen are always on same molecule,
		// So using one of them is sufficient.
		molecules.push_back(current->hydrogen->Molecule);
		molecules.push_back(current->acceptor->Molecule);

		current = current->Next;
	}

	return(CountUniqStr(molecules));
}

unsigned int ListOfHBonds::ForcefieldCount()
{
	std::vector< std::string > forcefields;

	struct HydrogenBond *current = Begin();

	while (current != NULL)
	{
		forcefields.push_back(current->donor->ForceField);
		forcefields.push_back(current->hydrogen->ForceField);
		forcefields.push_back(current->acceptor->ForceField);

		current = current->Next;
	}

	return(CountUniqStr(forcefields));
}

int ListOfHBonds::AddAtEnd(struct HydrogenBond *NewItem)
{
	NewItem->Next     = NULL;
	NewItem->Previous = End();

	if ( NewItem->Previous != NULL )
		NewItem->Previous->Next = NewItem;

	++HBsize;

	if( Begin() == NULL )
		size += 3;
	else
		size += 2;

	// This is actually the first element.
	if( NewItem->Previous == NULL )
		Beginning = NewItem;


	return (size);
}

Point ListOfHBonds::Round (Point p, double f)
{
	Point r( Round(p.x(), f),
	         Round(p.y(), f),
	         Round(p.z(), f) );

	return(r);
}

double ListOfHBonds::Round (double r,double f)
{
	return (r > 0.0) ? floor(r*f + 0.5)/f : ceil(r*f - 0.5)/f;
}

// Minimum Image vector of atom A from coordinate r.
Point ListOfHBonds::MinimumImage( struct thbAtom *A,
                                  unsigned int TrjIdx,
                                  Point r,
                                  struct PBC Cell)
{
	return r.minimumImage( A->p.at(TrjIdx), Cell.p.at(TrjIdx) ) ;
}

double ListOfHBonds::PrintAll( std::ostream *out,
                               struct PBC Cell,
                               unsigned int TrjIdx,
                               unsigned char flags )
{
	struct HydrogenBond *current;
	double EndToEndLength;

	current = Begin();

	Point initial( current->donor->p.at(TrjIdx).x(),
	               current->donor->p.at(TrjIdx).y(),
	               current->donor->p.at(TrjIdx).z() );


	// TODO: add --incell flag

	// Use for non-wrapped
	// Point r( current->donor->p.at(TrjIdx).x(),
	//          current->donor->p.at(TrjIdx).y(),
	//          current->donor->p.at(TrjIdx).z() );

	// Use for wrapped (in-cell)
	// Put the first element of this chain inside the PBC cell
	Point r = MinimumImage(current->donor, TrjIdx, Point(0,0,0), Cell) +
	          Cell.p.at(TrjIdx)/2.0;
//	Point r(0.0,0.0,0.0);

	if ( flags & POVRAY )
		*out << "sphere_sweep {\n\tlinear_spline\n\t" << AtomCount()
		          << "," << "\n";
	else if ( flags & JSON )
		*out << "{\n    \"hbchain\": [\n";

	OFmt colX(9,4);
	OFmt colY(9,4);
	OFmt colZ(9,4);
	while (current != NULL )
	{
		Point A = Round(r, 10000.0);

		r = r + MinimumImage(current->hydrogen, TrjIdx, current->donor->p.at(TrjIdx), Cell);
		Point B = Round(r, 10000.0);

		r = r + MinimumImage(current->acceptor, TrjIdx, current->hydrogen->p.at(TrjIdx),Cell);
		Point C = Round(r, 10000.0);

		if ( flags & POVRAY )
		{
			*out << "\t<";
			*out << colX << A.x() << ", ";
			*out << colY << A.y() << ", ";
			*out << colZ << A.z() << ">,ChainRadius" << "\n";

			*out << "\t<";
			*out << colX << B.x() << ", ";
			*out << colY << B.y() << ", ";
			*out << colZ << B.z() << ">,ChainRadius" << "\n";

			if ( current == End() )
			{
				*out << "\t<";
				*out << colX << C.x() << ", ";
				*out << colY << C.y() << ", ";
				*out << colZ << C.z() << ">,ChainRadius" << "\n";
			}
		}
		else if ( flags & JSON )
		{
			*out << "        { \"location\": [ "
			     << colX << A.x() << ", "
			     << colY << A.y() << ", "
			     << colZ << A.z() << " ],\n"
			     << "          "
			     << "\"element\": \""     << current->donor->Type << "\", "
			     << " \"molecule\": \""   << current->donor->Molecule << "\",\n"
			     << "          "
			     << "\"name\": \""        << current->donor->Name << "\","
			     << " \"forcefield\": \"" << current->donor->ForceField
				 << "\" },\n";

			*out << "        { \"location\": [ "
			     << colX << B.x() << ", "
			     << colY << B.y() << ", "
			     << colZ << B.z() << " ],\n"
			     << "          "
			     << "\"element\": \""     << current->hydrogen->Type << "\", "
			     << " \"molecule\": \""   << current->hydrogen->Molecule << "\",\n"
			     << "          "
			     << "\"name\": \""        << current->hydrogen->Name << "\","
			     << " \"forcefield\": \"" << current->hydrogen->ForceField
				 << "\" },\n";

			if ( current == End() )
			{
				*out << "        { \"location\": [ "
				     << colX << C.x() << ", "
				     << colY << C.y() << ", "
				     << colZ << C.z() << " ],\n"
				     << "          "
				     << "\"element\": \""     << current->acceptor->Type << "\", "
				     << " \"molecule\": \""   << current->acceptor->Molecule << "\",\n"
				     << "          "
				     << "\"name\": \""        << current->acceptor->Name << "\","
				     << " \"forcefield\": \"" << current->acceptor->ForceField
				     << "\" }\n"
				     << "    ]\n},\n";
			}
		}
		else
		{
			*out << colX << A.x() << " ";
			*out << colY << A.y() << " ";
			*out << colZ << A.z();

			*out << " [" << current->donor->Type << "]";
			*out << "  " << current->donor->Molecule;
			*out << "  " << current->donor->Name;
			*out << "  " << current->donor->ForceField << "\n";

			*out << colX << B.x() << " ";
			*out << colY << B.y() << " ";
			*out << colZ << B.z();

			*out << " [" << current->hydrogen->Type << "]";
			*out << "  " << current->hydrogen->Molecule;
			*out << "  " << current->hydrogen->Name;
			*out << "  " << current->hydrogen->ForceField << "\n";

			if ( current == End() )
			{
				*out << colX << C.x() << " ";
				*out << colY << C.y() << " ";
				*out << colZ << C.z();

				*out << " [" << current->acceptor->Type << "]";
				*out << "  " << current->acceptor->Molecule;
				*out << "  " << current->acceptor->Name;
				*out << "  " << current->acceptor->ForceField << "\n";
			}
		}
		current = current->Next;
		// counter++;
	}

	if ( flags & POVRAY )
		*out << "\ttolerance 0.07\n\ttexture{ChainLength" << AtomCount()
		          << "}\n}" << "\n";


	EndToEndLength = initial.distance(r);
	return(EndToEndLength);
}

int ListOfHBonds::AddAtBegin(struct HydrogenBond *NewItem)
{
	NewItem->Next     = Beginning;
	NewItem->Previous = NULL;

	if ( NewItem->Next != NULL )
		(NewItem->Next)->Previous = NewItem;

	++HBsize;

	if ( Begin() == NULL )
		size += 3;
	else
		size += 2;

	Beginning = NewItem;

	return (size);
}

// If the donor, hydrogen, and acceptor pointer of current and Item
// point to the same, then we found it.
bool ListOfHBonds::Find( struct HydrogenBond *Item )
{
	struct HydrogenBond *current;

	if ( Item == NULL )
		return(false);

	for( current = Begin(); current != NULL; current = current->Next )
	{
		if ( current == Item )
			return(true);
	}

	return(false);
}

// If the acceptor of Item is the same atom as the donor of the first
// element in the list, then Item should come before Begin().
bool ListOfHBonds::linksAtBegin( struct HydrogenBond *Item )
{
	if ( Begin()->donor == Item->acceptor )
		return(true);

	return(false);
}

// If the donor of Item is the same atom the the acceptor of the last
// element in the list, then Item should come after End().
bool ListOfHBonds::linksAtEnd( struct HydrogenBond *Item )
{
	if ( End()->acceptor == Item->donor )
		return(true);

	return(false);
}

bool ListOfHBonds::ClosedLoop(void)
{
	return ( linksAtEnd(Begin()) );
}

bool ListOfHBonds::DeleteList()
{
	Beginning = NULL;
	size = 0;
	HBsize = 0;

	return(true);
}
