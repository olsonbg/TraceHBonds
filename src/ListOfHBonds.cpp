#include "ListOfHBonds.h"

ListOfHBonds::ListOfHBonds()
  : size(0), Beginning(NULL)
{
}

unsigned int ListOfHBonds::AtomCount()
{
	return size;
}


struct HydrogenBond *ListOfHBonds::Begin()
{
	return(Beginning);
}

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
	struct HydrogenBond *Item;

	Item           = NewItem;
	Item->Next     = NULL;
	Item->Previous = End();

	if ( Item->Previous != NULL )
		Item->Previous->Next = Item;

	if( Begin() == NULL )
		size += 3;
	else
		size += 2;

	// This is actually the first element.
	if( Item->Previous == NULL )
		Beginning = Item;


	return (size);
}

void ListOfHBonds::PrintAll(unsigned int TrjIdx)
{
	struct HydrogenBond *current;

	current = Begin();

	while (current != NULL )
	{
		OFmt col(9,4);
		std::cout << col << current->donor->x.at(TrjIdx) << " ";
		std::cout << col << current->donor->y.at(TrjIdx) << " ";
		std::cout << col << current->donor->z.at(TrjIdx);
		std::cout << " [" << current->donor->Type << "]";
		std::cout << "  " << current->donor->Molecule;
		std::cout << "  " << current->donor->ForceField << "\n";

		std::cout << col << current->hydrogen->x.at(TrjIdx) << " ";
		std::cout << col << current->hydrogen->y.at(TrjIdx) << " ";
		std::cout << col << current->hydrogen->z.at(TrjIdx);
		std::cout << " [" << current->hydrogen->Type << "]";
		std::cout << "  " << current->hydrogen->Molecule;
		std::cout << "  " << current->hydrogen->ForceField << "\n";

		if ( current == End() )
		{
			std::cout << col << current->acceptor->x.at(TrjIdx) << " ";
			std::cout << col << current->acceptor->y.at(TrjIdx) << " ";
			std::cout << col << current->acceptor->z.at(TrjIdx);
			std::cout << " [" << current->acceptor->Type << "]";
			std::cout << "  " << current->acceptor->Molecule;
			std::cout << "  " << current->acceptor->ForceField << "\n";
		}
		current = current->Next;
	}
	return;
}

void ListOfHBonds::PrintAllPovRay(unsigned int TrjIdx)
{
	struct HydrogenBond *current;

	current = Begin();

	std::cout << "sphere_sweep {\n\tlinear_spline\n\t11\n" << "\n";
	while (current != NULL )
	{
		OFmt col(9,4);
		std::cout << "\t<";
		std::cout << col << current->donor->x.at(TrjIdx) << ", ";
		std::cout << col << current->donor->y.at(TrjIdx) << ",  ";
		std::cout << col << current->donor->z.at(TrjIdx) << ">,0.7" << "\n";

		std::cout << "\t<";
		std::cout << col << current->hydrogen->x.at(TrjIdx) << ", ";
		std::cout << col << current->hydrogen->y.at(TrjIdx) << ",  ";
		std::cout << col << current->hydrogen->z.at(TrjIdx) << ">,0.7" << "\n";

		if ( current == End() )
		{
			std::cout << "\t<";
			std::cout << col << current->acceptor->x.at(TrjIdx) << ", ";
			std::cout << col << current->acceptor->y.at(TrjIdx) << ",  ";
			std::cout << col << current->acceptor->z.at(TrjIdx) << ">,0.7" << "\n";
		}

		current = current->Next;
	}
	std::cout << "\ttolerance 0.07\n\ttexture{ChainLength11}" << "\n";
	return;
}

double ListOfHBonds::Round (double r,double f=1.0)
{
	return (r > 0.0) ? floor(r*f + 0.5)/f : ceil(r*f - 0.5)/f;
}

// Minimum Image distance of atom A from coordinate r.
// returns x,y,z as a vector of doubles.
std::vector<double> ListOfHBonds::MinimumImage( struct thbAtom *A,
                                                unsigned int TrjIdx,
                                                std::vector<double> r,
                                                struct PBC Cell)
{
	std::vector<double> d;
	double Nx, Ny, Nz;
	
		// Take care of periodic boundary conditions.
		// Minimum Image calculation.
		Nx = Round( (A->x.at(TrjIdx) - r[0])/Cell.x.at(TrjIdx));
		Ny = Round( (A->y.at(TrjIdx) - r[1])/Cell.y.at(TrjIdx));
		Nz = Round( (A->z.at(TrjIdx) - r[2])/Cell.z.at(TrjIdx));

		d.push_back( A->x.at(TrjIdx) - Nx*Cell.x.at(TrjIdx) );
		d.push_back( A->y.at(TrjIdx) - Ny*Cell.y.at(TrjIdx) );
		d.push_back( A->z.at(TrjIdx) - Nz*Cell.z.at(TrjIdx) );

		return(d);
}

double ListOfHBonds::PrintAll( std::ostream *out,
                               struct PBC Cell,
                               unsigned int TrjIdx,
                               bool POVRAY )
{
	struct HydrogenBond *current;
	std::vector<double>r;
	double initial_x, initial_y, initial_z;
	double EndToEndLength;

	current = Begin();

	initial_x = current->donor->x.at(TrjIdx);
	initial_y = current->donor->y.at(TrjIdx);
	initial_z = current->donor->z.at(TrjIdx);

	r.push_back(current->donor->x.at(TrjIdx));
	r.push_back(current->donor->y.at(TrjIdx));
	r.push_back(current->donor->z.at(TrjIdx));

	if (POVRAY) 
		*out << "sphere_sweep {\n\tlinear_spline\n\t" << AtomCount() 
		          << "," << "\n";

	OFmt colX(9,4);
	OFmt colY(9,4);
	OFmt colZ(9,4);
	while (current != NULL )
	{
		if (POVRAY)
		{
			r = MinimumImage( current->donor, TrjIdx, r, Cell );
			*out << "\t<";
			*out << colX << Round(r[0],10000.0) << ", ";
			*out << colY << Round(r[1],10000.0) << ", ";
			*out << colZ << Round(r[2],10000.0) << ">,0.7" << "\n";

			r = MinimumImage( current->hydrogen, TrjIdx, r, Cell );
			*out << "\t<";
			*out << colX << Round(r[0],10000.0) << ", ";
			*out << colY << Round(r[1],10000.0) << ", ";
			*out << colZ << Round(r[2],10000.0) << ">,0.7" << "\n";

			if ( current == End() )
			{
				r = MinimumImage( current->acceptor, TrjIdx, r, Cell );
				*out << "\t<";
				*out << colX << Round(r[0],10000.0) << ", ";
				*out << colY << Round(r[1],10000.0) << ", ";
				*out << colZ << Round(r[2],10000.0) << ">,0.7" << "\n";
			}
		}
		else
		{
			r = MinimumImage( current->donor, TrjIdx, r, Cell );
			*out << colX << Round(r[0],10000.0) << " ";
			*out << colY << Round(r[1],10000.0) << " ";
			*out << colZ << Round(r[2],10000.0);

			*out << " [" << current->donor->Type << "]";
			*out << "  " << current->donor->Molecule;
			*out << "  " << current->donor->Name;
			*out << "  " << current->donor->ForceField << "\n";

			r = MinimumImage( current->hydrogen, TrjIdx, r, Cell );
			*out << colX << Round(r[0],10000.0) << " ";
			*out << colY << Round(r[1],10000.0) << " ";
			*out << colZ << Round(r[2],10000.0);

			*out << " [" << current->hydrogen->Type << "]";
			*out << "  " << current->hydrogen->Molecule;
			*out << "  " << current->hydrogen->Name;
			*out << "  " << current->hydrogen->ForceField << "\n";

			if ( current == End() )
			{
				r = MinimumImage( current->acceptor, TrjIdx, r, Cell );
				*out << colX << Round(r[0],10000.0) << " ";
				*out << colY << Round(r[1],10000.0) << " ";
				*out << colZ << Round(r[2],10000.0);

				*out << " [" << current->acceptor->Type << "]";
				*out << "  " << current->acceptor->Molecule;
				*out << "  " << current->acceptor->Name;
				*out << "  " << current->acceptor->ForceField << "\n";
			}
		}
		current = current->Next;
		// counter++;
	}

	if (POVRAY)
		*out << "\ttolerance 0.07\n\ttexture{ChainLength" << AtomCount() 
		          << "}\n}" << "\n";

	EndToEndLength = sqrt( pow(initial_x-r[0],2) +
	                       pow(initial_y-r[1],2) +
	                       pow(initial_z-r[2],2) );

	return(EndToEndLength);
	// return(counter);
}

int ListOfHBonds::AddAtBegin(struct HydrogenBond *NewItem)
{
	struct HydrogenBond *Item;

	Item           = NewItem;
	Item->Next     = Begin();
	Item->Previous = NULL;

	if ( Item->Next != NULL )
		(Item->Next)->Previous = Item;


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
	// if ( (Item == NULL) || (size == 0) )
	//     return(false);

	if ( Begin()->donor == Item->acceptor )
		return(true);

	return(false);
}

// If the donor of Item is the same atom the the acceptor of the last
// element in the list, then Item should come after End().
bool ListOfHBonds::linksAtEnd( struct HydrogenBond *Item )
{
	// if ( (Item == NULL) || (size == 0) )
	//     return(false);

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

	return(true);
}
