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

void ListOfHBonds::PrintAll()
{
	struct HydrogenBond *current;

	current = Begin();

	while (current != NULL )
	{
		OFmt col(9,4);
		std::cout << col << current->donor->x << " ";
		std::cout << col << current->donor->y << " ";
		std::cout << col << current->donor->z;
		std::cout << " [" << current->donor->Type << "]";
		std::cout << "  " << current->donor->Molecule;
		std::cout << "  " << current->donor->ForceField << std::endl;

		std::cout << col << current->hydrogen->x << " ";
		std::cout << col << current->hydrogen->y << " ";
		std::cout << col << current->hydrogen->z;
		std::cout << " [" << current->hydrogen->Type << "]";
		std::cout << "  " << current->hydrogen->Molecule;
		std::cout << "  " << current->hydrogen->ForceField << std::endl;

		if ( current == End() )
		{
			std::cout << col << current->acceptor->x << " ";
			std::cout << col << current->acceptor->y << " ";
			std::cout << col << current->acceptor->z;
			std::cout << " [" << current->acceptor->Type << "]";
			std::cout << "  " << current->acceptor->Molecule;
			std::cout << "  " << current->acceptor->ForceField << std::endl;
		}
		current = current->Next;
	}
	return;
}

void ListOfHBonds::PrintAllPovRay()
{
	struct HydrogenBond *current;

	current = Begin();

	std::cout << "sphere_sweep {\n\tlinear_spline\n\t11\n" << std::endl;
	while (current != NULL )
	{
		OFmt col(9,4);
		std::cout << "\t<";
		std::cout << col << current->donor->x << ", ";
		std::cout << col << current->donor->y << ",  ";
		std::cout << col << current->donor->z << ">,0.7" << std::endl;

		std::cout << "\t<";
		std::cout << col << current->hydrogen->x << ", ";
		std::cout << col << current->hydrogen->y << ",  ";
		std::cout << col << current->hydrogen->z << ">,0.7" << std::endl;

		if ( current == End() )
		{
			std::cout << "\t<";
			std::cout << col << current->acceptor->x << ", ";
			std::cout << col << current->acceptor->y << ",  ";
			std::cout << col << current->acceptor->z << ">,0.7" << std::endl;
		}

		current = current->Next;
	}
	std::cout << "\ttolerance 0.07\n\ttexture{ChainLength11}" << std::endl;
	return;
}

double ListOfHBonds::Round (double r,double f=1.0)
{
	return (r > 0.0) ? floor(r*f + 0.5)/f : ceil(r*f - 0.5)/f;
}

// Minimum Image distance of atom A from coordinate r.
// returns x,y,z as a vector of doubles.
std::vector<double> ListOfHBonds::MinimumImage( struct thbAtom *A,
                                                std::vector<double> r,
                                                struct PBC Cell)
{
	std::vector<double> d;
	double Nx, Ny, Nz;
	
		// Take care of periodic boundary conditions.
		// Minimum Image calculation.
		Nx = Round( (A->x - r[0])/Cell.x);
		Ny = Round( (A->y - r[1])/Cell.y);
		Nz = Round( (A->z - r[2])/Cell.z);

		d.push_back( A->x - Nx*Cell.x );
		d.push_back( A->y - Ny*Cell.y );
		d.push_back( A->z - Nz*Cell.z );

		return(d);
}

double ListOfHBonds::PrintAll( std::ostream *out, struct PBC Cell, bool POVRAY )
{
	struct HydrogenBond *current;
	std::vector<double>r;
	double initial_x, initial_y, initial_z;
	double EndToEndLength;

	current = Begin();

	initial_x = current->donor->x;
	initial_y = current->donor->y;
	initial_z = current->donor->z;

	r.push_back(current->donor->x);
	r.push_back(current->donor->y);
	r.push_back(current->donor->z);

	if (POVRAY) 
		*out << "sphere_sweep {\n\tlinear_spline\n\t" << AtomCount() 
		          << "," << std::endl;

	OFmt colX(9,4);
	OFmt colY(9,4);
	OFmt colZ(9,4);
	while (current != NULL )
	{
		if (POVRAY)
		{
			r = MinimumImage( current->donor, r, Cell );
			*out << "\t<";
			*out << colX << Round(r[0],10000.0) << ", ";
			*out << colY << Round(r[1],10000.0) << ", ";
			*out << colZ << Round(r[2],10000.0) << ">,0.7" << std::endl;

			r = MinimumImage( current->hydrogen, r, Cell );
			*out << "\t<";
			*out << colX << Round(r[0],10000.0) << ", ";
			*out << colY << Round(r[1],10000.0) << ", ";
			*out << colZ << Round(r[2],10000.0) << ">,0.7" << std::endl;

			if ( current == End() )
			{
				r = MinimumImage( current->acceptor, r, Cell );
				*out << "\t<";
				*out << colX << Round(r[0],10000.0) << ", ";
				*out << colY << Round(r[1],10000.0) << ", ";
				*out << colZ << Round(r[2],10000.0) << ">,0.7" << std::endl;
			}
		}
		else
		{
			r = MinimumImage( current->donor, r, Cell );
			*out << colX << Round(r[0],10000.0) << " ";
			*out << colY << Round(r[1],10000.0) << " ";
			*out << colZ << Round(r[2],10000.0);

			*out << " [" << current->donor->Type << "]";
			*out << "  " << current->donor->Molecule;
			*out << "  " << current->donor->Name;
			*out << "  " << current->donor->ForceField << std::endl;

			r = MinimumImage( current->hydrogen, r, Cell );
			*out << colX << Round(r[0],10000.0) << " ";
			*out << colY << Round(r[1],10000.0) << " ";
			*out << colZ << Round(r[2],10000.0);

			*out << " [" << current->hydrogen->Type << "]";
			*out << "  " << current->hydrogen->Molecule;
			*out << "  " << current->hydrogen->Name;
			*out << "  " << current->hydrogen->ForceField << std::endl;

			if ( current == End() )
			{
				r = MinimumImage( current->acceptor, r, Cell );
				*out << colX << Round(r[0],10000.0) << " ";
				*out << colY << Round(r[1],10000.0) << " ";
				*out << colZ << Round(r[2],10000.0);

				*out << " [" << current->acceptor->Type << "]";
				*out << "  " << current->acceptor->Molecule;
				*out << "  " << current->acceptor->Name;
				*out << "  " << current->acceptor->ForceField << std::endl;
			}
		}
		current = current->Next;
		// counter++;
	}

	if (POVRAY)
		*out << "\ttolerance 0.07\n\ttexture{ChainLength" << AtomCount() 
		          << "}\n}" << std::endl;

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

// TODO: Rename this.
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

// TODO: Rename this.
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
