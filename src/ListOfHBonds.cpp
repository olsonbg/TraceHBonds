#include "ListOfHBonds.h"

ListOfHBonds::ListOfHBonds()
  : size(0), HBsize(0), Beginning(NULL)
{
}

unsigned int ListOfHBonds::AtomCount() {
	return size; }

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

std::vector< double >
ListOfHBonds::donorAcceptorDistances()
{
	std::vector< double >distances;
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

	std::vector< Point > *d;
	std::vector< Point > *a = &(current->acceptor->p);

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
	struct HydrogenBond *Item;

	Item           = NewItem;
	Item->Next     = NULL;
	Item->Previous = End();

	if ( Item->Previous != NULL )
		Item->Previous->Next = Item;

	++HBsize;

	if( Begin() == NULL )
		size += 3;
	else
		size += 2;

	// This is actually the first element.
	if( Item->Previous == NULL )
		Beginning = Item;


	return (size);
}

#ifdef DEBUG
void ListOfHBonds::PrintAll(unsigned int TrjIdx)
{
	struct HydrogenBond *current;

	current = Begin();

	while (current != NULL )
	{
		OFmt col(9,4);
		std::cout << col << current->donor->p.at(TrjIdx).x() << " ";
		std::cout << col << current->donor->p.at(TrjIdx).y() << " ";
		std::cout << col << current->donor->p.at(TrjIdx).z();
		std::cout << " [" << current->donor->Type << "]";
		std::cout << "  " << current->donor->Molecule;
		std::cout << "  " << current->donor->ForceField << "\n";

		std::cout << col << current->hydrogen->p.at(TrjIdx).x() << " ";
		std::cout << col << current->hydrogen->p.at(TrjIdx).y() << " ";
		std::cout << col << current->hydrogen->p.at(TrjIdx).z();
		std::cout << " [" << current->hydrogen->Type << "]";
		std::cout << "  " << current->hydrogen->Molecule;
		std::cout << "  " << current->hydrogen->ForceField << "\n";

		if ( current == End() )
		{
			std::cout << col << current->acceptor->p.at(TrjIdx).x() << " ";
			std::cout << col << current->acceptor->p.at(TrjIdx).y() << " ";
			std::cout << col << current->acceptor->p.at(TrjIdx).z();
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
		std::cout << col << current->donor->p.at(TrjIdx).x() << ", ";
		std::cout << col << current->donor->p.at(TrjIdx).y() << ",  ";
		std::cout << col << current->donor->p.at(TrjIdx).z() << ">,0.7" << "\n";

		std::cout << "\t<";
		std::cout << col << current->hydrogen->p.at(TrjIdx).x() << ", ";
		std::cout << col << current->hydrogen->p.at(TrjIdx).y() << ",  ";
		std::cout << col << current->hydrogen->p.at(TrjIdx).z() << ">,0.7" << "\n";

		if ( current == End() )
		{
			std::cout << "\t<";
			std::cout << col << current->acceptor->p.at(TrjIdx).x() << ", ";
			std::cout << col << current->acceptor->p.at(TrjIdx).y() << ",  ";
			std::cout << col << current->acceptor->p.at(TrjIdx).z() << ">,0.7" << "\n";
		}

		current = current->Next;
	}
	std::cout << "\ttolerance 0.07\n\ttexture{ChainLength11}" << "\n";
	return;
}
#endif

double ListOfHBonds::Round (double r,double f=1.0)
{
	return (r > 0.0) ? floor(r*f + 0.5)/f : ceil(r*f - 0.5)/f;
}

// Minimum Image vector of atom A from coordinate r.
Point ListOfHBonds::MinimumImage( struct thbAtom *A,
                                  unsigned int TrjIdx,
                                  Point r,
                                  struct PBC Cell)
{
	// std::vector<double> d;
	// double Nx, Ny, Nz;

	// Take care of periodic boundary conditions.
	// Minimum Image calculation.
	// Nx = Round( (A->x.at(TrjIdx) - r[0])/Cell.x.at(TrjIdx));
	// Ny = Round( (A->y.at(TrjIdx) - r[1])/Cell.y.at(TrjIdx));
	// Nz = Round( (A->z.at(TrjIdx) - r[2])/Cell.z.at(TrjIdx));

	// d.push_back( A->x.at(TrjIdx) - Nx*Cell.x.at(TrjIdx) );
	// d.push_back( A->y.at(TrjIdx) - Ny*Cell.y.at(TrjIdx) );
	// d.push_back( A->z.at(TrjIdx) - Nz*Cell.z.at(TrjIdx) );

	// return A->p.at(TrjIdx).minimumImage(r, Cell.p.at(TrjIdx) );
	return r.minimumImage( A->p.at(TrjIdx), Cell.p.at(TrjIdx) ) ;
}

double ListOfHBonds::PrintAll( std::ostream *out,
                               struct PBC Cell,
                               unsigned int TrjIdx,
                               bool POVRAY )
{
	struct HydrogenBond *current;
	double EndToEndLength;

	current = Begin();

	Point initial( current->donor->p.at(TrjIdx).x(),
	               current->donor->p.at(TrjIdx).y(),
	               current->donor->p.at(TrjIdx).z() );

	
	Point r( current->donor->p.at(TrjIdx).x(),
	         current->donor->p.at(TrjIdx).y(),
	         current->donor->p.at(TrjIdx).z() );

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
			if ( current != Begin() )
				r = r + MinimumImage( current->donor, TrjIdx, r, Cell );

			*out << "\t<";
			*out << colX << Round(r.x(),10000.0) << ", ";
			*out << colY << Round(r.y(),10000.0) << ", ";
			*out << colZ << Round(r.z(),10000.0) << ">,0.7" << "\n";

			r = r + MinimumImage( current->hydrogen, TrjIdx, r, Cell );
			*out << "\t<";
			*out << colX << Round(r.x(),10000.0) << ", ";
			*out << colY << Round(r.y(),10000.0) << ", ";
			*out << colZ << Round(r.z(),10000.0) << ">,0.7" << "\n";

			if ( current == End() )
			{
				r = r + MinimumImage( current->acceptor, TrjIdx, r, Cell );
				*out << "\t<";
				*out << colX << Round(r.x(),10000.0) << ", ";
				*out << colY << Round(r.y(),10000.0) << ", ";
				*out << colZ << Round(r.z(),10000.0) << ">,0.7" << "\n";
			}
		}
		else
		{
			if ( current != Begin() )
				r = r + MinimumImage( current->donor, TrjIdx, r, Cell );

			*out << colX << Round(r.x(),10000.0) << " ";
			*out << colY << Round(r.y(),10000.0) << " ";
			*out << colZ << Round(r.z(),10000.0);

			*out << " [" << current->donor->Type << "]";
			*out << "  " << current->donor->Molecule;
			*out << "  " << current->donor->Name;
			*out << "  " << current->donor->ForceField << "\n";

			r = r + MinimumImage( current->hydrogen, TrjIdx, r, Cell );
			*out << colX << Round(r.x(),10000.0) << " ";
			*out << colY << Round(r.y(),10000.0) << " ";
			*out << colZ << Round(r.z(),10000.0);

			*out << " [" << current->hydrogen->Type << "]";
			*out << "  " << current->hydrogen->Molecule;
			*out << "  " << current->hydrogen->Name;
			*out << "  " << current->hydrogen->ForceField << "\n";

			if ( current == End() )
			{
				r = r + MinimumImage( current->acceptor, TrjIdx, r, Cell );
				*out << colX << Round(r.x(),10000.0) << " ";
				*out << colY << Round(r.y(),10000.0) << " ";
				*out << colZ << Round(r.z(),10000.0);

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

	// EndToEndLength = sqrt( pow(initial_x-r[0],2) +
	//                        pow(initial_y-r[1],2) +
	//                        pow(initial_z-r[2],2) );

	EndToEndLength = initial.distance(r);
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
	HBsize = 0;

	return(true);
}
