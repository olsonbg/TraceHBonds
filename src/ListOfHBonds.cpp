#include "ListOfHBonds.h"

ListOfHBonds::ListOfHBonds()
  : size(0), Start(NULL)
{
}

unsigned int ListOfHBonds::Count()
{
	return size;
}


struct HBondAtom *ListOfHBonds::First()
{
	return(Start);
}

struct HBondAtom *ListOfHBonds::Last()
{
	struct HBondAtom *current = Start;

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
	char *theMolecule;
	struct HBondAtom *current = Start;

	if (current == NULL)
		return(0);

	theMolecule = current->dendrimer;
	while (current->Next != NULL)
	{
		current = current->Next;
		if ( strcmp(theMolecule,current->dendrimer) != 0 ) 
		{
			theMolecule = current->dendrimer;
			Switches++;
		}
	}

	return(Switches);
}

unsigned int ListOfHBonds::CountUniqStr( std::vector<char *>s )
{
	std::map<std::string, int>word_count;

	for(unsigned int i=0; i < s.size(); ++i)
		word_count[ s[i] ]++;

//	typedef std::map<std::string, int>::iterator iter;
//	iter end = word_count.end();
//	for(iter it=word_count.begin(); it != end; ++it)
//		std::cout << it->first << ", count =" << it->second << std::endl;

	return(word_count.size());
}

unsigned int ListOfHBonds::MoleculeCount()
{
	std::vector<char *>mols;
	struct HBondAtom *current = Start;

	while (current != NULL)
	{
		mols.push_back(current->dendrimer);
		current = current->Next;
	}

	return(CountUniqStr(mols));
}

unsigned int ListOfHBonds::ForcefieldCount()
{
	std::vector<char *>ff;
	struct HBondAtom *current = Start;

	while (current != NULL)
	{
		ff.push_back(current->ffType);
		current = current->Next;
	}

	return(CountUniqStr(ff));
}

int ListOfHBonds::AddAtEnd(struct HBondAtom *NewItem)
{
	struct HBondAtom *Item;

	Item           = NewItem;
	Item->Next     = NULL;
	Item->Previous = Last();

	if ( Item->Previous != NULL )
		Item->Previous->Next = Item;

	if( Item->Previous == NULL )
		Start = Item;

	return size++;
}

void ListOfHBonds::PrintAll()
{
	struct HBondAtom *current;

	current = First();

	while (current != NULL )
	{
		OFmt col(9,4);
		// printf("% 9.4f % 9.4f % 9.4f", current->x, current->y, current->z);
		std::cout << col << current->x << " ";
		std::cout << col << current->y << " ";
		std::cout << col << current->z;
		if ( current->Hydrogen )
			std::cout << " [H]";
		else
			std::cout << " [O]";

		std::cout << "  " << current->dendrimer;
		std::cout << "  " << current->ffType << std::endl;
		current = current->Next;
	}
	return;
}

void ListOfHBonds::PrintAllPovRay()
{
	struct HBondAtom *current;

	current = First();

	std::cout << "sphere_sweep {\n\tlinear_spline\n\t11\n" << std::endl;
	while (current != NULL )
	{
		// printf("\t<% 9.4f, % 9.4f, % 9.4f>,0.7", current->x,
		//                                          current->y,
		//                                          current->z);
		OFmt col(9,4);
		std::cout << "\t<";
		std::cout << col << current->x << ", ";
		std::cout << col << current->y << ",  ";
		std::cout << col << current->z << ">,0.7";
		current = current->Next;
	}
	std::cout << "\ttolerance 0.07\n\ttexture{ChainLength11}" << std::endl;
	return;
}

double ListOfHBonds::Round (double r)
{
	return (r > 0.0) ? floor(r + 0.5) : ceil(r - 0.5);
}

double ListOfHBonds::Round (double r,double f)
{
	return (r > 0.0) ? floor(r*f + 0.5)/f : ceil(r*f - 0.5)/f;
}

double ListOfHBonds::PrintAll( std::ostream *out, struct PBC Cell, bool POVRAY )
{
	struct HBondAtom *current;
	double tmp_x,tmp_y, tmp_z;
	double initial_x, initial_y, initial_z;
	double Nx, Ny, Nz;
	double EndToEndLength;

	// double counter=0;
	current = First();

	initial_x = tmp_x = current->x;
	initial_y = tmp_y = current->y;
	initial_z = tmp_z = current->z;

	if (POVRAY) 
		*out << "sphere_sweep {\n\tlinear_spline\n\t" << Count() 
		          << "," << std::endl;

	OFmt colX(9,4);
	OFmt colY(9,4);
	OFmt colZ(9,4);
	while (current != NULL )
	{
		Nx = Round( (current->x - tmp_x)/Cell.x);
		Ny = Round( (current->y - tmp_y)/Cell.y);
		Nz = Round( (current->z - tmp_z)/Cell.z);

		tmp_x = current->x - Nx*Cell.x;
		tmp_y = current->y - Ny*Cell.y;
		tmp_z = current->z - Nz*Cell.z;

		if (POVRAY)
		{
			// printf("\t<% 9.4lf, % 9.4lf, % 9.4lf>,0.7", Round(tmp_x,10000.0),
			//                                             Round(tmp_y,10000.0),
			//                                             Round(tmp_z,10000.0));
			*out << "\t<";
			*out << colX << Round(tmp_x,10000.0) << ", ";
			*out << colY << Round(tmp_y,10000.0) << ", ";
			*out << colZ << Round(tmp_z,10000.0) << ">,0.7";
		}
		else
		{
			*out << colX << Round(tmp_x,10000.0) << " ";
			*out << colY << Round(tmp_y,10000.0) << " ";
			*out << colZ << Round(tmp_z,10000.0);
			// printf("% 9.4lf % 9.4lf % 9.4lf", Round(tmp_x,10000.0),
			//                                   Round(tmp_y,10000.0),
			//                                   Round(tmp_z,10000.0));

			if ( current->Hydrogen )
				*out << " [H]";
			else
				*out << " [O]";

			*out << "  " << current->dendrimer;
			*out << "  " << current->Name;
			*out << "  " << current->ffType;
		}
		*out << std::endl;
		current = current->Next;
		// counter++;
	}

	if (POVRAY)
		*out << "\ttolerance 0.07\n\ttexture{ChainLength" << Count() 
		          << "}\n}" << std::endl;

	EndToEndLength = sqrt( pow(initial_x-tmp_x,2) +
	                       pow(initial_y-tmp_y,2) +
	                       pow(initial_z-tmp_z,2) );

	return(EndToEndLength);
	// return(counter);
}

int ListOfHBonds::AddAtStart(struct HBondAtom *NewItem)
{
	struct HBondAtom *Item;

	Item           = NewItem;
	Item->Next     = First();
	Item->Previous = NULL;

	if ( Item->Next != NULL )
		(Item->Next)->Previous = Item;

	Start = NewItem;

	return size++;
}

bool ListOfHBonds::Find( struct HBondAtom *Item )
{
	struct HBondAtom *current;

	if ( Item == NULL )
		return(false);

	for( current = First(); current != NULL; current = current->Next )
	{
		if ( (current->x == Item->x) &&
		     (current->y == Item->y) &&
		     (current->z == Item->z) )
			return(true);
	}

	return(false);
}

bool ListOfHBonds::IsSameAsFirst( struct HBondAtom *Item )
{
	struct HBondAtom *current;

	if ( (Item == NULL) || (size == 0) )
		return(false);

	current = First();

	if ( (current->x == Item->x) &&
	     (current->y == Item->y) &&
	     (current->z == Item->z) )
		return(true);

	return(false);
}

bool ListOfHBonds::IsSameAsLast( struct HBondAtom *Item )
{
	struct HBondAtom *current;

	if ( (Item == NULL) || (size == 0) )
		return(false);

	current = Last();

	if ( (current->x == Item->x) && 
	     (current->y == Item->y) && 
	     (current->z == Item->z) )
		return(true);

	return(false);
}

bool ListOfHBonds::ClosedLoop(void)
{
	return (IsSameAsLast(Start));
}

bool ListOfHBonds::DeleteList()
{
	Start = NULL;
	size = 0;

	return(true);
}
