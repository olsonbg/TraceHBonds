#include "ReadData.h"

double Vector_Length( struct HBondAtom *a, struct HBondAtom *b );

int ReadData( const char *filename,
              std::vector<struct HBondAtom *> *Atoms1,
              std::vector<struct HBondAtom *> *Atoms2,
              std::vector<struct HBondAtom *> *Atoms3,
              struct PBC *Cell)
{
	char line[256];
	struct HBondAtom *Atom1, *Atom2, *Atom3;

	// Donor --- Hydrogen ... Acceptor
	// ... Denotes the Hydrogen bond.
//	Atom1 = new struct HBondAtom; // Hydrogen
//	Atom2 = new struct HBondAtom; // Acceptor Oxygen
//	Atom3 = new struct HBondAtom; // Donor Oxygen
//	Atom1->Length = Hydrogen - Donor Length.
//	Atom2->Length = Hydrogen - Acceptor Length (HBond length).

	// store the 1st 6 bytes of file to determine the filetype (magic number).
	uint8_t magic[6];
	bool STDIN=false;
	FILE *fp ;

	if ( !strncmp(filename,"-",1) )
		STDIN=true;

	if ( !STDIN )
	{
		fp = fopen(filename,"r");

		if( fp == NULL )
		{
			std::perror(filename);
			return(1) ;
		}
		if ( !STDIN && (fread(magic, 1, 6, fp) != 6) )
		{
			std::cerr << "Error reading magic number." << std::endl;
			fclose(fp);
			return(1);
		}
		fclose(fp);
	}

	boost::iostreams::filtering_stream<boost::iostreams::input> in;
	std::ifstream ifp;
	if (ifp == NULL)
		return(1);

	// Determine filetype, then use appropriate boost filter.
#ifdef USE_ZLIB
	if ( !STDIN &&
	     (magic[0] == 0x1f) && (magic[1] == 0x8b) && (magic[2] == 0x08) )
	{
		// This is a gzip file.
		in.push(boost::iostreams::gzip_decompressor());
		ifp.open(filename,std::ios::in|std::ios::binary);
	}
#endif
#ifdef USE_BZIP2
	if ( !STDIN &&
	     (magic[0] == 0x42) && (magic[1] == 0x5a) && (magic[2] == 0x68) )
	{
		// This is a bzip2 file.
		in.push(boost::iostreams::bzip2_decompressor());
		ifp.open(filename,std::ios::in|std::ios::binary);
	}
#endif
	if ( !ifp.is_open() && !STDIN ) // Assume a plain text file if 
		                            // ifp is not open.
		ifp.open(filename,std::ios::in);

	if (STDIN)
		in.push(std::cin);
	else
		in.push(ifp);

	in.getline(line,255);

	int n;
	if ( !in.eof() )
	{
		n = sscanf( line,
		            "%le %le %le %le %le %le",
		            &Cell->x, &Cell->y, &Cell->z,
		            &Cell->alpha, &Cell->beta, &Cell->gamma);
		if ( n != 6 )
		{
			std::cerr << "Error on line 1. Expected 6 arguments, read " << n;
			std::cerr << "." <<std::endl;
			std::cerr << " L1: '" << line << "'" << std::endl;
			ifp.close();
			return(1);
		}
	}

	in.getline(line,255);

	int linen=1;
	while ( !in.eof() )
	{
		// Donor --- Hydrogen ... Acceptor
		// ... Denotes the Hydrogen bond.
		Atom1 = new struct HBondAtom; // Hydrogen
		Atom2 = new struct HBondAtom; // Acceptor Oxygen
		Atom3 = new struct HBondAtom; // Donor Oxygen
		n = sscanf( line,
		            "%78s\t%78s\t%78s\t%78s\t{%le %le %le}\t{%le %le %le}\t{%le %le %le}\t%le\t%le\t%78s\t%78s",
		            Atom1->Name, Atom2->Name, Atom3->Name,
		            Atom3->ffType,
		            &Atom1->x, &Atom1->y, &Atom1->z,
		            &Atom2->x, &Atom2->y, &Atom2->z,
		            &Atom3->x, &Atom3->y, &Atom3->z,
		            &Atom2->length, &Atom2->angle,
		            Atom1->dendrimer, Atom2->dendrimer);
		linen++;
		if ( n == 17 )
		{
			strncpy(Atom3->dendrimer,Atom1->dendrimer,80) ;
			Atom1->Hydrogen = true;
			Atom2->Hydrogen = Atom3->Hydrogen = false;
			Atom1->length = Vector_Length( Atom1, Atom3 );
			strcpy(Atom1->ffType,"h1o");
			strcpy(Atom2->ffType,"o2h");
			// Save atoms to the vector.
			Atoms1->push_back(Atom1);
			Atoms2->push_back(Atom2);
			Atoms3->push_back(Atom3);
		}
		else // Not an Hbond line.
		{
			// Check if this line is another PBC line.
			n = sscanf( line,
			            "%le %le %le %le %le %le",
			            &Cell->x, &Cell->y, &Cell->z,
			            &Cell->alpha, &Cell->beta, &Cell->gamma);
			if ( n != 6 ) // Unknown line, abort.
			{
			std::cerr << "Error on line " << linen;
			std::cerr << ". Expected either 6, or 17 arguments, read " <<n<< "."<<std::endl;
			std::cerr << " L" << linen << ": '" << line << "'" << std::endl;
			return(1);
			}

			delete Atom1;
			delete Atom2;
			delete Atom3;
		}
		in.getline(line,255);
	}
	ifp.close();

	return(0);
}

double Vector_Length( struct HBondAtom *a, struct HBondAtom *b )
{
	double length;

	length = pow(a->x - b->x,2) +
	         pow(a->y - b->y,2) +
	         pow(a->z - b->z,2);

	length = sqrt(length);

	return(length);
}
