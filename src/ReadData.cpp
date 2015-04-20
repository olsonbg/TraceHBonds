#include "ReadData.h"

struct thbAtom *AlreadyRecordedAtom( std::vector<struct thbAtom *> atom,
                                     struct thbAtom *find);

int ReadData( const char *filename,
              std::vector<struct HydrogenBond *> *hb,
              std::vector<struct thbAtom *> *atom,
              struct PBC *Cell )
{
	char line[256];
	struct HydrogenBond *HBond;
	struct thbAtom *hydrogen, *acceptor, *donor;
	unsigned int TrjIdx = 0;

	// Donor --- Hydrogen ... Acceptor
	// ... Denotes the Hydrogen bond.

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
		double CellX, CellY, CellZ;
		double CellAlpha, CellBeta, CellGamma;

		n = sscanf( line,
		            "%le %le %le %le %le %le",
		            &CellX, &CellY, &CellZ,
		            &CellAlpha, &CellBeta, &CellGamma);
		if ( n != 6 )
		{
			std::cerr << "Error on line 1. Expected 6 arguments, read " << n;
			std::cerr << "." <<std::endl;
			std::cerr << " L1: '" << line << "'" << std::endl;
			ifp.close();
			return(1);
		}
		if ( alloc_vector(Cell, 0.0, TrjIdx+1) )
		{
			Cell->x.at(TrjIdx) = CellX;
			Cell->y.at(TrjIdx) = CellY;
			Cell->z.at(TrjIdx) = CellZ;

			Cell->alpha.at(TrjIdx) = CellAlpha;
			Cell->beta.at(TrjIdx)  = CellBeta;
			Cell->gamma.at(TrjIdx) = CellGamma;
		}
	}

	in.getline(line,255);

	int linen=1;
	char acceptorForceField[80];
	char hydrogenName[80];
	char acceptorName[80];
	char donorName[80];
	char acceptorMolecule[80];
	char hydrogenMolecule[80];
	double donorX, donorY, donorZ;
	double hydrogenX, hydrogenY, hydrogenZ;
	double acceptorX, acceptorY, acceptorZ;
	while ( !in.eof() )
	{
		// Donor --- Hydrogen ... Acceptor
		// ... Denotes the Hydrogen bond.
		hydrogen = new struct thbAtom; // Atom1 of datafile.
		acceptor = new struct thbAtom; // Atom2 of datafile.
		donor    = new struct thbAtom; // Atom3 of datafile.
		HBond    = new struct HydrogenBond;

		n = sscanf( line,
		            "%78s\t%78s\t%78s\t%78s\t{%le %le %le}\t{%le %le %le}\t{%le %le %le}\t%le\t%le\t%78s\t%78s",
		            hydrogenName, acceptorName, donorName,
		            acceptorForceField,
		            &hydrogenX, &hydrogenY, &hydrogenZ,
		            &acceptorX, &acceptorY, &acceptorZ,
		            &donorX   , &donorY   , &donorZ,
		            &HBond->length, &HBond->angle,
		            hydrogenMolecule, acceptorMolecule);
		linen++;
		if ( n == 17 )
		{
			donor->Name        = donorName;
			hydrogen->Name     = hydrogenName;
			acceptor->Name     = acceptorName;

			donor->Type        = "O";
			hydrogen->Type     = "H";
			acceptor->Type     = "O";
			
			donor->Molecule    = hydrogenMolecule;
			hydrogen->Molecule = hydrogenMolecule;
			acceptor->Molecule = acceptorMolecule;
			
			acceptor->ForceField  = acceptorForceField;
			// The next two are an assumption, since we can not determine the
			// correct values from the datafile.
			hydrogen->ForceField = "h1o";
			donor->ForceField = "o2h";

			// Save atoms to the vector, unless the same atom has already been
			// saved.
			struct thbAtom *duplicate;
			
			duplicate = AlreadyRecordedAtom(*atom, donor);
			if ( duplicate == NULL )
			{
				if ( alloc_vector(donor, 0.0, TrjIdx+1) )
				{
					donor->x.at(TrjIdx)    = donorX;
					donor->y.at(TrjIdx)    = donorY;
					donor->z.at(TrjIdx)    = donorZ;
				}

				atom->push_back(donor);
				HBond->donor = donor;
			}
			else
			{
				if ( alloc_vector(duplicate, 0.0, TrjIdx+1) )
				{
					duplicate->x.at(TrjIdx)    = donorX;
					duplicate->y.at(TrjIdx)    = donorY;
					duplicate->z.at(TrjIdx)    = donorZ;
				}
				HBond->donor = duplicate;
				delete donor;
			}

			duplicate = AlreadyRecordedAtom(*atom, hydrogen);
			if ( duplicate == NULL )
			{
				if ( alloc_vector(hydrogen, 0.0, TrjIdx+1) )
				{
					hydrogen->x.at(TrjIdx) = hydrogenX;
					hydrogen->y.at(TrjIdx) = hydrogenY;
					hydrogen->z.at(TrjIdx) = hydrogenZ;
				}

				atom->push_back(hydrogen);
				HBond->hydrogen = hydrogen;
			}
			else
			{
				if ( alloc_vector(duplicate, 0.0, TrjIdx+1) )
				{
					duplicate->x.at(TrjIdx) = hydrogenX;
					duplicate->y.at(TrjIdx) = hydrogenY;
					duplicate->z.at(TrjIdx) = hydrogenZ;
				}

				HBond->hydrogen = duplicate;
				delete hydrogen;
			}

			duplicate = AlreadyRecordedAtom(*atom, acceptor);
			if ( duplicate == NULL )
			{
				if ( alloc_vector(acceptor, 0.0, TrjIdx+1) )
				{
					acceptor->x.at(TrjIdx) = acceptorX;
					acceptor->y.at(TrjIdx) = acceptorY;
					acceptor->z.at(TrjIdx) = acceptorZ;
				}
				atom->push_back(acceptor);
				HBond->acceptor = acceptor;
			}
			else
			{
				if ( alloc_vector(duplicate, 0.0, TrjIdx+1) )
				{
					duplicate->x.at(TrjIdx) = acceptorX;
					duplicate->y.at(TrjIdx) = acceptorY;
					duplicate->z.at(TrjIdx) = acceptorZ;
				}

				HBond->acceptor = duplicate;
				delete acceptor;
			}

			// Save the hydrogen bond to the vector.
			HBond->TrajIdx = TrjIdx;
			hb->push_back(HBond);
		}
		else // Not an Hbond line.
		{
			// Check if this line is another PBC line.
			double CellX, CellY, CellZ;
			double CellAlpha, CellBeta, CellGamma;

			n = sscanf( line,
			            "%le %le %le %le %le %le",
			            &CellX, &CellY, &CellZ,
			            &CellAlpha, &CellBeta, &CellGamma);
			if ( n != 6 ) // Unknown line, abort.
			{
			std::cerr << "Error on line " << linen;
			std::cerr << ". Expected either 6, or 17 arguments, read " <<n<< "."<<std::endl;
			std::cerr << " L" << linen << ": '" << line << "'" << std::endl;
			return(1);
			}

			if ( alloc_vector(Cell, 0.0, TrjIdx+1) )
			{
				// Start another frame of the trajectory.
				if (  (TrjIdx+1)%50==0 )
					std::cout << "\tframe " << TrjIdx+1 << "\r" << std::flush;

				TrjIdx++;

				if ( alloc_vector(Cell, 0.0, TrjIdx+1) )
				{
					Cell->x.at(TrjIdx) = CellX;
					Cell->y.at(TrjIdx) = CellY;
					Cell->z.at(TrjIdx) = CellZ;

					Cell->alpha.at(TrjIdx) = CellAlpha;
					Cell->beta.at(TrjIdx)  = CellBeta;
					Cell->gamma.at(TrjIdx) = CellGamma;
				}
			}
			delete donor;
			delete hydrogen;
			delete acceptor;
			delete HBond;
		}
		in.getline(line,255);
	}
	ifp.close();

	return(TrjIdx+1);
}

struct thbAtom *AlreadyRecordedAtom( std::vector<struct thbAtom *> atom,
                                     struct thbAtom *find)
{
	std::vector<struct thbAtom *>::iterator ai;

	for(ai=atom.begin(); ai < atom.end(); ++ai)
	{
		if ( ( (*ai)->Name == find->Name ) &&
		     ( (*ai)->Molecule == find->Molecule) )
			return( *ai );
	}
	return(NULL);
}
