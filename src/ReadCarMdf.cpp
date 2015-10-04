#include "MagicNumber.h"
#include "ReadCarMdf.h"
#include <time.h>
#ifdef USE_LZMA
#include "lzma.h"
#endif

extern bool THB_VERBOSE;

#ifdef PTHREADS
extern Queue<struct worker_data_s> inQueue;
extern Queue<struct worker_data_s> outQueue;
#endif

unsigned int toUInt(std::string s)
{
	std::stringstream in(s, std::ios_base::in);
	unsigned int result;

	in >> result;

	return result;
};

float toFloat(std::string s)
{
	std::stringstream in(s, std::ios_base::in);
	float result;

	in >> result;


	return result;
};

std::vector< std::string >connectedAtoms( const char *connected)
{
	// First split string by spaces.
	std::stringstream ss(connected);

	// Copy from the stream into a vector.
	std::vector< std::string >tokens;
	std::copy(std::istream_iterator< std::string >(ss),
              std::istream_iterator< std::string >(),
              std::back_inserter(tokens) );

	// Each element of vector can have a string like:
	// 1) xxxx_XX:YYYY/zzzz
	// 2) YYYY/zzzz
	// 3) xxxx_xx:YYYY
	// 4) YYYY
	//
	// xxxx_XX, and zzzz are optional
	//
	// The connection vectors with be in groups of 4 for the xxxx ("" if not
	// present), _XX (0 if not present), YYYY (always present), and zzzz (1.0
	// is no present)
	std::vector< std::string > connections;
	for( unsigned int i=0; i< tokens.size(); ++i)
	{
		std::string residue="";
		std::string residueNum="0";
		std::string atomName="";
		std::string bondOrder="1.0";

		std::size_t foundSlash      = tokens.at(i).find("/");
		std::string remaining1;
		if( foundSlash != std::string::npos )
		{
			bondOrder.assign( tokens.at(i), foundSlash+1, std::string::npos );
			remaining1.assign( tokens.at(i), 0, foundSlash );
		}
		else
			remaining1.assign( tokens.at(i) );

		// remaining1 can have one of the following formats:
		// 1) xxxx_XX:YYYY
		// 2) YYYY
		std::size_t foundColon      = remaining1.find(":");
		std::string remaining2 = "";
		if( foundColon != std::string::npos )
		{
			atomName.assign( tokens.at(i), foundColon, std::string::npos);
			remaining2.assign( remaining1, 0, foundColon );
		}
		else
		{
			atomName = remaining1;
			remaining2 = "";
		}

		// remaining2 can have one of the following formats:
		// 1) xxxx_XX
		// 2) <empty string>
		std::size_t foundUnderscore = tokens.at(i).find("_");
		if ( foundUnderscore != std::string::npos )
		{
			residue.assign(remaining2, 0, foundUnderscore );
			residueNum.assign(remaining2, foundUnderscore+1, std::string::npos);
		}


		connections.push_back(residue);
		connections.push_back(residueNum);
		connections.push_back(atomName);
		connections.push_back(bondOrder);
	}

	return(connections);
}

bool openfile( const char *filename,
               boost::iostreams::filtering_stream<boost::iostreams::input> *in,
               std::ifstream *ifp)
{
	int magicNum;

	magicNum = getMagicNumber(filename);

	if ( magicNum == -1 )
	{
		std::cout << "Error opening " << filename << "\n";
		return(false);
	}


#ifdef USE_ZLIB
	if ( magicNum == MAGICNUMBER_GZIP )
	{
		in->push(boost::iostreams::gzip_decompressor());
		ifp->open(filename,std::ios::in|std::ios::binary);
	}
#endif
#ifdef USE_BZIP2
	if ( magicNum == MAGICNUMBER_BZIP2 )
	{
		in->push(boost::iostreams::bzip2_decompressor());
		ifp->open(filename,std::ios::in|std::ios::binary);
	}
#endif
#ifdef USE_LZMA
	if ( magicNum == MAGICNUMBER_LZMA )
	{
		in->push(lzma_input_filter());
		ifp->open(filename,std::ios::in|std::ios::binary);
	}
#endif

	if ( magicNum == MAGICNUMBER_UNKNOWN )
		ifp->open(filename,std::ios::in);

	if ( !ifp->is_open() )
	{
		std::cerr << "Error: can not open this type of file." << "\n";
		return(false);
	}

	in->push(*ifp);
	return(true);
}

bool ReadMdf(boost::iostreams::filtering_stream<boost::iostreams::input> *in,
             std::vector<struct thbAtom *> *atom)
{
	struct thbAtom *NewAtom;
	int lineno=0; // Line number of datafile, used for error messages.

	char line[256];
	char molecule[80];
	char residue[4];
	unsigned int residueNum;
	unsigned int ID=0;
	char name[80];
	char type[4];
	char forcefield[4];
	char connected[256];


	in->getline(line,255); lineno++;
	DEBUG_MSG("line[0]: " << line);
	bool ReadingMolecule = false;
	while ( !in->eof() )
	{
		if( (line[0] == '!') ||
		    (line[0] == ' ') )
		{
			in->getline(line,255); lineno++;
			continue;
		}

		if( line[0] == '@' ) //A Parameter line.
		{
			//Check for '@molecule' parameter.
			if ( !strncmp(line,"@molecule ", 10) )
			{
				ReadingMolecule = true;
				if( sscanf(line, "%*s %s", molecule) != 1 )
				{
					std::cerr << "Error datermining molecule name!" << "\n";
					return(false);
				}
			}
			else
				ReadingMolecule = false;
		}
		else if( (line[0] != '#') && ReadingMolecule ) //Must be an atom line.
		{
			int n = sscanf( line,
			                "%4s_%d:%12s %2s %3s %*s %*s %*s %*s %*s %*s %*s %*s %*s %[^\n]",
			                residue, &residueNum, name, type, forcefield,
			                connected);

			if ( n != 6 )
			{
				std::cerr << "Error. Expected 6 arguments, read "
				          << n << "." << "\n"
						  << " L" << lineno <<": '" << line << "'" << "\n";

				// ifp.close();
				return(1);
			}
			std::vector<std::string> cA = connectedAtoms(connected);

			NewAtom = new struct thbAtom;
			NewAtom->ID         = ID; ID++;
			NewAtom->Name       = name;
			NewAtom->Type       = type;
			NewAtom->Molecule   = molecule;
			NewAtom->Residue    = residue;
			NewAtom->ResidueNum = residueNum;
			NewAtom->ForceField = forcefield;

			for (unsigned int i=0; i< cA.size(); i=i+4)
			{
				NewAtom->ConnectedAtomBondOrder.push_back (toFloat(cA.at(i+3)));
				NewAtom->ConnectedAtomName.push_back      ( cA.at(i+2));

				if( cA.at(i) == "" )
				{ // Its the same as the main atom.
					NewAtom->ConnectedAtomResidue.push_back(NewAtom->Residue);
					NewAtom->ConnectedAtomResidueNum.push_back(NewAtom->ResidueNum);
				}
				else
				{
					NewAtom->ConnectedAtomResidue.push_back( cA.at(i) );
					NewAtom->ConnectedAtomResidueNum.push_back(toUInt(cA.at(i+1)));
				}
			}
			// Save the NewAtom to the list.
			atom->push_back(NewAtom);

		}

		in->getline(line,255); lineno++;
	}

	return(true);
}
// Is b connected to a
bool isConnectedAtom( struct thbAtom *a, struct thbAtom *b, unsigned int i)
{


	// Compare the first characters of Name initially, to speed this up a
	// little bit.
	if( (a->ConnectedAtomName.at(i).at(0) == b->Name[0] ) &&
	    (a->ConnectedAtomName.at(i)       == b->Name ) &&
	    (a->Molecule                      == b->Molecule) &&
	    (a->ConnectedAtomResidue.at(i)    == b->Residue) &&
	    (a->ConnectedAtomResidueNum.at(i) == b->ResidueNum) )
		return(true);

	return(false);
}
void doAtomConnections( std::vector<struct thbAtom *> *atom )
{
	std::vector<struct thbAtom *>::iterator it_a;
	std::vector<struct thbAtom *>::iterator it_b;

	// All atoms.
	for(it_a = atom->begin() ; it_a < atom->end(); ++it_a)
	{
		// All connections of atom *it_a.
		for(int i=(*it_a)->ConnectedAtomName.size()-1; i >=0 ; --i)
		{
			// Search all atoms for a match to connection.
			for(it_b = atom->begin(); it_b < atom->end(); ++it_b)
			{
				if ( isConnectedAtom( *it_a, *it_b, i ) )
				{
					// Connect these two atoms.
					(*it_a)->ConnectedAtom.push_back(*it_b);
					(*it_b)->ConnectedAtom.push_back(*it_a);

					// We don't need these anymore.
					(*it_a)->ConnectedAtomName.pop_back();
					(*it_a)->ConnectedAtomResidue.pop_back();
					(*it_a)->ConnectedAtomResidueNum.pop_back();
					(*it_a)->ConnectedAtomBondOrder.pop_back();

					// Need to find appropriate element to delete from *it_b->
					// ConnectedAtom* vectors. I am losing the BondOrder value
					// with this methond, but I don't need it.
					int j = (*it_b)->ConnectedAtomName.size()-1;
					for(; j>=0; --j)
					{
						if ( isConnectedAtom( *it_b, *it_a, j ) )
						{
							(*it_b)->ConnectedAtomName.erase(
							        (*it_b)->ConnectedAtomName.begin()+j);
							(*it_b)->ConnectedAtomResidue.erase(
							        (*it_b)->ConnectedAtomResidue.begin()+j);
							(*it_b)->ConnectedAtomResidueNum.erase(
							        (*it_b)->ConnectedAtomResidueNum.begin()+j);
							(*it_b)->ConnectedAtomBondOrder.erase(
							        (*it_b)->ConnectedAtomBondOrder.begin()+j);
							break;
						}
					}
					break;
				}
			}
		}
	}
}

bool ConnectionsMDF(const char *filename,
                    std::vector<struct thbAtom *> *atom)
{
	//
	// Read atoms from MDF files.
	//

	std::string MDFfile = filename;

	// The user may specify either a .arc, or .car file. They both use a .mdf
	// file for connections.
	size_t tag = MDFfile.rfind(".arc");
	if ( tag == std::string::npos )
		tag = MDFfile.rfind(".car");

	if ( tag != std::string::npos )
		MDFfile.replace(tag, 4, ".mdf");

	std::ifstream MDFifp;
	if (MDFifp == NULL)
		return(false);
    boost::iostreams::filtering_stream<boost::iostreams::input> MDFin;

    if ( !openfile( MDFfile.c_str(), &MDFin, &MDFifp ) ) {
	    return(false); }
    else {
		VERBOSE_MSG("Reading atom configuration from " << MDFfile);
		ReadMdf( &MDFin, atom );
		MDFifp.close();
    }

	// Free up some space, if possible.
	if ( atom->size() != atom->capacity() )
		std::vector<struct thbAtom *>(*atom).swap(*atom);
	DEBUG_MSG("Capacity/size of atom: " << atom->capacity() << "/" << atom->size());

	VERBOSE_MSG("Determining atom connections.");
	doAtomConnections( atom );

	return(true);
}

bool PositionsCAR(const char *filename,
                  unsigned int everyNth,
                  std::vector<struct thbAtom *> *atom,
                  struct PBC *Cell,
                  std::vector<struct thbAtom *> *hydrogens,
                  std::vector<struct thbAtom *> *acceptors,
                  double rCutoff, double angleCutoff, bool SaveMemory)
{
	// Read coordinates from CAR files.
	std::string CARfile = filename;

	std::ifstream CARifp;
	if (CARifp == NULL)
		return(false);
    boost::iostreams::filtering_stream<boost::iostreams::input> CARin;

	if ( !openfile( CARfile.c_str(), &CARin, &CARifp ) ) {
		return(false); }
	else {
		VERBOSE_MSG("Reading atom coordinates from " << CARfile);

		time_t timer = time(NULL);

		while ( 1 ) {

			struct worker_data_s wd;
			if ( SaveMemory ) {
				wd.coordinates = new std::vector<Point>;
				wd.coordinates->reserve(atom->size());
			} else {
				wd.coordinates = NULL;
			}

			/** Read CAR, or next frame of ARC. */
			if ( !ReadCar(&CARin, everyNth,
			              atom, Cell, wd.coordinates, SaveMemory) ) {
				if ( SaveMemory ) {
					delete wd.coordinates; }
				break;
			}

			if ( (Cell->framesInFile-1) % everyNth ) {
				// Skip this frame (--every)
				if ( SaveMemory ) {
					delete wd.coordinates; }
				continue;
			}

			// Update the message every second.
			if ( difftime(time(NULL), timer) > 1.0 )
			{
				VERBOSE_RMSG("Frames : " << Cell->frames);
				timer = time(NULL);
			}

			if ( SaveMemory ) {
				wd.jobtype = THREAD_JOB_HBS2;

			} else {
				wd.jobtype = THREAD_JOB_HBS;
			}
			wd.jobnum = Cell->frames;
			wd.num_threads = NumberOfCPUs();
			wd.cell = Cell->p.back();
			wd.hydrogens = hydrogens;
			wd.acceptors = acceptors;
			wd.TrjIdx = Cell->frames - 1;
			wd.rCutoff = rCutoff;
			wd.angleCutoff = angleCutoff;
			wd.hb = new HBVec;
			wd.hb->reserve(acceptors->size()*2);

			inQueue.push(wd);
		}
		CARifp.close();
	}

	return(true);
}

bool ReadCar(boost::iostreams::filtering_stream<boost::iostreams::input> *in,
             unsigned int everyNth,
             std::vector<struct thbAtom *> *atom,
             struct PBC *Cell,
             std::vector<Point> *Coordinates,
             bool SaveMemory)
{
	char line[83];
	char CarEND[] = "end                                                                             ";

	int n; // To check number of assigned values in sscanf
	int lineno=0; // Line number of datafile, used for error messages.

	/*
	 * For each frame, start counting from atomNum=0. car/arc/mdf files
	 * guarantee the the order of the atoms will not change and the sequence of
	 * atoms in the car/arc files are the same as in the mdf file.
	 */
	unsigned int atomNum = 0;

	if ( in->eof() ) { return(false);}
	in->getline(line,82); lineno++;
	if ( in->eof() ) {
		DEBUG_MSG("line[0]: <EOF>" << line);
		return(false);
	}

	// User requested to skip this frame, so read to end of this frame
	// then return.
	// TODO: Fix this every line!
	if ( Cell->framesInFile % everyNth ) {
		while ( ! in->eof() ) {
			if ( *line == *CarEND) {
				in->getline(line,82);
				if ( *line == *CarEND ) {
					Cell->framesInFile++;
					return(true); // at end of this frame
				}
			}
			else
				in->getline(line,82);
		}
		Cell->framesInFile++;
		return(false); // End of file.
	}

	while ( ! in->eof() )
	{
		// VERBOSE_MSG(":" << line);
		if ( line[0] == '!' )
		{}//	std::cout << "# Found a comment." << "\n";
		else if ( ! strncmp(line, "Materials Studio Generated CAR File", 35) )
		{}//	std::cout << "# Found Materials Studio comment." << "\n";
		else if ( *line == *CarEND )
		{
			in->getline(line,82); lineno++;
			// Two 'end' statements in a row indicate a new frame of the
			// trajectory.
			if ( *line == *CarEND ) {
				atomNum = 0;
				// return after reading a full frame.
				return(true);
			}
			continue;
		}
		else if ( ! strncmp(line, "PBC=ON", 6) ) {  }
		else if ( ! strncmp(line, "PBC=OFF", 7) ) {  }
		else if ( ! strncmp(line, "Configurations", 14) ) {}
		else if ( ! strncmp(line, "PBC ", 4) )
		{
			double CellX, CellY, CellZ;
			double CellAlpha, CellBeta, CellGamma;

			n = sscanf( line,
			            "%*s %le %le %le %le %le %le",
			            &CellX, &CellY, &CellZ,
			            &CellAlpha, &CellBeta, &CellGamma);

			if ( n != 6 )
			{
				std::cerr << "Error on line 1. Expected 6 arguments, read "
				          << n << "." << "\n"
				          << " L1: '" << lineno << "'" << "\n";

				// ifp.close();
				return(false);
			}
			Cell->p.push_back( Point(CellX, CellY, CellZ) );

			Cell->angles.push_back( Point(CellAlpha, CellBeta, CellGamma) );
			Cell->frames++;
			Cell->framesInFile++;
		}
		else // Must be a new atom.
		{
			double x;
			double y;
			double z;

			n = sscanf( line, "%*s %le %le %le", &x, &y, &z);
			if ( n != 3 )
			{
				std::cerr << "Error on line 1. Expected 3 coordinates, read "
				          << n << "." << "\n"
				          << " L1: '" << lineno << "'" << "\n";
				// ifp.close();
				return(false);
			}
			if ( SaveMemory )
				Coordinates->push_back( Point(x, y, z) );
			else
				atom->at(atomNum)->p.push_back( Point(x,y,z) );

			// Update atomNum counter.
			atomNum++;
		}
		in->getline(line,82); lineno++;
	}
	// ifp.close();
	// in.close();

	return(false);
}
