#include <vector>
#include <time.h>
#include "MessageDefines.h"
#include "deleteVectorPointers.h"
#include "OpenFile.h"
#include "ReadLAMMPS.h"

#define GETLINE_BUF    1023

bool ReadLAMMPSBonds(boost::iostreams::filtering_stream<boost::iostreams::input> *in,
                     std::vector<struct thbBond *> *bonds,
                     int *lineno) {

	char line[GETLINE_BUF+1];
	int n; // To check number of assigned values in sscanf
	struct thbBond *NewBond;
	unsigned int A_id; // Atom ID number
	unsigned int B_id; // Atom ID number

	if ( in->eof() ) { return(false);}
	in->getline(line,GETLINE_BUF); (*lineno)++;
	// Skip mandatory blank line after the 'Bonds' keyword
	if ( in->eof() ) { return(false);}
	in->getline(line,GETLINE_BUF); (*lineno)++;
	// Make sure line is null terminated.
	line[GETLINE_BUF]='\0';

	while ( ! in->eof() )
	{
		n = sscanf( line, "%*u %*u %u %u", &A_id, &B_id);

		if (n <= 0)
			return(true); // End of Bonds section.

		if (n != 2)
			return(false); // Wrong number of arguments.

		NewBond = new struct thbBond;
		NewBond->A_ID = A_id;
		NewBond->B_ID = B_id;
		// LAMMPS does not store bond order, so set to 1 for now.
		NewBond->Order = 1;

		bonds->push_back(NewBond);

		in->getline(line,GETLINE_BUF); (*lineno)++;
		// Make sure line is null terminated.
		line[GETLINE_BUF]='\0';
	}
	(*lineno)--;

	return(true);
}
bool ReadLAMMPSAtoms(boost::iostreams::filtering_stream<boost::iostreams::input> *in,
                     std::vector<struct thbAtom   *> *atom,
                     int *lineno) {

	char line[GETLINE_BUF+1];
	int n; // To check number of assigned values in sscanf
	struct thbAtom *NewAtom;
	unsigned int id; // Atom ID number
	unsigned int resnum; // Residue number (molecule-ID in LAMMPS).
	unsigned int atomtypeid; // Atom Type ID

	if ( in->eof() ) { return(false);}
	in->getline(line,GETLINE_BUF); (*lineno)++;
	// Skip mandatory blank line after the 'Atoms' keyword
	if ( in->eof() ) { return(false);}
	in->getline(line,GETLINE_BUF); (*lineno)++;
	// Make sure line is null terminated.
	line[GETLINE_BUF]='\0';

	while ( ! in->eof() )
	{
		n = sscanf( line, "%u %u %u", &id, &resnum, &atomtypeid );

		if (n <= 0)
			return(true); // End of Atoms section.

		if (n != 3)
			return(false); // Wrong number of values.

		NewAtom = new struct thbAtom;
		NewAtom->ID = id;
		NewAtom->ResidueNum = resnum;
		NewAtom->AtomTypeID = atomtypeid;

		atom->push_back(NewAtom);

		in->getline(line,GETLINE_BUF); (*lineno)++;
		// Make sure line is null terminated.
		line[GETLINE_BUF]='\0';
	}
	(*lineno)--;

	return(true);
}

bool ReadLAMMPSMasses(boost::iostreams::filtering_stream<boost::iostreams::input> *in,
                      std::vector<struct AtomTypes *> *atomtype,
                      int *lineno) {

	char line[GETLINE_BUF+1];
	int n; // To check number of assigned values in sscanf
	struct AtomTypes *NewMass;
	int id; // Masses ID number
	double mass; // Mass
	char *comment;

	if ( in->eof() ) { return(false);}
	in->getline(line,GETLINE_BUF); (*lineno)++;
	// Skip mandatory blank line after the 'Masses' keyword
	if ( in->eof() ) { return(false);}
	in->getline(line,GETLINE_BUF); (*lineno)++;
	// Make sure line is null terminated.
	line[GETLINE_BUF]='\0';

	while ( ! in->eof() )
	{
		n = sscanf( line, "%u %lf", &id, &mass );
		if (n <= 0)
			return(true); // End of Masses section.

		if (n != 2)
			return(false); // Wrong number of values.

		NewMass = new struct AtomTypes;
		NewMass->id = id;
		NewMass->mass = mass;

		// Check for a comment.
		if ( (comment = strstr(line,"#")) != NULL ){
			comment++; // Skip the '#' character

			// Skip leading whitespace
			for(; isspace(*comment) && *comment != '\0'; comment++ );

			strncpy(NewMass->comment, comment, 1023);
		} else {
			NewMass->comment[0] = '\0';
		}
		atomtype->push_back(NewMass);

		in->getline(line,GETLINE_BUF); (*lineno)++;
		// Make sure line is null terminated.
		line[GETLINE_BUF]='\0';
	}

	(*lineno)--;

	return(true);
}

struct AtomTypes *ReadLAMMPSAtomTypeWithID(std::vector<struct AtomTypes *> *atomtype,
                                 unsigned int ID) {

	for(auto it = std::begin(*atomtype) ; it < std::end(*atomtype); ++it)
		if ( (*it)->id == ID ) return((*it));

	return(NULL);
}

struct thbAtom *ReadLAMMPSAtomWithID(std::vector<struct thbAtom *> *atom,
                           unsigned int ID) {

	// std::vector<struct theAtom *>::iterator iter = atom->begin();
	for(auto it = std::begin(*atom) ; it < std::end(*atom); ++it)
		if ( (*it)->ID == ID ) return((*it));

    return(NULL);
}

bool ReadLAMMPSBondsToMolecule(std::vector<struct thbBond      *> *bonds) {

	DEBUG_MSG("Starting ReadLAMMPSBondsToMolecule()");
	for(auto it = std::begin(*bonds) ; it < std::end(*bonds); ++it) {
		if ( (*it)->A->Molecule == NULL )
		{
			std::cout << "Error: No molecule defined for atom "
				      << (*it)->A->ID
				      << "!\n";
			return(false);
		}

		(*it)->A->Molecule->bonds.push_back(*it);
	}
	DEBUG_MSG("Ending ReadLAMMPSBondsToMolecule()");

    return(true);
}

bool ReadLAMMPSAssignMolecules(std::vector<struct thbAtom      *> *atom,
                     std::vector<struct thbMolecule  *> *molecules,
                     std::vector<struct MoleculeDefs *> *moldefs ) {

	struct thbMolecule *NewMolecule;
	struct thbAtom *A;

	DEBUG_MSG("Starting ReadLAMMPSAssignMolecules()");
	for(auto it = std::begin(*moldefs) ; it < std::end(*moldefs); ++it)
	{
		NewMolecule = new struct thbMolecule;

		NewMolecule->Name = std::string("Molecule") + std::to_string((*it)->id);

		for(unsigned int i=(*it)->atomFirst; i <= (*it)->atomLast; ++i)
		{
			A = ReadLAMMPSAtomWithID(atom, i);
			if ( A == NULL ) {
				std::cerr << "Error, could not find atom with ID = " << i << " when assigning molecules!\n";
				return(false);
			}
			NewMolecule->atoms.push_back(A);
			A->Molecule = NewMolecule;
		}

		molecules->push_back(NewMolecule);
	}

	DEBUG_MSG("Ending ReadLAMMPSAssignMolecules()");
	return(true);
}

bool ReadLAMMPSProcessData(std::vector<struct thbAtom      *> *atom,
                       std::vector<struct thbBond      *> *bonds,
                       std::vector<struct AtomTypes    *> *atomtype) {

	for(auto it = std::begin(*atom) ; it < std::end(*atom); ++it)
	{
		(*it)->Name = "Atom" + std::to_string((*it)->ID);
		(*it)->Residue = std::string("Residue");
		(*it)->ForceField = ReadLAMMPSAtomTypeWithID(atomtype, (*it)->AtomTypeID)->comment;

	}

	for(auto it = std::begin(*bonds); it < std::end(*bonds); ++it) {
		struct thbAtom *A = ReadLAMMPSAtomWithID(atom,(*it)->A_ID);
		struct thbAtom *B = ReadLAMMPSAtomWithID(atom,(*it)->B_ID);

		(*it)->A = A;
		(*it)->B = B;

		if ( A == NULL ) {
			std::cerr << "Error, could not find atom with ID = " << (*it)->A_ID << "\n";
			return(false);
		}
		if ( B == NULL ) {
			std::cerr << "Error, could not find atom with ID = " << (*it)->B_ID << "\n";
			return(false);
		}

		// If ForceField of either atom of this bond contains a '=', then
		// assume this bond is a second order.
		if ( (A->ForceField.find("=") != std::string::npos) ||
			 (B->ForceField.find("=") != std::string::npos) ) {
			(*it)->Order = 2;
		}

		A->Bonds.push_back(*it);
		B->Bonds.push_back(*it);

		A->ConnectedAtom.push_back(B);
		B->ConnectedAtom.push_back(A);

		A->ConnectedAtomName.push_back(B->Name);
		B->ConnectedAtomName.push_back(A->Name);

		A->ConnectedAtomResidue.push_back(B->Residue);
		B->ConnectedAtomResidue.push_back(A->Residue);

		A->ConnectedAtomResidueNum.push_back(B->ResidueNum);
		B->ConnectedAtomResidueNum.push_back(A->ResidueNum);

		A->ConnectedAtomBondOrder.push_back((*it)->Order);
		B->ConnectedAtomBondOrder.push_back((*it)->Order);
	}

	/*
	 * for(auto it = std::begin(*bonds) ; it < std::end(*bonds); ++it)
	 *     std::cout << ": "
	 *               << (*it)->A->Name
	 *               << " ("
	 *               << (*it)->A->ForceField
	 *               << ") -- "
	 *               << (*it)->B->Name
	 *               << " ("
	 *               << (*it)->B->ForceField
	 *               << ")"
	 *               << "\n";
	 */

    return(true);
}

bool ReadLAMMPSData(boost::iostreams::filtering_stream<boost::iostreams::input> *in,
                    std::vector<struct thbAtom      *> *atom,
                    std::vector<struct thbBond      *> *bonds) {

	char line[GETLINE_BUF+1];
	int n; // To check number of assigned values in sscanf
	int lineno=0; // Line number of datafile, used for error messages.
	unsigned int val; // temporary value to be stored later.
	std::vector<struct AtomTypes *> atomtype;

	if ( in->eof() ) { return(false);}
	in->getline(line,GETLINE_BUF); lineno++;
	// Make sure line is null terminated.
	line[GETLINE_BUF]='\0';

	while ( ! in->eof() )
	{
		if ( strstr(line, "atoms") != NULL ) {
			n = sscanf( line, "%u", &val );
			if (n != 1) {
				std::cerr << "Error on line " << lineno << ". Expected integer for # of atoms.\n";
				return(false);
			}
			std::cout << "Atoms: " << line << "\n";
			atom->reserve(val);
		} else if ( strstr(line, "bonds") != NULL ) {
			n = sscanf( line, "%u", &val );
			if (n != 1) {
				std::cerr << "Error on line " << lineno << ". Expected integer for # of bonds\n";
				return(false);
			}
			std::cout << "Bonds: " << line << "\n";
			bonds->reserve(val);

		} else if ( strstr(line, "atom types") != NULL ) {
			n = sscanf( line, "%u", &val );
			if (n != 1) {
				std::cerr << "Error on line " << lineno << ". Expected integer for # of bond types\n";
				return(false);
			}
			std::cout << "Atom types: " << line << "\n";
			atomtype.reserve(val);
			DEBUG_MSG("Capacity/size of atomtypes: " << atomtype.capacity() << "/" << atomtype.size());
		} else if ( strstr(line, "Masses") != NULL ) {
			DEBUG_MSG("Masses section.");
			if (!ReadLAMMPSMasses( in, &atomtype, &lineno )) {
				std::cerr << "Error on line " << lineno << ". Expected 2 values\n";
				DeleteVectorPointers( atomtype ); atomtype.clear();
				return(false);
			}
			DEBUG_MSG("Capacity/size of atomtypes: " << atomtype.capacity() << "/" << atomtype.size());
		} else if ( strstr(line, "Atoms") != NULL ) {
			DEBUG_MSG("Atoms section.");
			if ( strstr(line,"# full") == NULL ) {
				std::cerr << "Error: Data file needs to be in full atom style.\n";
				DeleteVectorPointers( atomtype ); atomtype.clear();
				return(false);
			}
			if (!ReadLAMMPSAtoms(in, atom, &lineno)) {
				std::cerr << "Error on line " << lineno << ". Expected 3 values\n";
				return(false);
			}
		} else if ( strstr(line, "Bonds") != NULL ) {
			DEBUG_MSG("Bonds section.");
			if (!ReadLAMMPSBonds( in, bonds, &lineno)) {
				std::cerr << "Error on line " << lineno << ". Expected 2 values\n";
				DeleteVectorPointers( atomtype ); atomtype.clear();
				return(false);
			}
		}

		in->getline(line,GETLINE_BUF); lineno++;
		// Make sure line is null terminated.
		line[GETLINE_BUF]='\0';
	}

	DEBUG_MSG("Read " << lineno-1 << " lines.");

	if ( !ReadLAMMPSProcessData( atom, bonds, &atomtype)) {
		std::cerr << "Error: Can not process LAMMPS data.\n";
		DeleteVectorPointers( atomtype ); atomtype.clear();
		return(false);
	}

	DeleteVectorPointers( atomtype ); atomtype.clear();

    return(true);
}

bool ReadLAMMPSMolDefs(boost::iostreams::filtering_stream<boost::iostreams::input> *in,
                 std::vector<struct MoleculeDefs *> *moldefs) {

	char line[GETLINE_BUF+1];
	int n; // To check number of assigned values in sscanf
	int lineno=0; // Line number of datafile, used for error messages.
	unsigned int molID, atomA, atomB;
	struct MoleculeDefs *NewMolDef;

	if ( in->eof() ) { return(false);}
	in->getline(line,82); lineno++;
	// Make sure line is null terminated.
	line[GETLINE_BUF]='\0';

	while ( ! in->eof() )
	{
		if ( line[0] == '#' )
		{}//	std::cout << "# Found a comment." << "\n";
		else {
			n = sscanf( line, "%u %u %u", &molID, &atomA, &atomB );

			if ( n != 3 )
			{
				std::cerr << "Error on line " << lineno << ". Expected 3 arguments, read "
				          << n << ".\n";
				return(false);
			}
			NewMolDef = new struct MoleculeDefs;
			NewMolDef->id = molID;
			NewMolDef->atomFirst = atomA;
			NewMolDef->atomLast = atomB;

			moldefs->push_back(NewMolDef);

			// std::cout << "L: " << NewMolDef->id << " / " << atomA << " / " << atomB << "\n";
		}
		in->getline(line,82); lineno++;
		// Make sure line is null terminated.
		line[GETLINE_BUF]='\0';
	}

    return(true);
}

bool ReadLAMMPSFrameCoordinates(boost::iostreams::filtering_stream<boost::iostreams::input> *in,
                                std::vector<struct thbAtom   *> *atom,
                                std::vector<Point> *Coordinates,
                                bool SaveMemory,
                                class Point CellDimensions,
                                class Point Offset,
                                unsigned int N,
                                unsigned int *lineno) {

	char line[GETLINE_BUF+1];
	unsigned int n;
	unsigned int ID;
	double x, y, z;
	Point theOffset = Offset;

	bool isScaled = false;

	// If CellDimensions is not all zeros, then we must be reading a file
	// with scaled coordinates. Also, use a 0,0,0 offset when reading
	// scaled coordinates.
	if ( (CellDimensions.x() != 0.0) &&
	     (CellDimensions.y() != 0.0) &&
	     (CellDimensions.z() != 0.0) )
	{
		isScaled = true;
		theOffset = Point(0,0,0);
	}

	for (unsigned int i=0; (i < N) && !in->eof(); ++i)
	{
		in->getline(line,GETLINE_BUF); (*lineno)++;
		// Make sure line is null terminated.
		line[GETLINE_BUF]='\0';

		n = sscanf(line, "%u %*u %lf %lf %lf", &ID, &x, &y, &z);
		if (n != 4)
		{
			std::cerr << "Error on line " << *lineno << ". Expected at least 5 numbers.\n";
			return(false);
		}

		// Convert to real coordinate by multiplying by the cell dimensions.
		if ( isScaled ) {
			x = x*CellDimensions.x();
			y = y*CellDimensions.y();
			z = z*CellDimensions.z();
		}

		if ( SaveMemory )
			Coordinates->push_back( Point(x,y,z) - theOffset );
		else
			ReadLAMMPSAtomWithID(atom, ID)->p.push_back(Point(x,y,z) - theOffset);
	}

	return(true);
}
bool ReadLAMMPSFrame(boost::iostreams::filtering_stream<boost::iostreams::input> *in,
                          std::vector<struct thbAtom   *> *atom,
                          struct PBC *Cell,
                          std::vector<Point> *Coordinates,
                          bool SaveMemory) {

	char line[GETLINE_BUF+1];
	int n; // To check number of assigned values in sscanf
	double Xmin, Xmax;
	double Ymin, Ymax;
	double Zmin, Zmax;
	class Point PointOffset;
	unsigned int lineno = 0;
	unsigned int NumberOfAtoms = 0;
	unsigned int timestep = 0;

	if ( in->eof() ) { return(false);}
	in->getline(line,GETLINE_BUF); lineno++;
	// Make sure line is null terminated.
	line[GETLINE_BUF]='\0';


	while ( ! in->eof() )
	{
		if( !strncmp(line, "ITEM: TIMESTEP", 14) ) {
			in->getline(line,GETLINE_BUF); lineno++;
			// Make sure line is null terminated.
			line[GETLINE_BUF]='\0';
			n = sscanf(line, "%u", &timestep);
			if (n != 1)
			{
				std::cerr << "Error on line " << lineno
				          << ". Expected 1 number.";
				return(false);
			}
		}
		else if( !strncmp(line, "ITEM: NUMBER OF ATOMS", 21) ) {
			in->getline(line,GETLINE_BUF); lineno++;
			// Make sure line is null terminated.
			line[GETLINE_BUF]='\0';
			n = sscanf(line, "%u", &NumberOfAtoms);
			if (n != 1)
			{
				std::cerr << "Error on line " << lineno
				          << ". Expected 1 number.";
				return(false);
			}
		}
		else if( !strncmp(line, "ITEM: BOX BOUNDS", 16) ) {
			in->getline(line,GETLINE_BUF); lineno++;
			// Make sure line is null terminated.
			line[GETLINE_BUF]='\0';
			n = sscanf(line, "%lf %lf", &Xmin, &Xmax);
			if (n != 2)
			{
				std::cerr << "Error on line " << lineno
				          << ". Expected 2 numbers.";
				return(false);
			}

			in->getline(line,GETLINE_BUF); lineno++;
			// Make sure line is null terminated.
			line[GETLINE_BUF]='\0';
			n = sscanf(line, "%lf %lf", &Ymin, &Ymax);
			if (n != 2)
			{
				std::cerr << "Error on line " << lineno
				          << ". Expected 2 numbers.";
				return(false);
			}

			in->getline(line,GETLINE_BUF); lineno++;
			// Make sure line is null terminated.
			line[GETLINE_BUF]='\0';
			n = sscanf(line, "%lf %lf", &Zmin, &Zmax);
			if (n != 2)
			{
				std::cerr << "Error on line " << lineno
				          << ". Expected 2 numbers.\n";
				return(false);
			}

			// This will be subtracted from all coordinates.
			PointOffset = Point(Xmin, Ymin, Zmin);

			// Shift PBC box to start at origin.
			Cell->p.push_back(Point(Xmax,Ymax,Zmax) - PointOffset);
			Cell->angles.push_back(Point(90.0, 90.0, 90.0));
			Cell->frames++;

		} else if (!strncmp(line, "ITEM: ATOMS", 11)) {

			if (!strncmp(line, "ITEM: ATOMS id mol x y z", 24)) {
				return(ReadLAMMPSFrameCoordinates(in,
				                                  atom,
				                                  Coordinates,
				                                  SaveMemory,
				                                  Point(0,0,0),
				                                  PointOffset,
				                                  NumberOfAtoms,
				                                  &lineno));
			} else if (!strncmp(line, "ITEM: ATOMS id mol xs ys zs", 24)) {
				return(ReadLAMMPSFrameCoordinates(in,
				                                  atom,
				                                  Coordinates,
				                                  SaveMemory,
				                                  Cell->p.back(),
				                                  PointOffset,
				                                  NumberOfAtoms,
				                                  &lineno));
			} else {
				std::cerr << "Error on line " << lineno
				          << ". Expected 'id mol x y z'.\n";
				return(false);
			}
		}

		in->getline(line,GETLINE_BUF); lineno++;
		// Make sure line is null terminated.
		line[GETLINE_BUF]='\0';
	}

	// Nothing left to read in this file, so signal done with a return false.
	return(false);
}


bool ReadLAMMPSConnections( char *fileData,
                            char *fileMols,
                            std::vector<struct thbAtom     *> *atom,
                            std::vector<struct thbMolecule *> *molecules,
                            std::vector<struct thbBond     *> *bonds ) {


	std::vector<struct MoleculeDefs *> moldefs;
	struct PBC Cell;

	// Open LAMMPS molecular definitions file.
	std::ifstream ifpMols;
	boost::iostreams::filtering_stream<boost::iostreams::input> inMols;

	if ( !openfile(fileMols, &inMols, &ifpMols) ) {
		return(false);
	} else {
		VERBOSE_MSG("Reading molecule configuration from " << fileMols);
		// ReadMdf( &MDFin, atom, molecules );
		if ( !ReadLAMMPSMolDefs(&inMols, &moldefs) ) {
			ifpMols.close();
			return(false);
		}
		VERBOSE_MSG("Read " << moldefs.size() << " molecule configuations");
		ifpMols.close();
	}

	DEBUG_MSG("Capacity/size of moldefs: " << moldefs.capacity() << "/" << moldefs.size());

	// Open LAMMPS data file.
	std::ifstream ifpData;
	boost::iostreams::filtering_stream<boost::iostreams::input> inData;

	if ( !openfile(fileData, &inData, &ifpData) ) {
		DeleteVectorPointers( moldefs ); moldefs.clear();
		return(false);
	} else {
		VERBOSE_MSG("Reading LAMMPS data file from " << fileData);

		if ( !ReadLAMMPSData(&inData, atom, bonds) ) {
			ifpData.close();
			// Do not need moldefs anymore.
			DeleteVectorPointers( moldefs ); moldefs.clear();
			return(false);
		}

		ifpData.close();

		if ( !ReadLAMMPSAssignMolecules(atom, molecules, &moldefs ) ) {
			DeleteVectorPointers( moldefs ); moldefs.clear();
			return(false);
		}
		if ( !ReadLAMMPSBondsToMolecule(bonds) ) {
			DeleteVectorPointers( moldefs ); moldefs.clear();
			return(false);
		}

	}

	DEBUG_MSG("Capacity/size of atoms: " << atom->capacity() << "/" << atom->size());
	DEBUG_MSG("Capacity/size of bonds: " << bonds->capacity() << "/" << bonds->size());


	DeleteVectorPointers( moldefs ); moldefs.clear();

	return(true);
}

bool ReadLAMMPSPositions(const char *fileTrj,
                         std::vector<struct thbAtom *> *atom,
                         struct PBC *Cell,
                         std::vector<struct thbAtom *> *hydrogens,
                         std::vector<struct thbAtom *> *acceptors,
                         double rCutoff, double angleCutoff, bool SaveMemory) {

	// Open LAMMPS trajectory file
	std::ifstream ifpTrj;
	boost::iostreams::filtering_stream<boost::iostreams::input> inTrj;

	if ( !openfile(fileTrj, &inTrj, &ifpTrj) ) {
		return(false);
	}

	VERBOSE_MSG("Reading atom coordinates from " << fileTrj);

	time_t timer = time(NULL);

	while ( 1 ) {

		struct worker_data_s wd;
		if ( SaveMemory ) {
			wd.coordinates = new std::vector<Point>;
			wd.coordinates->reserve(atom->size());
		} else {
			wd.coordinates = NULL;
		}

		/** Read next frame */
		if ( !ReadLAMMPSFrame(&inTrj, atom, Cell, wd.coordinates, SaveMemory) ) {
			if ( SaveMemory ) { delete wd.coordinates; }
			break;
		}

		// Do not update message more often then once per second.
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

	ifpTrj.close();

	return(true);
}

