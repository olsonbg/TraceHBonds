#include "Structure.h"

void
saveStructure(unsigned int NumFramesInTrajectory,
              char *ofPrefix, char *ofSuffix,
              struct PBC *Cell,
              std::vector<struct thbMolecule *> *molecules)
{
	time_t timer = time(NULL);
	unsigned int TrjIdx;

	std::vector<struct thbMolecule *>::iterator mol_it;
	std::vector<struct thbAtom *>::iterator atom_it;
	std::vector<struct thbBond *>::iterator bond_it;

	OFmt colX(9,4);
	OFmt colY(9,4);
	OFmt colZ(9,4);
	for( TrjIdx = 0 ; TrjIdx != NumFramesInTrajectory; ++TrjIdx )
	{
		std::stringstream ofilename;
		ofilename << ofPrefix << "-Structure" << TrjIdx+1 << ofSuffix;

		std::ofstream out;
		out.open(ofilename.str().c_str(),std::ios::out);
		if ( out.is_open() )
		{
			if ( (difftime(time(NULL), timer) > 1.0) ||
				 (TrjIdx == NumFramesInTrajectory-1) )
			{
				BRIEF_RMSG("Saving chemical structures, frame " << TrjIdx+1 << "/" << NumFramesInTrajectory << " (" << ofilename.str() << ")");
				timer = time(NULL);
			}

			out << "[{\n    \"PBC\": {\n"
			    << "      \"xyz\": ["
			    << Cell->p.at(TrjIdx).x()      << ", "
			    << Cell->p.at(TrjIdx).y()      << ", "
			    << Cell->p.at(TrjIdx).z()      << " "
			    << "],\n"
			    << "      \"angles\": ["
			    << Cell->angles.at(TrjIdx).x() << ", "
			    << Cell->angles.at(TrjIdx).y() << ", "
			    << Cell->angles.at(TrjIdx).z()
			    << "]\n"
			    << "    }\n},";

			out << "\n{   \"molecules\": [\n";
			mol_it = molecules->begin();
			for ( ; mol_it < molecules->end(); ++mol_it)
			{
				out << "    { \"name\": "
				    << "\"" << (*mol_it)->Name << "\",\n"
				    << "      \"numAtoms\": "
				    << (*mol_it)->atoms.size() << ",";

				out << "\n      \"atoms\": [\n";
				atom_it = (*mol_it)->atoms.begin();
				for ( ; atom_it < (*mol_it)->atoms.end(); ++atom_it)
				{
					out << "              { \"location\": [ "
					    << colX << (*atom_it)->p.at(TrjIdx).x() << ", "
					    << colY << (*atom_it)->p.at(TrjIdx).y() << ", "
					    << colZ << (*atom_it)->p.at(TrjIdx).z() << " ],\n"
					    << "                "
					    << "\"name\": \""     << (*atom_it)->Name << "\", "
					    << "\"element\": \""     << (*atom_it)->Type << "\", "
					    << "\"forcefield\": \"" << (*atom_it)->ForceField
					    << "\", \"ID\": " << (*atom_it)->ID << "}";

					if (atom_it != (*mol_it)->atoms.end()-1) {
						out << ",\n";
					}
				}
				out << " ],\n";
				out << "      \"bonds\": [\n";
				bond_it = (*mol_it)->bonds.begin();
				for ( ; bond_it < (*mol_it)->bonds.end(); ++bond_it)
				{
					out << "              { \"IDs\": [ "
					    << (*bond_it)->A->ID << ", "
					    << (*bond_it)->B->ID << "], "
					    << "\"order\": " << (*bond_it)->Order
					    << "}";
					if ( bond_it != (*mol_it)->bonds.end()-1 ) {
						out << ",\n";
					}
				}
				if ( mol_it != molecules->end()-1 ) {
					out << " ] },\n";
				}
			}
			out << "]} ]\n}]\n";
			out.close();
		} else {
			BRIEF_MSG("ERROR: Can not save " <<ofilename.str() << "!" );
		}
	}
	BRIEF_MSG("");
}
