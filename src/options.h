#ifndef _options_h
#define _options_h

#include "TraceHBonds.h"

struct useroptions
{
	unsigned int NumBins;
	unsigned int EveryNth; // Load every frame by default
	char *fArc;            // arc filename (without the .arc extension).
	char *ofPrefix;        // Prefix for output filename.
	char *ofSuffix;        // Suffix for output filename.
	double rCutoff;        // Distance Cutoff for Hydrogen bonds.
	double angleCutoff;    // Angle Cutoff for Hydrogen bonds.
	unsigned char flags;
	/**
	 * Forcefield (thbAtom::ForceField) of hydrogens that may contribute to
	 * hydrogen bonding
	 */
	std::vector<std::string>Hydrogens;
	/**
	 * Forcefield  (thbAtom::ForceField) of atoms which may act as hydrogen
	 * bond acceptors
	 */
	std::vector<std::string>Acceptors;
	bool SaveMemory;

	public:
	useroptions() {
		NumBins     = 0;
		EveryNth    = 1;
		fArc        = NULL;
		ofPrefix    = NULL;
		ofSuffix    = NULL;
		rCutoff     =   0.0;
		angleCutoff = 180.0;
		flags       = 0;
		SaveMemory  = true;
	}
};
#endif // _options_h
