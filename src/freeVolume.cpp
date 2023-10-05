#include <random>
#include "freeVolume.h"
#include "MessageDefines.h"

// #ifdef PTHREADS
// extern Queue<struct worker_data_s> inQueue;
// extern Queue<struct worker_data_s> outQueue;
// #endif

// Calculate fractional free volume.
void freeVolume(std::vector<thbAtom *> *atom,
                Point cell,
                unsigned int TrjIdx)
{
	// Will be used to obtain a seed for the random number engine
	std::random_device rd;

	// Standard mersenne_twister_engine seeded with rd()
	std::mt19937 gen(rd());

	std::uniform_real_distribution<> disX(0, cell.x());
	std::uniform_real_distribution<> disY(0, cell.y());
	std::uniform_real_distribution<> disZ(0, cell.z());

	// std::cout << cell.x() << ", "
	//           << cell.y() << ", "
	//           << cell.z() << ".\n";


	// Repeatedly pick a random location and see if it is considered
	// either occupied or free volume
	unsigned int throwsPerCheck = 10000;
	unsigned int throwsTotal = 0;
	unsigned int freeVolumeCount = 0;
	double freeVolumePrevious = 2.0;
	bool goodEnough = false;
	
	// For the PIB file I am currently using for testing.
	std::vector<double>rVDW = {
	                           2.0595/2.0,
	                           2.9816/2.0,
	                           2.9816/2.0,
	                           2.9816/2.0,
	                           2.9816/2.0,
	                          };

	// Radius of free-volume probe.
	double rProbe = 0.75;

	// Distance from particle center, plus probe radius, that is considered
	// occupied.
	std::vector<double>rOccupiedSquared;



	
	// std::vector<double>rVDW = { 0.5, 0.5, 0.5, 0.5, 0.5};
	// For dendrimer used for testing.
	// std::vector<double>rVDW = {
	//                            2.7685,
	//                            2.8124,
	//                            2.3201,
	//                            1.7744,
	//                            2.8697,
	//                            2.8124,
	//                            2.3201,
	//                            2.3201,
	//                            0,
	//                           };

	// Build rOccupiedSquared vector
	for(auto it = std::begin(rVDW) ; it < std::end(rVDW); ++it)
		rOccupiedSquared.push_back( ((*it)+rProbe)*((*it)+rProbe) );

	while ( !goodEnough ) {
	for(unsigned int i = 0; i < throwsPerCheck; i++)
	{
		// Each call to dis(gen) generates a new random double.
		Point rXYZ(disX(gen), disY(gen), disZ(gen));
		bool occupiedSite = false;
		double dSquared;
		for(auto it = std::begin(*atom) ; it < std::end(*atom); ++it)
		{
			dSquared=(*it)->p.at(TrjIdx).minimumImage(rXYZ, cell).magnitudeSquared();
 
			// std::cout << "AtomTypeID: " << (*it)->AtomTypeID << '\n';
			if ( (*it)->AtomTypeID > rOccupiedSquared.size() )
			{
				std::cerr << "Error: Atom Type ID " << (*it)->AtomTypeID
				          << ", in frame " << TrjIdx+1 << ","
				          << " not in rVWD vector.\n";
				return;
			}

			if ( dSquared < rOccupiedSquared.at((*it)->AtomTypeID-1) )
			{
				occupiedSite = true;
				// Found an atom with a specific distance,
				// so this site is considered occupied. Break out of
				// this loop since it is not a free volume site.
				break;
			}
		}
		if ( !occupiedSite )
			freeVolumeCount++;
	}
	throwsTotal += throwsPerCheck;

	double freeVolumeCurrent = freeVolumeCount/(double)throwsTotal;

	double fracChange = (freeVolumeCurrent - freeVolumePrevious)/freeVolumePrevious;

	// std::cout << throwsTotal << " - " << fracChange*100.0 << "\n";

	if ( fabs(fracChange) < 0.01 )
		goodEnough = true;

	freeVolumePrevious = freeVolumeCurrent;

	}

	// Frame 1: Free-volume found 1530 times out of 20000 samples -- 7.65%
	std::cout << "Frame " << TrjIdx+1 << ": "
	          << "Free-volume found "
	          << freeVolumeCount << " times out of "
	          << throwsTotal << " samples -- "
	          << 100.0*freeVolumeCount/(double)throwsTotal << "%\n";
	// T/C = vbox/vfree
	// vfree/vobx = C/T
}
