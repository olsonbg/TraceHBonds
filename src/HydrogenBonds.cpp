#include "HydrogenBonds.h"

extern bool THB_VERBOSE;

#ifdef PTHREADS
extern Queue<struct worker_data_s> inQueue;
extern Queue<struct worker_data_s> outQueue;
#endif

typedef std::vector<struct HydrogenBond *> HBVec;
typedef std::vector<struct HydrogenBondIterator_s> HBVecIter;

void getHydrogenBondElements( std::vector<struct thbAtom *> *atom,
                              std::vector<struct thbAtom *> *hydrogendonors,
                              std::vector<struct thbAtom *> *acceptors,
                              struct HydrogenBondMatching *match)
{
	std::vector<struct thbAtom *>::iterator it_a1;

	// The user may have specified a match more than once, when an acceptor may
	// hydrogen bond nore than once. I don't think a hydrogen can hydrogen bond
	// more than once, but I'll leave the option here.
	std::map<std::string,unsigned int>H;
	std::map<std::string,unsigned int>A;

	for(unsigned int i=0; i < match->Hydrogens.size(); ++i) {
		H[ match->Hydrogens.at(i) ]++; }

	for(unsigned int i=0; i < match->Acceptors.size(); ++i) {
		A[ match->Acceptors.at(i) ]++; }

	for( it_a1 = atom->begin(); it_a1 < atom->end(); ++it_a1)
	{
		if ( H.find((*it_a1)->ForceField) != H.end())
		{
			(*it_a1)->HydrogenBondMax = H[(*it_a1)->ForceField];
			hydrogendonors->push_back( *it_a1 );
			// hydrogendonors->insert(hydrogendonors->end(),
			//                        H[(*it_a1)->ForceField],
			//                        *it_a1);
		}
		else if( A.find((*it_a1)->ForceField) != A.end())
		{
			(*it_a1)->HydrogenBondMax = A[(*it_a1)->ForceField];
			acceptors->push_back( *it_a1 );
			// acceptors->insert(acceptors->end(),
			//                   A[(*it_a1)->ForceField],
			//                   *it_a1);
		}
	}
	DEBUG_MSG("Capacity/size of hydrogens: " << hydrogendonors->capacity() << "/" << hydrogendonors->size());
	if ( hydrogendonors->size() != hydrogendonors->capacity() ) {
		std::vector<struct thbAtom *>(*hydrogendonors).swap(*hydrogendonors);
		DEBUG_MSG("Capacity/size of hydrogens: " << hydrogendonors->capacity() << "/" << hydrogendonors->size());
	}

	DEBUG_MSG("Capacity/size of acceptors: " << acceptors->capacity() << "/" << acceptors->size());
	if ( acceptors->size() != acceptors->capacity() ) {
		std::vector<struct thbAtom *>(*acceptors).swap(*acceptors);
		DEBUG_MSG("Capacity/size of acceptors: " << acceptors->capacity() << "/" << acceptors->size());
	}

	if (THB_VERBOSE)
	{
		std::vector<std::string>::iterator it;

		VERBOSE_CMSG("Hydrogen donors : " << hydrogendonors->size() << " [ ");
		for(it=match->Hydrogens.begin(); it != match->Hydrogens.end();++it)
			VERBOSE_CMSG(*it << " ");
		VERBOSE_MSG("]");

		VERBOSE_CMSG("Acceptors       : " << acceptors->size() << " [ ");
		for(it=match->Acceptors.begin(); it != match->Acceptors.end();++it)
			VERBOSE_CMSG(*it << " ");
		VERBOSE_MSG("]");
	}
}

// Savemem version.
void HBs( HBVec *hb,
          Point cell,
          std::vector<struct thbAtom *>*hydrogens,
          std::vector<struct thbAtom *>*acceptors,
          std::vector<Point> *Coordinates,
          double TrjIdx,
          double rCutoff, double angleCutoff)
{

	std::vector<struct thbAtom *>::iterator it_h;
	std::vector<struct thbAtom *>::iterator it_a;

	Point r;

	double rCutoff2 = pow(rCutoff,2.0);
	double r2;

	// Location of the Hydrogen (H) atom, the Acceptor (A) atom, and the Donor
	// (D) covalently bonded to the Hydrogen.
	Point H,A,D;

	for( it_h = hydrogens->begin(); it_h < hydrogens->end(); ++it_h)
	{
		// location of the hydrogen of interest.
		H = Coordinates->at( (*it_h)->ID );

		for( it_a = acceptors->begin(); it_a < acceptors->end(); ++it_a)
		{
			// Make sure this acceptor is not covalently bonded to the hydrogen
			// we are looking at.  The angle check would catch this, but this
			// will skip a few calculations.
			if ( (*it_a) == (*it_h)->ConnectedAtom.at(0) )
				continue;

			// location of the acceptor atom of interest.
			A = Coordinates->at( (*it_a)->ID );

			r = H.minimumImage( A, cell );
			r2 = r.magnitudeSquared();

			if ( r2 < rCutoff2)
			{
				// Distance cutoff is good, now check the angle.
				// location of the donor atom connected to the Hydrogen.
				D = Coordinates->at( (*it_h)->ConnectedAtom.at(0)->ID );

				double angle = H.angle(A,D);
				if ( angle > angleCutoff )
				{
					struct HydrogenBond *NewHB;
					NewHB = new struct HydrogenBond;

					NewHB->length   = sqrt(r2);
					NewHB->angle    = angle;
					NewHB->hydrogen = *it_h;
					NewHB->acceptor = *it_a;
					NewHB->donor    = (*it_h)->ConnectedAtom.at(0);
					NewHB->TrajIdx  = TrjIdx;
					// NewHB->acceptorDonorDistance=D.minimumImageDistance(A,cell);

					hb->push_back(NewHB);
				}
			}
		}
	}
}

void HBs( HBVec *hb,
          Point cell,
          std::vector<struct thbAtom *>*hydrogens,
          std::vector<struct thbAtom *>*acceptors,
          double TrjIdx,
          double rCutoff, double angleCutoff)
{

	std::vector<struct thbAtom *>::iterator it_h;
	std::vector<struct thbAtom *>::iterator it_a;

	Point r;

	double rCutoff2 = pow(rCutoff,2.0);
	double r2;

	// Location of the Hydrogen (H) atom, the Acceptor (A) atom, and the Donor
	// (D) covalently bonded to the Hydrogen.
	Point H,A,D;

	for( it_h = hydrogens->begin(); it_h < hydrogens->end(); ++it_h)
	{
		// location of the hydrogen of interest.
		H = (*it_h)->p.at(TrjIdx);

		for( it_a = acceptors->begin(); it_a < acceptors->end(); ++it_a)
		{
			// Make sure this acceptor is not covalently bonded to the hydrogen
			// we are looking at.  The angle check would catch this, but this
			// will skip a few calculations.
			if ( (*it_a) == (*it_h)->ConnectedAtom.at(0) )
				continue;

			// location of the acceptor atom of interest.
			A = (*it_a)->p.at(TrjIdx);

			r = H.minimumImage( A, cell );
			r2 = r.magnitudeSquared();

			if ( r2 < rCutoff2)
			{
				// Distance cutoff is good, now check the angle.
				// location of the donor atom connected to the Hydrogen.
				D = (*it_h)->ConnectedAtom.at(0)->p.at(TrjIdx);

				double angle = H.angle(A,D);
				if ( angle > angleCutoff )
				{
					struct HydrogenBond *NewHB;
					NewHB = new struct HydrogenBond;

					NewHB->length   = sqrt(r2);
					NewHB->angle    = angle;
					NewHB->hydrogen = *it_h;
					NewHB->acceptor = *it_a;
					NewHB->donor    = (*it_h)->ConnectedAtom.at(0);
					NewHB->TrajIdx  = TrjIdx;

					NewHB->acceptorDonorDistance=D.minimumImageDistance(A,cell);

					hb->push_back(NewHB);
				}
			}
		}
	}
}

void AtomNeighbors( HBVec *hb,
                    struct PBC *Cell,
                    std::vector<struct thbAtom *>*hydrogens,
                    std::vector<struct thbAtom *>*acceptors,
                    double rCutoff, double angleCutoff )
{
	VERBOSE_MSG("Finding hydrogen bonds with:\n\n\tRc    < " << rCutoff << " Angstroms, and \n\tangle > " << angleCutoff << " degrees.\n");

	unsigned int NumFramesInTrajectory = Cell->frames;

	for( unsigned int TrjIdx=0; TrjIdx < NumFramesInTrajectory; ++TrjIdx)
	{
#ifdef PTHREADS
		struct worker_data_s wd;
		wd.jobtype = THREAD_JOB_HBS;
		wd.jobnum = TrjIdx;
		wd.num_threads = NumberOfCPUs();
		wd.cell = Cell->p.at(TrjIdx);
		wd.hydrogens = hydrogens;
		wd.acceptors = acceptors;
		wd.TrjIdx = TrjIdx;
		wd.rCutoff = rCutoff;
		wd.angleCutoff = angleCutoff;

		wd.hb = new HBVec;
		wd.hb->reserve(acceptors->size()*2);

		inQueue.push(wd);
#else
		if (  ((TrjIdx+1)%10==0) || ((TrjIdx+1)==NumFramesInTrajectory)  )
			VERBOSE_RMSG("Processing frame " << TrjIdx+1 <<"/"<< Cell->frames << ". Hydrogen-acceptor pairs found: " << hb->size() << ".");
		HBs( hb, cell, hydrogens, acceptors, TrjIdx, rCutoff, angleCutoff);
#endif
	}

#ifdef PTHREADS
	// Get the results back from the worker threads.
	for( unsigned int TrjIdx=0; TrjIdx < NumFramesInTrajectory; ++TrjIdx)
	{
		if (  ((TrjIdx+1)%10==0) || ((TrjIdx+1)==NumFramesInTrajectory)  )
			VERBOSE_RMSG("Processing frame " << TrjIdx+1 <<"/"<< Cell->frames << ". Hydrogen-acceptor pairs found: " << hb->size() << ".");

		struct worker_data_s wd = outQueue.pop();
		hb->reserve( hb->size() + wd.hb->size() );
		hb->insert(hb->end(),
		           wd.hb->begin(),
		           wd.hb->end() );

		wd.hb->clear();
		delete wd.hb;
	}
#endif // PTHREADS

	VERBOSE_MSG("Processing frame " << Cell->frames <<"/"<< Cell->frames << ". Hydrogen-acceptor pairs found: " << hb->size() << ".");
}
