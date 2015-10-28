#ifndef _sizehist_h
#define _sizehist_h
void sizehist(unsigned int NumFramesInTrajectory,
              char *ofPrefix, char *ofSuffix,
              HBVec *hb,
              struct PBC *Cell,
              int NumBins,
              unsigned char flags,
              std::vector<ListOfHBonds *>*HBStrings);

#endif // _sizehist_h
