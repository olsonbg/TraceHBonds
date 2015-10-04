#ifndef _sizehist_h
#define _sizehist_h
void sizehist(unsigned int NumFramesInTrajectory,
              unsigned int everyNth,
              char *ofPrefix, char *ofSuffix,
              HBVec *hb,
              struct PBC *Cell,
              int NumBins,
              bool povray,
              std::vector<ListOfHBonds *>*HBStrings);

#endif // _sizehist_h
