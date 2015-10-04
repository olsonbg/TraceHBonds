#ifndef _neighborhist_h
#define _neighborhist_h
void neighborhist(unsigned int NumFramesInTrajectory,
                  unsigned int everyNth,
                  char *ofPrefix, char *ofSuffix,
                  struct PBC *Cell,
                  std::vector<ListOfHBonds *>*HBStrings);

#endif // _neighborhist_h
