#ifndef _Trace_h
#define _Trace_h

extern bool THB_VERBOSE;

bool TraceThread( std::vector<ListOfHBonds *> *HBStrings,
                  struct HydrogenBondIterator_s *HBit);

void Trace( std::vector<ListOfHBonds *> *HBStrings,
            std::vector<struct HydrogenBondIterator_s> *TrjIdx_iter );
#endif // _Trace_h
