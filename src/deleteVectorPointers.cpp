#include <algorithm>
#include "deleteVectorPointers.h"
#include "ListOfHBonds.h"

// Used with remove_if()
template <typename T> bool deleteVectorPointers( T* v ) {
	delete v;
	return true;
}

template<class T> void DeleteVectorPointers( std::vector<T *> v )
{
	remove_if(v.begin(),v.end(),deleteVectorPointers<T>);
}

// Explicit template instantiation

// Used in: sizehist.cpp
template void DeleteVectorPointers<ListOfHBonds>(std::vector<ListOfHBonds *>);
template bool deleteVectorPointers<ListOfHBonds>(ListOfHBonds *);

// Used in: TraceHBonds.cpp
template void DeleteVectorPointers<thbAtom>(std::vector<thbAtom *>);
template bool deleteVectorPointers<thbAtom>(thbAtom *);

// Used in: TraceHBonds.cpp
template void DeleteVectorPointers<thbMolecule>(std::vector<thbMolecule *>);
template bool deleteVectorPointers<thbMolecule>(thbMolecule *);

// Used in: TraceHBonds.cpp
template void DeleteVectorPointers<thbBond>(std::vector<thbBond *>);
template bool deleteVectorPointers<thbBond>(thbBond *);

// Used in: TraceHBonds.cpp
template void DeleteVectorPointers<HydrogenBond>(std::vector<HydrogenBond *>);
template bool deleteVectorPointers<HydrogenBond>(HydrogenBond *);
