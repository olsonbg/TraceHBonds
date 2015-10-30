#include "flags.h"
#include <type_traits>

unsigned int operator& (unsigned int l, Flags r) {
	typedef std::underlying_type<Flags>::type unint;
	return( l & unint(r) );
}

unsigned int operator| (unsigned int l, Flags r) {
	typedef std::underlying_type<Flags>::type unint;
	return( l | unint(r) );
}

Flags operator& (Flags l, Flags r) {
	typedef std::underlying_type<Flags>::type unint;
	return Flags( unint(l) & unint(r) );
}

Flags operator| (Flags l, Flags r) {
	typedef std::underlying_type<Flags>::type unint;
	return Flags( unint(l) | unint(r) );
}
