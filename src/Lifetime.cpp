#include "Lifetime.h"

std::vector< std::vector<bool> >
Lifetime( HBVecIter *TrjIdx_iter )
{
	// Look to see if this hydrogen bond exists in other frames.

	HBVec::iterator iter_hb;
	HBVec::iterator iter_begin;
	HBVec::iterator iter_end;

	unsigned int NumFrames = TrjIdx_iter->size();
	iter_begin = TrjIdx_iter->at( 0 ).begin;
	iter_end   = TrjIdx_iter->at( 0 ).end;

	HBVec::iterator iter_hbmain = iter_begin;

	unsigned int NumHBsInFrameZero=0;

	for(iter_hb = iter_begin ; iter_hb < iter_end; ++iter_hb ) {
		NumHBsInFrameZero++; }
	std::cout << " Number of hydrogen bonds in initial frame: " << NumHBsInFrameZero << "\n";

	std::vector< std::vector<bool> >b(NumHBsInFrameZero, std::vector<bool>(NumFrames,false));

	if (NumHBsInFrameZero == 0)
		return b;

	for( unsigned int h=0; h < NumHBsInFrameZero; ++h)
	{
		b.at(h).at(0) = true;
		for( unsigned int f=1; f < NumFrames; ++f )
		{
			// The range of hydrogen bonds in frame f.
			iter_begin = TrjIdx_iter->at( f ).begin;
			iter_end   = TrjIdx_iter->at( f ).end;

			bool found = false;
			for(iter_hb = iter_begin ; iter_hb < iter_end; ++iter_hb )
			{
				if ( ((*iter_hb)->hydrogen == (*(iter_hbmain+h))->hydrogen) &&
				     ((*iter_hb)->acceptor == (*(iter_hbmain+h))->acceptor) )
				{
					// This matched the hydrogen bond in the first frame.
					found = true;
					break;
				}
			}
			if (found == true) {
				b.at(h).at(f) = true; }
		}
	}
	// for( std::vector< std::vector<bool> >::iterator vbit=b.begin(); vbit < b.begin()+1; ++vbit)
	// {
	//     for( std::vector<bool>::iterator bit=vbit->begin(); bit < vbit->end(); ++bit)
	//     {
	//         std::cout << (*bit==true?"|":" ");
	//     }
	// std::cout << "\n";
	// }

	return b;
}
