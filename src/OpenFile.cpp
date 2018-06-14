#include "MagicNumber.h"
#include "OpenFile.h"


bool openfile( const char *filename,
               boost::iostreams::filtering_stream<boost::iostreams::input> *in,
               std::ifstream *ifp)
{
	int magicNum;

	magicNum = getMagicNumber(filename);

	if ( magicNum == -1 )
	{
		std::cout << "Error opening " << filename << "\n";
		return(false);
	}


#ifdef USE_ZLIB
	if ( magicNum == MAGICNUMBER_GZIP )
	{
		in->push(boost::iostreams::gzip_decompressor());
		ifp->open(filename,std::ios::in|std::ios::binary);
	}
#endif
#ifdef USE_BZIP2
	if ( magicNum == MAGICNUMBER_BZIP2 )
	{
		in->push(boost::iostreams::bzip2_decompressor());
		ifp->open(filename,std::ios::in|std::ios::binary);
	}
#endif
#ifdef USE_LZMA
	if ( magicNum == MAGICNUMBER_LZMA )
	{
		in->push(lzma_input_filter());
		ifp->open(filename,std::ios::in|std::ios::binary);
	}
#endif

	if ( magicNum == MAGICNUMBER_UNKNOWN )
		ifp->open(filename,std::ios::in);

	if ( !ifp->is_open() )
	{
		std::cerr << "Error: can not open this type of file." << "\n";
		return(false);
	}

	in->push(*ifp);
	return(true);
}
