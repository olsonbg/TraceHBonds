#include "MagicNumber.h"

int getMagicNumber( const char *filename )
{
	uint8_t magic[6];
	FILE *fp ;

	fp = fopen(filename,"r");

	if( fp == NULL )
	{
		std::perror(filename);
		return(-1);
	}

	if ( fread(magic, 1, 6, fp) != 6 )
	{
		std::cerr << "Error reading magic number." << "\n";
		fclose(fp);
		return(MAGICNUMBER_UNKNOWN);
	}
	fclose(fp);

	unsigned int mn = 0;
	while( MagicNumberKeys[mn][0] != MAGICNUMBER_UNKNOWN )
	{
		bool FoundMatch=true;
		for (unsigned int i = 0; i < 6; i++)
		{
			if( MagicNumberKeys[mn][i+1] == -1 )
				break;
			else if( magic[i] != MagicNumberKeys[mn][i+1] )
			{
				FoundMatch=false;
				break;
			}
		}
		if ( FoundMatch )
			return(MagicNumberKeys[mn][0]);
		++mn;
	}

	return(MAGICNUMBER_UNKNOWN);
}
