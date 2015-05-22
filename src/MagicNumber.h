#ifndef _MagicNumber_h
#define _MagicNumber_h

#include<iostream>
#include<fstream>
#include<stdint.h>

const int MAGICNUMBER_UNKNOWN =  0;
const int MAGICNUMBER_GZIP    =  1;
const int MAGICNUMBER_BZIP2   =  2;
const int MAGICNUMBER_LZMA    =  3;

const int MagicNumberKeys[][7] =
{
	{MAGICNUMBER_GZIP    , 0x1f, 0x8b, 0x08,   -1,   -1    -1},
	{MAGICNUMBER_BZIP2   ,  'B',  'Z',  'h',   -1,   -1,   -1},
	{MAGICNUMBER_LZMA    , 0xfd,  '7',  'z',  'X',  'Z', 0x00},
	{MAGICNUMBER_UNKNOWN ,   -1,   -1,   -1,   -1,   -1, 0x00}
};

int getMagicNumber(const char *filename);
#endif
