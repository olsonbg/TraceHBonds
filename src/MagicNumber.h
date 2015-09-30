/**
 *
 * \file
 * \date   21 April 2015
 * \brief  Determine magic number of file
 *
 **/
#ifndef _MagicNumber_h
#define _MagicNumber_h

#include<iostream>
#include<fstream>
#include<stdint.h>

const int MAGICNUMBER_UNKNOWN =  0; /**< Unknown file format */
const int MAGICNUMBER_GZIP    =  1; /**< GZip file           */
const int MAGICNUMBER_BZIP2   =  2; /**< Bzip2 file          */
const int MAGICNUMBER_LZMA    =  3; /**< Lzma file (xz)      */

/** Magic numbers. */
const int MagicNumberKeys[][7] =
{
	{MAGICNUMBER_GZIP    , 0x1f, 0x8b, 0x08,   -1,   -1    -1},
	{MAGICNUMBER_BZIP2   ,  'B',  'Z',  'h',   -1,   -1,   -1},
	{MAGICNUMBER_LZMA    , 0xfd,  '7',  'z',  'X',  'Z', 0x00},
	{MAGICNUMBER_UNKNOWN ,   -1,   -1,   -1,   -1,   -1, 0x00}
};

/**
 *
 * Currently supports, GZIP, BZIP2, and LZMA.
 *
 * \param filename   Name of files to determine magic number of
 *
 * \return Integer corresponding to type of file as determined from the magic
 * number.
 *
 **/
int getMagicNumber(const char *filename);
#endif
