#include <fstream>
#include <iostream>
#include <vector>
#include <string.h>
#include "SimpleMath.h"
#include "NeighborPrint.h"


void Print_AllFrames(std::ostream *out,
                     std::vector<struct Histograms_s> *frame)
{
	for(unsigned int f=0; f< frame->size();++f)
	{
		if (f!=0) *out << "\t";
		*out << "\t\tFrame " << frame->at(f).TrjIdx + 1;
	}
	*out << "\n";

	*out << "Atoms in chain\tNth n.n.";
	for(unsigned int f=0; f< frame->size();++f)
		*out << "\tCount\tAverage\tStdDev";
	*out << "\n";

	unsigned int MaxHydrogenBondCount;
	unsigned int MaxNeighborElements=0;
	// Find the MaxHydrogenBondCount across all frames.
	for( unsigned int f=0; f< frame->size(); ++f )
	{
		if ( MaxNeighborElements < frame->at(f).NearestNeighbors.size() )
			MaxNeighborElements = frame->at(f).NearestNeighbors.size();
	}
	MaxHydrogenBondCount = (-1 + sqrt(1+8*MaxNeighborElements))/2;
	unsigned int offset=0;
	for( unsigned int i = 0; i < MaxHydrogenBondCount; i++)
	{
		offset = i*(i+1)/2;
		*out <<  2*i+3 << "\t";
		for( unsigned int j=0; j < i+1; j++ )
		{
			if ( j!=0 ) *out << " \t";
			*out << j+1;
			for( unsigned int f=0; f < frame->size(); ++f)
			{
				double avg = 0.0;
				double err = 0.0;
				unsigned int count = 0.0;

				if ( ((offset+j) < frame->at(f).NearestNeighbors.size()) &&
					 (frame->at(f).NearestNeighbors.at(offset+j).at(0) != -1) )
				{
					avg = SimpleMath::average(frame->at(f).NearestNeighbors.at(offset+j));
					err = SimpleMath::stddev(frame->at(f).NearestNeighbors.at(offset+j),avg);
					count = frame->at(f).NearestNeighbors.at(offset+j).size();
				}
				*out << "\t" << count << "\t" << avg << "\t" << err;
			}
			*out << "\n";
		}
		*out << "\n";
	}
}

void Print_CombineFrames(std::ostream *out,
                         std::vector<struct Histograms_s> *frame)
{
	// Combine all frames for an overall average.
	*out << "Combining all columns:" << "\n";
	*out << "Atoms in chain\tNth n.n.";
	*out << "\tCount\tAverage\tStdDev" << "\n";

	unsigned int MaxHydrogenBondCount;
	unsigned int MaxNeighborElements=0;
	// Find the MaxHydrogenBondCount across all frames.
	for( unsigned int f=0; f< frame->size(); ++f )
	{
		if ( MaxNeighborElements < frame->at(f).NearestNeighbors.size() )
			MaxNeighborElements = frame->at(f).NearestNeighbors.size();
	}
	MaxHydrogenBondCount = (-1 + sqrt(1+8*MaxNeighborElements))/2;
	unsigned int offset=0;

	for( unsigned int i = 0; i < MaxHydrogenBondCount; i++)
	{
		offset = i*(i+1)/2;
		*out <<  2*i+3 << "\t";
		for( unsigned int j=0; j < i+1; j++ )
		{
			if ( j!=0 ) *out << " \t";
			*out << j+1;

			double sum=0.0;
			unsigned int count = 0;
			double SumSqDiff = 0.0;
			double avg = 0.0;
			double err = 0.0;

			for( unsigned int f=0; f < frame->size(); ++f)
			{
				if ( ((offset+j) < frame->at(f).NearestNeighbors.size()) &&
					 (frame->at(f).NearestNeighbors.at(offset+j).at(0) != -1) )
				{
					count += frame->at(f).NearestNeighbors.at(offset+j).size();
					sum += SimpleMath::Sum(frame->at(f).NearestNeighbors.at(offset+j));
				}
			}
			if (count != 0) avg = sum/count;

			for( unsigned int f=0; f < frame->size(); ++f)
			{
				if ( ((offset+j) < frame->at(f).NearestNeighbors.size()) &&
					 (frame->at(f).NearestNeighbors.at(offset+j).at(0) != -1) )
					SumSqDiff += SimpleMath::SumSquaredDifferences(frame->at(f).NearestNeighbors.at(offset+j),avg);
			}
			err = SimpleMath::stddev(SumSqDiff,count);

			*out << "\t" << count << "\t" << avg << "\t" << err;
			*out << "\n";
		}
		*out << "\n";
	}
	return;
}

void Print_CombineNeighbors(std::ostream *out,
                            std::vector<struct Histograms_s> *frame)
{
	// Combine all frames and all Nth nearest ChainLengths.
	*out << "Combining all columns and Chains:" << "\n";
	*out << "Nth n.n.";
	*out << "\tCount\tAverage\tStdDev" << "\n";

	unsigned int MaxHydrogenBondCount;
	unsigned int MaxNeighborElements=0;
	// Find the MaxHydrogenBondCount across all frames.
	for( unsigned int f=0; f< frame->size(); ++f )
	{
		if ( MaxNeighborElements < frame->at(f).NearestNeighbors.size() )
			MaxNeighborElements = frame->at(f).NearestNeighbors.size();
	}
	MaxHydrogenBondCount = (-1 + sqrt(1+8*MaxNeighborElements))/2;
	unsigned int offset=0;

	for( unsigned int j = 0; j < MaxHydrogenBondCount; j++)
	{
		offset = j*(j+1)/2;
		*out << j+1;

		double sum=0.0;
		unsigned int count = 0;
		double SumSqDiff = 0.0;
		double avg = 0.0;
		double err = 0.0;

		for( unsigned int f=0; f < frame->size(); ++f)
			for( unsigned int i = j; i < MaxHydrogenBondCount; i++)
			{
				offset = i*(i+1)/2;
				if ( ((offset+j) < frame->at(f).NearestNeighbors.size()) &&
					 (frame->at(f).NearestNeighbors.at(offset+j).at(0) != -1) )
				{
					count += frame->at(f).NearestNeighbors.at(offset+j).size();
					sum += SimpleMath::Sum(frame->at(f).NearestNeighbors.at(offset+j));
				}
			}
		if (count != 0) avg = sum/count;

		for( unsigned int f=0; f < frame->size(); ++f)
			for( unsigned int i = j; i < MaxHydrogenBondCount; i++)
			{
				offset = i*(i+1)/2;
				if ( ((offset+j) < frame->at(f).NearestNeighbors.size()) &&
					 (frame->at(f).NearestNeighbors.at(offset+j).at(0) != -1) )
					SumSqDiff += SimpleMath::SumSquaredDifferences(frame->at(f).NearestNeighbors.at(offset+j),avg);
			}
		err = SimpleMath::stddev(SumSqDiff,count);

		*out << "\t" << count << "\t" << avg << "\t" << err;
		*out << "\n";
	}
	*out << "\n";
}

// vim:tw=76:ts=4:sw=4
