#include <math.h>
#include <vector>
#include "SimpleMath.h"

namespace SimpleMath
{
	double distance( double x1, double y1, double z1,
	                 double x2, double y2, double z2)
	{
		return(sqrt( pow(x2-x1,2) + pow(y2-y1,2) + pow(z2-z1,2)));
	}

	double Sum( std::vector<double> v)
	{
		std::vector<double>::iterator it;
		double sum = 0.0;

		for(it = v.begin(); it < v.end(); ++it)
			sum += *it;

		return(sum);
	}

	double SumSquaredDifferences( std::vector<double> v, double avg)
	{
		std::vector<double>::iterator it;
		double sumSq = 0.0;

		for(it = v.begin(); it < v.end(); ++it)
			sumSq += pow(*it - avg,2.0);

		return( sumSq );
	}

	double SumSquaredDifferences( std::vector<double> v )
	{
		double avg = average(v);

		return(SumSquaredDifferences(v, avg));
	}

	double average( std::vector<double> v)
	{
		std::vector<double>::iterator it;
		double sum = 0.0;

		if ( v.size() == 0) return(0.0);

		for(it = v.begin(); it < v.end(); ++it)
			sum += *it;

		return(sum/v.size());
	}

	double stddev( double sum, unsigned int N)
	{
		unsigned int E = 1.0; // 1.0 == Corrected Sample standard deviation.
		                      // 0.0 == Uncorrected sample standard deviation
		                      //        (Standard deviation of the sample)

		if ( N < (E+1) ) return(0.0);

		return( sqrt(sum/double(N-E)) );
	}

	double stddev( std::vector<double> v, double avg)
	{
		std::vector<double>::iterator it;
		double sum = 0.0;

		for(it = v.begin(); it < v.end(); ++it)
			sum += pow(*it - avg,2.0);

		return( stddev(sum,v.size()) );
	}

	double stddev( std::vector<double> v)
	{
		double avg = average(v);

		return( stddev(v, avg) );
	}
}
