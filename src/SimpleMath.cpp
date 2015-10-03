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

	double Sum( vd v)
	{
		vd::iterator it;
		double sum = 0.0;

		for(it = v.begin(); it < v.end(); ++it)
			sum += *it;

		return(sum);
	}

	double SumSquaredDifferences( vd v, double avg)
	{
		vd::iterator it;
		double sumSq = 0.0;

		for(it = v.begin(); it < v.end(); ++it)
			sumSq += pow(*it - avg,2.0);

		return( sumSq );
	}

	double SumSquaredDifferences( vd v )
	{
		double avg = average(v);

		return(SumSquaredDifferences(v, avg));
	}

	double average( vd v)
	{
		vd::iterator it;
		double sum = 0.0;

		if ( v.size() == 0) return(0.0);

		for(it = v.begin(); it < v.end(); ++it)
			sum += *it;

		return(sum/v.size());
	}

	double stddev( double sumSq, unsigned int N)
	{
		unsigned int E = 1.0; // 1.0 == Corrected Sample standard deviation.
		                      // 0.0 == Uncorrected sample standard deviation
		                      //        (Standard deviation of the sample)
		// Return 0.0 if not enough data point to calculate standard deviation
		if ( N <= (E+1) ) return(0.0);

		return( sqrt(sumSq/double(N-E)) );
	}

	double stddev( vd v, double avg)
	{
		vd::iterator it;
		double sumSq = SumSquaredDifferences(v, avg);

		return( stddev(sumSq,v.size()) );
	}

	double stddev( vd v)
	{
		double avg = average(v);

		return( stddev(v, avg) );
	}
}
