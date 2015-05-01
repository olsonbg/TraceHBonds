#ifndef _SimpleMath_h
#define _SimpleMath_h

#include<vector>

namespace SimpleMath {
	double distance( double x1, double y1, double z1,
	                 double x2, double y2, double z2);
	double Sum( std::vector<double> v);
	double SumSquaredDifferences( std::vector<double> v,double avg);
	double SumSquaredDifferences( std::vector<double> v);
	double average( std::vector<double> v);
	double stddev( std::vector<double> v, double avg);
	double stddev( std::vector<double> v);
	double stddev( double avg, unsigned int N);
}
#endif
