/**
 * \file
 * \brief Simple math operations on vector of doubles
 */
#ifndef _SimpleMath_h
#define _SimpleMath_h

#include<vector>

/**
 * Simple Math namespace
 */
namespace SimpleMath {
	/**
	 * Calculate distance between two points given their coordinates in three
	 * dimensional Cartesian space.
	 *
	 * \param[in] x1   x-coordinate of point 1
	 * \param[in] y1   y-coordinate of point 1
	 * \param[in] z1   z-coordinate of point 1
	 * \param[in] x2   x-coordinate of point 2
	 * \param[in] y2   y-coordinate of point 2
	 * \param[in] z2   z-coordinate of point 2
	 *
	 * \return Distance between two points
	 */
	double distance( double x1, double y1, double z1,
	                 double x2, double y2, double z2);
	/**
	 * Sum all values
	 *
	 * \param v  Vector of doubles.
	 *
	 * \return Sum of all doubles.
	 */
	double Sum( std::vector<double> v);

	/**
	 * Sum of the squared differences
	 *
	 * \f[
	 * \sum\limits_{i=0}^N\left(v_i-\mathrm{avg}\right)^2
	 * \f]
	 *
	 * \param[in] v    Vector of doubles
	 * \param[in] avg  Average of all values in \p v
	 *
	 * \return Sum of squared differences
	 */
	double SumSquaredDifferences( std::vector<double> v,double avg);
	/**
	 * Sum of the squared differences
	 *
	 * \f[
	 * \sum\limits_{i=0}^N\left(v_i-\left<v\right>\right)^2
	 * \f]
	 *
	 * \param[in] v  Vector of doubles
	 *
	 * \return Sum of squared differences
	 */
	double SumSquaredDifferences( std::vector<double> v);

	/**
	 * Average
	 *
	 * \param[in] v  Vector of doubles
	 *
	 * \return Average of all values in \p v
	 */
	double average( std::vector<double> v);

	/**
	 * Corrected Sample standard deviation
	 *
	 * \f[
	 * \sqrt\frac{\sum\limits_{i=0}^N\left(v_i-\mathrm{avg}\right)^2}{N-2}
	 * \f]
	 *
	 * \param[in] v    Vector of doubles
	 * \param[in] avg  Average of all values in \p v
	 *
	 * \return Corrected Sample standard deviation
	 */

	double stddev( std::vector<double> v, double avg);
	/**
	 * Corrected Sample standard deviation
	 *
	 * \f[
	 * \sqrt\frac{\sum\limits_{i=0}^N\left(v_i-\left<v\right>\right)^2}{N-2}
	 * \f]
	 *
	 * \param[in] v   Vector of doubles
	 *
	 * \return Corrected Sample standard deviation
	 */
	double stddev( std::vector<double> v);

	/**
	 * Corrected Sample standard deviation
	 *
	 * \f[
	 * \sqrt\frac{\mathrm{sumSq}}{N-2}
	 * \f]
	 *
	 * \param[in] sumSq   Sum of the squared differences of a series of values
	 * \param[in] N       Number of values
	 *
	 * \return Corrected Sample standard deviation
	 */
	double stddev( double sumSq, unsigned int N);
}
#endif
