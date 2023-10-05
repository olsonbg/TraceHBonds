/**
 * \file
 * \date  30 April 2015
 * \brief Class for 3 dimensional coordinates and simple calculation
 **/
#ifndef _Point_h
#define _Point_h

#define _USE_MATH_DEFINES
#include <math.h>

//const double PI = 3.14159265358979323846;

/**
 *
 * Contains 3 dimensional coordinates, and performs simple calculations on them.
 *
 * \anchor uDef
 *  u is defined here/
 **/
class Point {
	private:
		double u1; /**< x-coordinate */
		double u2; /**< y-coordinate */
		double u3; /**< z-coordinate */
		double Round (double r, double f=1.0);

		// Round each element.
		Point Round (Point v, double f=1.0);

	public:
		/**
		 *
		 * Constructor: Can be called with zero, one, two, or three
		 * arguments.
		 *
		 * \param x  X-coordinate of point
		 * \param y  Y-coordinate of point
		 * \param z  Z-coordinate of point
		 *
		 **/
		Point(double x=0.0, double y = 0.0, double z = 0.0);

		double x(); /**< \return x-coordinate of point */
		double y(); /**< \return y-coordinate of point */
		double z(); /**< \return z-coordinate of point */

		/**
		* Magnitude of vector pointing from \p v to \p u, where \p u has
		* coordinates given by x(), y(), and z().
		*  \f[
		*  \sqrt{(u_x - v_x)^2 + (u_y - v_y)^2 + (u_z - v_z)^2}
		*  \f]
		*
		* \return Scalar distance between \p u and \p v
		*/
		double distance( Point v );
		/**
		* Squared magnitude of vector pointing from \p v to \p u, where \p u
		* has coordinates given by x(), y(), and z().
		*
		*  \f[
		*  (u_x - v_x)^2 + (u_y - v_y)^2 + (u_z - v_z)^2
		*  \f]
		*
		* \return Squared distance between \p u and \p v
		*/
		double distanceSquared( Point v );
		/**
		* Magnitude of vector from origin to u, where \p u has coordinates
		* given by x(), y(), and z().
		*
		* \see distance() with \f$v=(0,0,0)\f$
		*
		* \return Magnitude of vector \p u
		**/
		double magnitude();
		/**
		 * Squared Magnitude of vector from origin to u, where \p u has
		 * coordinates given by x(), y(), and z().
		 *
		 * \see distanceSquared() with \f$v=(0,0,0)\f$
		 *
		 * \return Squared magnitude of vector \p u
		 */
		double magnitudeSquared();

		/**
		 * Minimum Image vector pointing to from \p P to \p v.
		 *
		 * \param v   Point of interest
		 * \param c   Dimensions of periodic cell
		 *
		 * \return Minimum image vector (as \p Point) to \p v
		 */
		Point minimumImage( Point v, Point c );

		/**
		 * Minimum Image distance.
		 *
		 * \param v   Point of interest
		 * \param c   Dimensions of periodic cell
		 *
		 * \return Minimum image distance
		 */
		double minimumImageDistance( Point v, Point c );

		/**
		 * dot product: \f$u\cdot v\f$, where \p u has coordinates given by
		 * x(), y(), and z().
		 *
		 * \param[in] v  Point of interest
		 *
		 * \return Dot product with \p v
		 */
		double dot(Point v);

		/**
		 * cross product: \f$u\times v\f$, where \p u has coordinates given by
		 * x(), y(), and z().
		 *
		 * \param[in] v  Point of interest
		 *
		 * \return Cross product
		 */
		Point cross(Point v);

		/**
		 * Angle formed between two vectors, the first pointing from \p u to \p a, and the second pointing from \p u to \p b, where \p u has coordinates given by x(), y(), and z().
		 *
		 * The calculated angle will be between 0 and \f$\pi\f$, inclusive.
		 *
		 * \param[in] a  Point
		 * \param[in] b  Point
		 *
		 * \return Angle, where \f$0\le \mathrm{Angle} \le \pi\f$
		 */
		double angle(Point a, Point b);

		/**
		 * Addition operator
		 *
		 * Vector addition of  \p v and \p u, where \p u has coordinates given by x(), y(), and z().
		 *
		 * \param[in] v Point
		 *
		 * \return vector sum of \p u and \p v
		 **/
		Point operator+(Point v);

		/**
		 * Subtraction operator
		 *
		 * Vector subtraction of \p v from \p u, where \p u has coordinates given by x(), y(), and z().
		 *
		 * \param[in] v Point
		 *
		 * \return vector subtraction of \p v from \p u
		 **/
		Point operator-(Point v);

		/**
		 * Division operator
		 *
		 * Each component of \p u is divided by \p v, where \p u has coordinates
		 * given by x(), y(), and z().
		 *
		 * \param[in] v Double
		 *
		 * \return Point
		 **/
		Point operator/(double v);

		/**
		 * Multiply operator
		 *
		 * Each component of \p u is multiplied by \p v, where \p u has
		 * coordinates given by x(), y(), and z().
		 *
		 * \param[in] v Double
		 *
		 * \return Double
		 **/
		Point operator*(double v);

		/**
		 * Division operator
		 *
		 * Each component of \p u is divided by the respective component of \p
		 * v, where \p u has coordinates given by x(), y(), and z().
		 *
		 * \param[in] v Point
		 *
		 * \return Point
		 **/
		Point operator/(Point v);

		/**
		 * Multiply operator
		 *
		 * Each component of \p u is multiplied by the respective component of \p
		 * v, where \p u has coordinates given by x(), y(), and z().
		 *
		 * \param[in] v Point
		 *
		 * \return Point
		 **/
		Point operator*(Point v);

		/**
		 * Equivalent operator
		 *
		 * Compared each component of \p u with the respective component of \p
		 * v, where \p u has coordinates given by x(), y(), and z(). If all
		 * components are equivalent, the return \c TRUE.
		 *
		 * \param[in] v Point
		 *
		 * \return \c TRUE of all components of \p u and \p v are equivalent, \c FALSE otherwise
		 **/
		bool operator==(Point v);
};
#endif // _Point_h
