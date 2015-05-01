#ifndef _Point_h
#define _Point_h

#include <math.h>

const double PI = 3.14159265358979323846;
class Point {
	private:
		double u1;
		double u2;
		double u3;
		double Round (double r, double f=1.0);

		// Round each element.
		Point Round (Point v, double f=1.0);

	public:
		// Constructor: Can be called with zero, one, two, or three
		// arguments.
		Point(double x=0.0, double y = 0.0, double z = 0.0);

		double x();
		double y();
		double z();

		double distance( Point v );
		double distanceSquared( Point v );
		double magnitude();
		double magnitudeSquared();

		// Minimum Image vector. c is the size of the periodic cell.
		Point minimumImage( Point v, Point c );

		// Minimum Image distance. c is the size of the periodic cell.
		double minimumImageDistance( Point v, Point c );

		// dot product
		double dot(Point v);

		// cross product
		Point cross(Point v);

		// Angle between a-this-b
		double angle(Point a, Point b);
		// Addition operator
		Point operator+(Point v);

		// Subtraction operator
		Point operator-(Point v);

		// Divide operator (divide each component)
		Point operator/(Point v);

		// Multiply operator (multiply each component)
		Point operator*(Point v);

		// Equivalent operator
		bool operator==(Point v);
};
#endif // _Point_h
