#include <math.h>
#include "Point.h"

double Point::Round (double r, double f) {
	return (r > 0.0) ? floor(r*f + 0.5)/f : ceil(r*f - 0.5)/f; }

// Round each element.
Point Point::Round (Point v, double f) {
	return Point( Round( v.x(),f ),
	              Round( v.y(),f ),
	              Round( v.z(),f ) ); }

// Constructor: Can be called with zero, one, two, or three
// arguments.
Point::Point(double xx, double yy, double zz) {
	u1 = xx;
	u2 = yy;
	u3 = zz;
}

double Point::x() {return u1; }
double Point::y() {return u2; }
double Point::z() {return u3; }

double Point::distanceSquared( Point v ) {
	double xd = u1 - v.x();
	double yd = u2 - v.y();
	double zd = u3 - v.z();

	return( xd*xd + yd*yd + zd*zd );
}

double Point::distance( Point v ) {
	return sqrt( distanceSquared(v) ) ; }

double Point::magnitude() {
	Point dummy; // dummy is (0,0,0)
	return( distance(dummy) );
}

double Point::magnitudeSquared() {
	Point dummy; // dummy is (0,0,0)
	return( distanceSquared(dummy) );
}
// Minimum Image vector. c is the size of the periodic cell.
Point Point::minimumImage( Point v, Point c ) {
	Point d = v - *this;

	d = d - Round(d/c)*c;

	return d;
}

// Minimum Image distance. c is the size of the periodic cell.
double Point::minimumImageDistance( Point v, Point c) {
	return distance( minimumImage(v,c) ); }

// dot product
double Point::dot(Point v) {
	return ( u1 * v.x() + u2 * v.y() + u3 * v.z() ); }

// cross product
Point Point::cross(Point v) {
	return Point( u2*v.z() - u3*v.y(),
				  u3*v.x() - u1*v.z(),
				  u1*v.y() - u2*v.x() ); }

// Angle between a-this-b.
double Point::angle(Point a, Point b)
{
	Point u = a - *this;
	Point v = b - *this;
	// Point u( a.x() - u1, a.y() - u2, a.z() - u3);
	// Point v( b.x() - u1, b.y() - u2, b.z() - u3);

	return acos( u.dot(v)/u.magnitude()/v.magnitude() )*180.0/M_PI;
}

// Addition operator
Point Point::operator+(Point v) {
	return Point(u1 + v.x(), u2 + v.y(), u3 + v.z() ); }

// Subtraction operator
Point Point::operator-(Point v) {
	return Point(u1 - v.x(), u2 - v.y(), u3 - v.z() ); }

// Divide operator
Point Point::operator/(double v) {
	return Point(u1/v, u2/v, u3/v ); }

// Multiply operator
Point Point::operator*(double v) {
	return Point( u1*v, u2*v, u3*v); }

// Divide operator (divide each component)
Point Point::operator/(Point v) {
	return Point(u1/v.x(), u2/v.y(), u3/v.z() ); }

// Multiply operator (multiply each component)
Point Point::operator*(Point v) {
	return Point( u1*v.x(), u2*v.y(), u3*v.z()); }

// Equivalent operator
bool Point::operator==(Point v) {
	if ( (u1==v.x()) &&
	     (u2==v.y()) &&
	     (u3==v.z()) )
		return true;

	return false;
}
