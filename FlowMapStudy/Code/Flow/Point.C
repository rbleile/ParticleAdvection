#include "Point.h"	

Point::Point()
{
	x = 0.0;
	y = 0.0;
	z = 0.0;
	t = 0.0;
}

Point::Point( double ix, double iy, double iz )
{
	x = ix;
	y = iy;
	z = iz;
	t = 0.0;
}

Point::Point( double ix, double iy, double iz, double it )
{
	x = ix;
	y = iy;
	z = iz;
	t = it;
}

void Point::setPoint( double ix, double iy, double iz )
{
	x = ix;
	y = iy;
	z = iz;
}

void Point::setPoint( double ix, double iy, double iz, double it )
{
	x = ix;
	y = iy;
	z = iz;
	t = it;
}

