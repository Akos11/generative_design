#include "CircleSolid.h"


void CircleSolid::setCenter(Vector _center) {
	center = _center;
}
void CircleSolid::setAxis(Vector _axis) {
	axis = _axis;
}
void CircleSolid::setR(double _r) {
	r = _r;
}
bool CircleSolid::isInside(Vector point)
{
	if ((point - center).length() < r + accuracy)
		return true;
	else
		return false;
}
