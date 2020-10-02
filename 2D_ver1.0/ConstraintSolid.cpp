#include "ConstraintSolid.h"

void ConstraintSolid::setKeepIn(bool b) {
	keepIn = b;
}

void ConstraintSolid::setKeepOut(bool b) {
	keepOut = b;
}
bool ConstraintSolid::isInside(Vector point) {
	return false;
}