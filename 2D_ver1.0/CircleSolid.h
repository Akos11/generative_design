#pragma once
#include "ConstraintSolid.h"


class CircleSolid : public ConstraintSolid {
private:
	Vector center;
	Vector axis;
	double r;

public:

	void setCenter(Vector _center);
	void setAxis(Vector _axis);
	void setR(double _r);
	virtual bool isInside(Vector point) override;
	~CircleSolid() {
	}
};