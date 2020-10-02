#pragma once
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <QDebug>

class ConstraintSolid {
private:
	bool keepIn;
	bool keepOut;
protected:
	double accuracy = 1;
	//int idx;
public:
	using Vector = OpenMesh::VectorT<double, 3>;
	void setKeepIn(bool b);
	void setKeepOut(bool b);
	/*void setIdx(int idx);
	int getIdx();*/
	virtual bool isInside(OpenMesh::Vec3d point);
	virtual ~ConstraintSolid() {}
};