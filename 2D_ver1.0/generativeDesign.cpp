#include "MyViewer.h"
#include <OpenMesh/Tools/Smoother/JacobiLaplaceSmootherT.hh>
#include <iostream>


void MyViewer::resetFlags() {
	for (auto v : mesh.vertices()) {
		mesh.data(v).flags.tagged = false;
		mesh.data(v).flags.temporary_tagged = false;

	}
}
//Generative Design

/// <summary>
/// Calculates the incidences of the input  constraint meshes and the organic part
/// </summary>
void MyViewer::calculateIncidence() {
	//reset incidence data
	for (auto v : mesh.vertices()) {
		mesh.data(v).I.clear();
	}
	// Compute mean values using dihedral angles
	for (auto v : mesh.vertices()) {
		int i = 0;
		for (auto& solid : constraint_solids) {
			if (solid->isInside(mesh.point(v)))
				mesh.data(v).I.push_back(i);
			i++;


		}
	}
	updateMesh();
	update();
}

/// <summary>
/// Remeshing the border between the organic and incident part
/// </summary>
/// <param name="L">The desired edge length in teh remeshed area</param>
/// <param name="iterations">Number of iterations of the remeeshing</param>
void MyViewer::reMeshOrganicBoundaries(double L, int iterations) {
	

	//Calculate the average edge length
	//Just for debugging purposes
	double sum = 0.0;
	int n = 0;
	//Calculate avg edge length
	for (MyMesh::EdgeIter it = mesh.edges_begin(); it != mesh.edges_end(); ++it) {
		sum += mesh.calc_edge_length(*it);
		n++;
	}
	qDebug() << sum / n;

	


	//Do the given iteration number of remeshing
	for (size_t i = 0; i < iterations; i++)
	{
		calculateIncidence();
		reMeshTagvertices();
		reMeshEdgeLength(L);
		reMeshVertexValences();
		updateMesh();
		reMeshVertexPositions();
	}
	updateMesh();
	update();
	calculateIncidence();
	reMeshTagvertices();


}
/// <summary>
/// First step of remeshing
/// Split the edges that are too short for the given L parameter
/// COllapse the edges that are too long for the given L parameter
/// </summary>
/// <param name="L">The desired edge length</param>
void MyViewer::reMeshEdgeLength(double L) {
	bool split = true;
	bool collapse = true;

	//Splitting
	while (split) {
		split = false;
		//Find one edge that needs splitting
		//Then split the edge and start over again
		for (MyMesh::EdgeIter it = mesh.edges_begin(); it != mesh.edges_end(); ++it) {
			if (!checkIfEdgeTagged(it.handle()))
				continue;

			double length = mesh.calc_edge_length(*it);

			if (length > 4.0 / 3.0 * L) {
				const MyViewer::MyTraits::Point& newPoint =
					(mesh.point(mesh.from_vertex_handle(mesh.halfedge_handle(*it, 0)))
						+
						mesh.point(mesh.to_vertex_handle(mesh.halfedge_handle(*it, 0))))
					/ 2;
				MyMesh::VertexHandle tempVh = mesh.split_copy(*it, newPoint);
				mesh.data(tempVh).flags.tagged = true;
				split = true;
				mesh.garbage_collection();
				break;
			}
		}
	}
	//Find one edge that needs collapsing
	//Then collapse the edge and start over again
	while (collapse) {
		collapse = false;
		for (MyMesh::EdgeIter it = mesh.edges_begin(); it != mesh.edges_end(); ++it) {
			if (!checkIfEdgeTagged(it.handle()))
				continue;
			double length = mesh.calc_edge_length(*it);
			if (length < 4.0 / 5.0 * L && !mesh.is_boundary(*it)) {
				MyMesh::HalfedgeHandle h = mesh.halfedge_handle(*it, 0);
				MyMesh::VertexHandle fromVh = mesh.from_vertex_handle(h);
				MyMesh::VertexHandle toVh = mesh.to_vertex_handle(h);

				if (mesh.is_collapse_ok(h) && !mesh.status(h).deleted() && !mesh.is_boundary(fromVh) && !mesh.is_boundary(toVh)) {
					MyMesh::Point newPos = (mesh.point(fromVh) + mesh.point(toVh)) / 2.0f;
					mesh.set_point(fromVh, newPos);
					mesh.set_point(toVh, newPos);
					mesh.collapse(h);
					mesh.garbage_collection();

					collapse = true;
					break;
				}
			}
		}
	}
	mesh.garbage_collection();
}

/// <summary>
/// Tag the 2 ring vertices of the ones that are on the boundary of the organic and incident regions
/// </summary>
void MyViewer::reMeshTagvertices() {
	resetFlags();
	//Tag vertices on the boundary of the 2 regions
	for (auto v : mesh.vertices()) {
		//mesh.status(v).set_tagged2(true);
		if (mesh.data(v).I.size() > 0) {
			for (MyMesh::VertexVertexIter vv_iter = mesh.vv_iter(v); vv_iter.is_valid(); vv_iter++) {
				if (mesh.data(*vv_iter).I.size() == 0) {
					mesh.data(v).flags.tagged = true;
					break;
				}
			}

		}

	}
	//Tag the 2 ring neighbours of the previously tagged vertices
	for (size_t i = 0; i < 2; i++)
	{
		for (auto v : mesh.vertices()) {
			if (mesh.data(v).flags.tagged) {
				for (MyMesh::VertexVertexIter vv_iter = mesh.vv_iter(v); vv_iter.is_valid(); vv_iter++) {
					mesh.data(*vv_iter).flags.temporary_tagged = true;
				}
			}
		}
		for (auto v : mesh.vertices()) {
			if (mesh.data(v).flags.temporary_tagged) {
				mesh.data(v).flags.temporary_tagged = false;
				mesh.data(v).flags.tagged = true;
			}
		}
	}
}
/// <summary>
/// Second step of remeshing
/// For every wuad pefrofm a flip if we get better valences that way
/// </summary>
void MyViewer::reMeshVertexValences() {
	bool flipped = true;
	while (flipped) {

		flipped = false;
		for (MyMesh::EdgeIter it = mesh.edges_begin(); it != mesh.edges_end(); ++it) {
			if (!checkIfEdgeTagged(it.handle()))
				continue;
			if (!mesh.is_boundary(*it)) {
				double squaredDifference = 0;
				MyMesh::HalfedgeHandle h = mesh.halfedge_handle(*it, 0);

				//Calculate valance square difference
				for (MyMesh::FaceVertexIter fv_iter = mesh.fv_iter(mesh.face_handle(h)); fv_iter.is_valid(); fv_iter++)
				{
					double temp;
					if (mesh.is_boundary(fv_iter.handle())) {
						temp = 4.0 - mesh.valence(fv_iter.handle());
					}
					else {
						temp = 6.0 - mesh.valence(fv_iter.handle());
					}
					squaredDifference += temp * temp;
				}
				MyMesh::HalfedgeHandle ho = mesh.opposite_halfedge_handle(h);
				for (MyMesh::FaceVertexIter fv_iter = mesh.fv_iter(mesh.face_handle(ho)); fv_iter.is_valid(); fv_iter++)
				{
					if (fv_iter.handle() != mesh.from_vertex_handle(ho) && fv_iter.handle() != mesh.to_vertex_handle(ho)) {
						double temp;
						if (mesh.is_boundary(fv_iter.handle())) {
							temp = 4.0 - mesh.valence(fv_iter.handle());
						}
						else {
							temp = 6.0 - mesh.valence(fv_iter.handle());
						}
						squaredDifference += temp * temp;
					}

				}

				//Calculate valance square difference wuth flip
				double squaredDifferenceFlip = 0;

				//Calculate valance square difference
				for (MyMesh::FaceVertexIter fv_iter = mesh.fv_iter(mesh.face_handle(h)); fv_iter.is_valid(); fv_iter++)
				{
					double temp;
					double valence = mesh.valence(fv_iter.handle());
					if (mesh.from_vertex_handle(h) == fv_iter.handle() || mesh.to_vertex_handle(h) == fv_iter.handle()) {
						valence--;
					}
					else {
						valence++;
					}
					if (mesh.is_boundary(fv_iter.handle())) {
						temp = 4.0 - valence;
					}
					else {
						temp = 6.0 - valence;
					}

					squaredDifferenceFlip += temp * temp;
				}
				for (MyMesh::FaceVertexIter fv_iter = mesh.fv_iter(mesh.face_handle(ho)); fv_iter.is_valid(); fv_iter++)
				{
					if (fv_iter.handle() != mesh.from_vertex_handle(ho) && fv_iter.handle() != mesh.to_vertex_handle(ho)) {
						double temp;
						if (mesh.is_boundary(*fv_iter)) {
							temp = 4.0 - (mesh.valence(fv_iter.handle()) + 1);
						}
						else {
							temp = 6.0 - (mesh.valence(fv_iter.handle()) + 1);
						}

						squaredDifferenceFlip += temp * temp;
					}

				}
				//If better then flip
				if (squaredDifferenceFlip < squaredDifference) {
					mesh.flip(*it);
					flipped = true;
					break;
				}
			}
		}
	}
}
/// <summary>
/// Third step of remeshing
/// For every vertex move the vertex into the avg of their neighbours projected back to the tangent plane
/// </summary>
void MyViewer::reMeshVertexPositions() {
	for (auto vh : mesh.vertices())
	{
		if (!mesh.data(vh).flags.tagged || mesh.is_boundary(vh))
			continue;
		MyMesh::Point center = MyMesh::Point(0,0,0);
		MyMesh::VertexVertexIter    vv_it;
		int n = 0;
		for (vv_it = mesh.vv_iter(vh); vv_it.is_valid(); ++vv_it) {
			center += mesh.point(*vv_it);
			n++;
		}
		center = center / n;
		MyMesh::Point xminusc = mesh.point(vh) - center;
		double s = mesh.normal(vh)[0]* xminusc[0] + mesh.normal(vh)[1] * xminusc[1] + mesh.normal(vh)[2] * xminusc[2];
		//qDebug() << mesh.normal(vh).length();
		mesh.set_point(vh, center + s * mesh.normal(vh));
	}
}
bool MyViewer::checkIfEdgeTagged(MyMesh::EdgeHandle eh) {
	return checkIfEdgeTagged(mesh.halfedge_handle(eh, 0));
}
bool MyViewer::checkIfEdgeTagged(MyMesh::HalfedgeHandle hh) {
	return (mesh.data(mesh.from_vertex_handle(hh)).flags.tagged && mesh.data(mesh.to_vertex_handle(hh)).flags.tagged);
}
void MyViewer::reMeshSmoothing(int iterations) {
	for (size_t i = 0; i < iterations; i++)
	{
		reMeshSmoothingIteration();
	}
}
/// <summary>
/// Smooth the edges between the organic and incident regions
/// This way the angle differneces will decrease
/// </summary>
void MyViewer::reMeshSmoothingIteration() {
	for (auto v : mesh.vertices()) {
		mesh.data(v).flags.temporary_tagged = false;
	}
	MyMesh::VertexHandle vhPrev;
	bool prevSet = false;
	bool nextSet = false;
	MyMesh::VertexHandle vhNext;
	for (auto v : mesh.vertices()) {
		if (checkIfBetweenRegions(v)) {

			for (MyMesh::VertexVertexIter vv_iter = mesh.vv_iter(v); vv_iter.is_valid(); vv_iter++) {
				if (checkIfBetweenRegions(vv_iter.handle())) {
					if (!prevSet) {
						vhPrev = vv_iter.handle();
						prevSet = true;
					}
					else {
						vhNext = vv_iter.handle();
						nextSet = true;
					}
				}
			}
			if (prevSet && nextSet) {
				MyMesh::Point pPrev = mesh.point(vhPrev);
				MyMesh::Point pNext = mesh.point(vhNext);
				MyMesh::Point p = mesh.point(v);
				MyMesh::Point L = 0.5 * (pNext - p) + 0.5 * (pPrev - p);
				mesh.data(v).newPos = p + 0.5 * L;

				mesh.data(v).flags.temporary_tagged = true;
			}
			
		}
		nextSet = false;
		prevSet = false;

	}
	for (auto v : mesh.vertices()) {
		if (mesh.data(v).flags.temporary_tagged) {
			mesh.set_point(v, mesh.data(v).newPos);
		}
			
	}
	updateMesh();
	update();
}
bool MyViewer::checkIfBetweenRegions(MyMesh::VertexHandle vh) {
	bool incident = mesh.data(vh).I.size() > 0;
	for (MyMesh::VertexVertexIter vv_iter = mesh.vv_iter(vh); vv_iter.is_valid(); vv_iter++) {
		if (mesh.data(*vv_iter).I.size() == 0 && incident) {

			return true;
		}
	}
	return false;
}
void MyViewer::catmullClark() {
	catmul_clark_subdivider(1);
	updateMesh(false);
}
void MyViewer::quadrangulate() {
	makeEvenTriangles();
	makeQuadDominant();
	makePureQuad();
}
/// <summary>
/// For the quadrangulation we need an even number of triangles
/// If their number is not equal perform a split at the longest boundary edge
/// </summary>
void MyViewer::makeEvenTriangles() {
	qDebug() << mesh.n_faces();
	if (mesh.n_faces() % 2 != 0) {
		double longest = 0;
		MyMesh::EdgeHandle eh;
		for (MyMesh::EdgeIter it = mesh.edges_begin(); it != mesh.edges_end(); ++it) {
			if (mesh.calc_edge_length(*it) > longest && mesh.is_boundary(*it)) {
				longest = mesh.calc_edge_length(*it);
				eh = it.handle();
			}			
		}
		const MyViewer::MyTraits::Point& newPoint =
			(mesh.point(mesh.from_vertex_handle(mesh.halfedge_handle(eh, 0)))
				+
				mesh.point(mesh.to_vertex_handle(mesh.halfedge_handle(eh, 0))))
			/ 2;
		MyMesh::VertexHandle tempVh = mesh.split_copy(eh, newPoint);
		mesh.garbage_collection();
		updateMesh();
		update();
	}
}
/// <summary>
/// First step of quadrangulation
/// make a quad dominant mesh by deleting the best edges from the triangles
/// </summary>
void MyViewer::makeQuadDominant() {
	//Score edges by there "Squareness" (from dot products - angles)

}

void MyViewer::makePureQuad() {

}