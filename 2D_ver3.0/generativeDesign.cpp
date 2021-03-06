#include "MyViewer.h"
#include <OpenMesh/Tools/Smoother/JacobiLaplaceSmootherT.hh>
#include <iostream>



void MyViewer::resetFlags() {
	for (auto v : mesh.vertices()) {
		mesh.data(v).flags.tagged = false;
		mesh.data(v).flags.temporary_tagged = false;

	}
}

void MyViewer::resetEdgeProps() {
	for (auto e : mesh.edges()) {
		mesh.data(e).squarness = 0;
		mesh.data(e).tagged = false;

	}
}


void MyViewer::decimate() {
	OpenMesh::Decimater::DecimaterT<MyMesh> decimater(mesh);
	OpenMesh::Decimater::ModRoundnessT<MyMesh>::Handle mod;
	decimater.add(mod);
	decimater.initialize();
	decimater.decimate_to_faces(10);
	mesh.garbage_collection();
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
void MyViewer::reMeshAll(int iterations) {
	//Calculate the average edge length
	double sum = 0.0;
	int n = 0;
	//Calculate avg edge length
	for (MyMesh::EdgeIter it = mesh.edges_begin(); it != mesh.edges_end(); ++it) {
		sum += mesh.calc_edge_length(*it);
		n++;
	}
	double L = (sum / n);
	MyMesh::HalfedgeIter he_it, he_end = mesh.halfedges_end();
	//qDebug() << "L: " << L << "\n";
	//Do the given iteration number of remeshing
	for (size_t i = 0; i < iterations; i++)
	{
		//qDebug() << "Tag vertices \n";
		reMeshTagvertices(true);

		//qDebug() << "reMeshEdgeLength \n";
		reMeshEdgeLength(L);

		//qDebug() << "reMeshVertexValences \n";
		reMeshVertexValences();
		//qDebug() << "updateMesh \n";
		updateMesh();
		//qDebug() << "reMeshVertexPositions \n";
		reMeshVertexPositions(false);
	}

	resetFlags();
	//qDebug() << "updateMesh \n";
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
		//qDebug() << "TEST\n";
	}

	//Collapsing
	while (collapse) {
		collapse = false;
		//Find one edge that needs collapsing
		//Then collapse the edge and start over again
		for (MyMesh::EdgeIter it = mesh.edges_begin(); it != mesh.edges_end(); ++it) {
			if (!checkIfEdgeTagged(it.handle()))
				continue;
			double length = mesh.calc_edge_length(*it);
			MyMesh::HalfedgeHandle h = mesh.halfedge_handle(*it, 0);
			MyMesh::VertexHandle fromVh = mesh.from_vertex_handle(h);
			MyMesh::VertexHandle toVh = mesh.to_vertex_handle(h);

			if (length < 4.0 / 5.0 * L && (!mesh.is_boundary(*it)||!mesh.data(fromVh).flags.locked)) {


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
void MyViewer::reMeshTagvertices(bool all) {
	resetFlags();
	if (all) {
		for (auto v : mesh.vertices()) {
			mesh.data(v).flags.tagged = true;
		}
	}
	else {
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
void MyViewer::reMeshVertexPositions(bool remeshUnlockedBoundaries) {
	for (auto vh : mesh.vertices())
	{
		if (!mesh.data(vh).flags.tagged || (mesh.is_boundary(vh) && (!remeshUnlockedBoundaries || mesh.data(vh).flags.locked)))
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
void MyViewer::removeIncidentVertices(){
	resetFlags();
	for (auto v : mesh.vertices()) {
		if (mesh.data(v).I.size() > 0 && !checkIfBetweenRegions(v))
			mesh.data(v).flags.tagged = true;
	}
	for (auto v : mesh.vertices()) {
		if (mesh.data(v).flags.tagged)
			mesh.delete_vertex(v);
	}
	mesh.garbage_collection();
	updateMesh();
	update();
}
void MyViewer::catmullClark() {
	catmul_clark_subdivider(1);
	updateMesh(false);
}
void MyViewer::partition() {
	resetFlags();
	resetEdgeProps();
	for (auto e : mesh.edges()) {
		if (mesh.is_boundary(e)) {
			mesh.data(e).tagged = true;

		}
	}
	for (auto v : mesh.vertices()) {
		mesh.data(v).I.clear();
		if (mesh.valence(v) != 4 && !mesh.is_boundary(v)) {
			mesh.data(v).flags.tagged = true;
			for (MyMesh::VertexEdgeIter ve_iter = mesh.ve_iter(v); ve_iter.is_valid(); ve_iter++) {
				mesh.data(ve_iter).tagged = true;
				if (mesh.from_vertex_handle(mesh.halfedge_handle(ve_iter, 0)) != v &&
					!mesh.is_boundary(mesh.from_vertex_handle(mesh.halfedge_handle(ve_iter, 0))))
					mesh.data(mesh.from_vertex_handle(mesh.halfedge_handle(ve_iter, 0))).flags.temporary_tagged = true;
				if (mesh.from_vertex_handle(mesh.halfedge_handle(ve_iter, 1)) != v &&
					!mesh.is_boundary(mesh.from_vertex_handle(mesh.halfedge_handle(ve_iter, 1))))
					mesh.data(mesh.from_vertex_handle(mesh.halfedge_handle(ve_iter, 1))).flags.temporary_tagged = true;

			}
		}
	}
	int idx = 0;
	bool finish = false;
	while (!finish) {
		finish = true;
		for (auto v : mesh.vertices()) {
			if (mesh.data(v).flags.temporary_tagged && !mesh.data(v).flags.tagged) {
				std::vector<MyMesh::EdgeHandle>edges;
				for (MyMesh::VertexEdgeIter ve_iter = mesh.ve_iter(v); ve_iter.is_valid(); ve_iter++) {
					edges.push_back(ve_iter.handle());
				}
				if (mesh.data(edges[0]).tagged && !mesh.data(edges[2]).tagged &&
					!(mesh.data(edges[1]).tagged && mesh.data(edges[3]).tagged)) {
					mesh.data(edges[2]).tagged = true;
					if (!mesh.is_boundary(mesh.from_vertex_handle(mesh.halfedge_handle(edges[2], 0))))
						mesh.data(mesh.from_vertex_handle(mesh.halfedge_handle(edges[2], 0))).flags.temporary_tagged2 = true;
					if (!mesh.is_boundary(mesh.from_vertex_handle(mesh.halfedge_handle(edges[2], 1))))
						mesh.data(mesh.from_vertex_handle(mesh.halfedge_handle(edges[2], 1))).flags.temporary_tagged2 = true;
					mesh.data(v).flags.temporary_tagged2 = false;
					continue;
				}
				if (mesh.data(edges[2]).tagged && !mesh.data(edges[0]).tagged &&
					!(mesh.data(edges[1]).tagged && mesh.data(edges[3]).tagged)) {
					mesh.data(edges[0]).tagged = true;
					if (!mesh.is_boundary(mesh.from_vertex_handle(mesh.halfedge_handle(edges[0], 0))))
						mesh.data(mesh.from_vertex_handle(mesh.halfedge_handle(edges[0], 0))).flags.temporary_tagged2 = true;
					if (!mesh.is_boundary(mesh.from_vertex_handle(mesh.halfedge_handle(edges[0], 1))))
						mesh.data(mesh.from_vertex_handle(mesh.halfedge_handle(edges[0], 1))).flags.temporary_tagged2 = true;
					mesh.data(v).flags.temporary_tagged2 = false;
					continue;

				}
				if (mesh.data(edges[1]).tagged && !mesh.data(edges[3]).tagged &&
					!(mesh.data(edges[0]).tagged && mesh.data(edges[2]).tagged)) {
					mesh.data(edges[3]).tagged = true;
					if (!mesh.is_boundary(mesh.from_vertex_handle(mesh.halfedge_handle(edges[3], 0))))
						mesh.data(mesh.from_vertex_handle(mesh.halfedge_handle(edges[3], 0))).flags.temporary_tagged2 = true;
					if (!mesh.is_boundary(mesh.from_vertex_handle(mesh.halfedge_handle(edges[3], 1))))
						mesh.data(mesh.from_vertex_handle(mesh.halfedge_handle(edges[3], 1))).flags.temporary_tagged2 = true;
					mesh.data(v).flags.temporary_tagged2 = false;
					continue;

				}
				if (mesh.data(edges[3]).tagged && !mesh.data(edges[1]).tagged &&
					!(mesh.data(edges[0]).tagged && mesh.data(edges[2]).tagged)) {
					mesh.data(edges[1]).tagged = true;
					if (!mesh.is_boundary(mesh.from_vertex_handle(mesh.halfedge_handle(edges[1], 0))))
						mesh.data(mesh.from_vertex_handle(mesh.halfedge_handle(edges[1], 0))).flags.temporary_tagged2 = true;
					if (!mesh.is_boundary(mesh.from_vertex_handle(mesh.halfedge_handle(edges[1], 1))))
						mesh.data(mesh.from_vertex_handle(mesh.halfedge_handle(edges[1], 1))).flags.temporary_tagged2 = true;
					mesh.data(v).flags.temporary_tagged2 = false;
					continue;

				}
			}
		}
		int count = 0;
		for (auto v : mesh.vertices()) {
			mesh.data(v).flags.temporary_tagged = false;
			if (mesh.data(v).flags.temporary_tagged2) {
				mesh.data(v).flags.temporary_tagged = true;
				mesh.data(v).flags.temporary_tagged2 = false;
				finish = false;
				count++;
			}
		}
		qDebug() << count;
	}
	updateMesh();
	update();
}

void MyViewer::initiatePDE() {
	int i = 0;
	for (auto v : mesh.vertices()) {
		mesh.data(v).idx = i++;
		mesh.data(v).u = OpenMesh::Vec3d(0.0);
	}

	initiateBoundaryConstraints();
}

void MyViewer::initiateBoundaryConstraints() {
	for (auto v : mesh.vertices()) {
		if (mesh.is_boundary(v)) {
			MyMesh::Point neighbour1;
			MyMesh::Point neighbour2;
			bool found = false;
			for (MyMesh::VertexVertexIter vv_iter = mesh.vv_iter(v); vv_iter.is_valid(); vv_iter++) {
				if (mesh.is_boundary(vv_iter.handle())) {
					if (!found) {
						neighbour1 = mesh.point(vv_iter.handle());
						found = true;
					}
					else {
						neighbour2 = mesh.point(vv_iter.handle());
					}
				}
			}
			double randomNumber = ((double)rand()/ RAND_MAX )*2- 1;
			double randomNumber2 = ((double)rand() / RAND_MAX) * 2 - 1;
			mesh.data(v).u = OpenMesh::Vec3d(randomNumber, randomNumber2,0.0);
		}
	}
	//for (auto v : mesh.vertices()) {
	//	if (mesh.data(v).u.length() > 0.0) {
	//		double theta = (atan2(mesh.data(v).u[0], mesh.data(v).u[1]) + 3.14159)/ 4.0;
	//		printf("%f\n", theta);
	//	}
	//}
}