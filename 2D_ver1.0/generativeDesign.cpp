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

	//Collapsing
	while (collapse) {
		collapse = false;
		//Find one edge that needs collapsing
		//Then collapse the edge and start over again
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
	calculateSquarness();
}
void MyViewer::calculateSquarness() {
	resetEdgeProps();
	std::vector<MyMesh::EdgeHandle> edges;
	for (MyMesh::EdgeIter it = mesh.edges_begin(); it != mesh.edges_end(); ++it) {
		if (mesh.is_boundary(*it))
			continue;
		mesh.data(it.handle()).tagged = false;
		MyMesh::HalfedgeHandle h0 = mesh.halfedge_handle(*it, 0);
		MyMesh::HalfedgeHandle h1 = mesh.opposite_halfedge_handle(h0);


		MyMesh::VertexHandle v0 = mesh.from_vertex_handle(h0);
		MyMesh::VertexHandle v1 = mesh.to_vertex_handle(h0);

		MyMesh::FaceHandle f0 = mesh.face_handle(h0);
		MyMesh::FaceHandle f1 = mesh.face_handle(h1);

		MyMesh::VertexHandle p0;
		MyMesh::VertexHandle p1;


		for (MyMesh::FaceVertexIter fv_iter = mesh.fv_iter(f0); fv_iter.is_valid(); fv_iter++) {
			if (fv_iter.handle() != v0 && fv_iter.handle() != v1) {
				p0 = fv_iter.handle();
			}
		}

		for (MyMesh::FaceVertexIter fv_iter = mesh.fv_iter(f1); fv_iter.is_valid(); fv_iter++) {
			if (fv_iter.handle() != v0 && fv_iter.handle() != v1) {
				p1 = fv_iter.handle();
			}
		}

		OpenMesh::Vec3d p0v0 = (mesh.point(v0) - mesh.point(p0)).normalized();
		OpenMesh::Vec3d p0v1 = (mesh.point(v1) - mesh.point(p0)).normalized();
		OpenMesh::Vec3d p1v0 = (mesh.point(v0) - mesh.point(p1)).normalized();
		OpenMesh::Vec3d p1v1 = (mesh.point(v1) - mesh.point(p1)).normalized();

		double temp = abs(dot(p0v0, p0v1)) + abs(dot(p1v0, p1v1));
		mesh.data(it.handle()).squarness = temp;

		bool inserted = false;
		if (edges.size() > 0) {
			for (std::vector<MyMesh::EdgeHandle>::iterator vec_it = edges.begin(); vec_it != edges.end(); ++vec_it) {
				if (mesh.data(*vec_it).squarness > temp) {
					edges.insert(vec_it, it.handle());
					inserted = true;
					break;
				}
			}
			if (!inserted)
				edges.push_back(it.handle());
		}
		else {
			edges.push_back(it.handle());
		}


	}
	int count = 0;
	for (auto eh : edges)
	{
		if (canDelete(eh)) {
			mesh.data(eh).tagged = true;
			MyMesh::HalfedgeHandle h0 = mesh.halfedge_handle(eh, 0);


			MyMesh::VertexHandle v0 = mesh.from_vertex_handle(h0);
			MyMesh::VertexHandle v1 = mesh.to_vertex_handle(h0);

			mesh.data(v0).edgeTagged = true;
			mesh.data(v1).edgeTagged = true;
			count++;
		}
	}
	qDebug() << "Deleted " << count;
	counter = 0;
	//deleteEdges(edges);
}
void MyViewer::deleteEdges() {

	bool deleted = true;
	while (deleted) {
		deleted = false;
		for (MyMesh::EdgeIter it = mesh.edges_begin(); it != mesh.edges_end(); ++it) {
			if (mesh.data(it.handle()).tagged ) {
				MyMesh::HalfedgeHandle h0 = mesh.halfedge_handle(*it, 0);
				MyMesh::HalfedgeHandle h1 = mesh.opposite_halfedge_handle(h0);


				MyMesh::VertexHandle v0 = mesh.from_vertex_handle(h0);
				MyMesh::VertexHandle v1 = mesh.to_vertex_handle(h0);

				MyMesh::FaceHandle f0 = mesh.face_handle(h0);
				MyMesh::FaceHandle f1 = mesh.face_handle(h1);

				MyMesh::VertexHandle p0;
				MyMesh::VertexHandle p1;

				std::vector<MyMesh::HalfedgeHandle> temphh;

				for (MyMesh::FaceVertexIter fv_iter = mesh.fv_iter(f0); fv_iter.is_valid(); fv_iter++) {
					if (fv_iter.handle() != v0 && fv_iter.handle() != v1) {
						p0 = fv_iter.handle();
					}
				}
				for (MyMesh::FaceVertexIter fv_iter = mesh.fv_iter(f1); fv_iter.is_valid(); fv_iter++) {
					if (fv_iter.handle() != v0 && fv_iter.handle() != v1) {
						p1 = fv_iter.handle();
					}
				}
				for (MyMesh::FaceHalfedgeIter fh_iter = mesh.fh_iter(f0); fh_iter.is_valid(); fh_iter++) {
					if (fh_iter.handle() != h0 && fh_iter.handle() != h1)
						temphh.push_back(fh_iter.handle());
				}
				for (MyMesh::FaceHalfedgeIter fh_iter = mesh.fh_iter(f1); fh_iter.is_valid(); fh_iter++) {
					if (fh_iter.handle() != h0 && fh_iter.handle() != h1)
						temphh.push_back(fh_iter.handle());
				}
				//if (counter == 138) {
				//	qDebug() << "ASD";
					delete_edge(it.handle(),false);
					mesh.garbage_collection();
				//}
				//else {
				//	mesh.data(it.handle()).tagged = false;
				//}


				//qDebug() << temphh.size();
				//std::vector<MyMesh::VertexHandle>  face_vhandles;
				//face_vhandles.push_back(v0);
				//face_vhandles.push_back(p1);
				//face_vhandles.push_back(v1);
				//face_vhandles.push_back(p0);
				//mesh.add_face(face_vhandles);
				//MyMesh::FaceHandle temp_fh = mesh.new_face();
				//for (MyMesh::FaceHalfedgeIter fh_iter = mesh.fh_iter(temp_fh); fh_iter.is_valid(); fh_iter++) {
				//	qDebug() << fh_iter.handle().idx();
				//}
				//qDebug() << "OGNIRONGRL";
				//for (auto h : temphh) {
				//	mesh.set_face_handle(h, temp_fh);
				//	mesh.set_halfedge_handle(temp_fh, h);
				//}
				//for (MyMesh::FaceHalfedgeIter fh_iter = mesh.fh_iter(temp_fh); fh_iter.is_valid(); fh_iter++) {
				//	qDebug() << fh_iter.handle().idx();
				//}
				//MyMesh::FaceHandle temp =  mesh.add_face(face_vhandles);
				/*for (MyMesh::FaceVertexIter fv_iter = mesh.fv_iter(temp); fv_iter.is_valid(); fv_iter++) {
					qDebug() << fv_iter.handle().idx();
				}*/
				deleted = true;
				//qDebug() <<counter++;
				break;
			}
		}

	}
	updateMesh();
	update();
}
bool MyViewer::canDelete(MyMesh::EdgeHandle eh) {
	MyMesh::HalfedgeHandle h0 = mesh.halfedge_handle(eh, 0);
	MyMesh::HalfedgeHandle h1 = mesh.opposite_halfedge_handle(h0);

	MyMesh::FaceHandle f0 = mesh.face_handle(h0);
	MyMesh::FaceHandle f1 = mesh.face_handle(h1);

	for (MyMesh::FaceEdgeIter fe_iter = mesh.fe_iter(f0); fe_iter.is_valid(); fe_iter++) {
		if (mesh.data(fe_iter.handle()).tagged)
			return false;
	}
	for (MyMesh::FaceEdgeIter fe_iter = mesh.fe_iter(f1); fe_iter.is_valid(); fe_iter++) {
		if (mesh.data(fe_iter.handle()).tagged)
			return false;
	}
	return true;
}
void MyViewer::makePureQuad() {

}
void MyViewer::delete_edge(MyMesh::EdgeHandle _eh, bool debugData) {
	MyMesh::HalfedgeHandle h0 = mesh.halfedge_handle(_eh, 0);
	MyMesh::HalfedgeHandle h1 = mesh.opposite_halfedge_handle(h0);

	MyMesh::HalfedgeHandle h0prev = mesh.prev_halfedge_handle(h0);
	MyMesh::HalfedgeHandle h0next = mesh.next_halfedge_handle(h0);
	MyMesh::HalfedgeHandle h1prev = mesh.prev_halfedge_handle(h1);
	MyMesh::HalfedgeHandle h1next = mesh.next_halfedge_handle(h1);

	MyMesh::VertexHandle v0 = mesh.from_vertex_handle(h0);
	MyMesh::VertexHandle v1 = mesh.to_vertex_handle(h0);

	MyMesh::FaceHandle f0 = mesh.face_handle(h0);
	MyMesh::FaceHandle f1 = mesh.face_handle(h1);

	MyMesh::VertexHandle p0;
	MyMesh::VertexHandle p1;

	for (MyMesh::FaceVertexIter fv_iter = mesh.fv_iter(f0); fv_iter.is_valid(); fv_iter++) {
		if (fv_iter.handle() != v0 && fv_iter.handle() != v1) {
			p0 = fv_iter.handle();
		}
	}

	for (MyMesh::FaceVertexIter fv_iter = mesh.fv_iter(f1); fv_iter.is_valid(); fv_iter++) {
		if (fv_iter.handle() != v0 && fv_iter.handle() != v1) {
			p1 = fv_iter.handle();
		}
	}


	mesh.status(f0).set_deleted(true);
	mesh.status(h0).set_deleted(true);
	mesh.status(h1).set_deleted(true);
	mesh.status(_eh).set_deleted(true);


	//set face attributes
	mesh.set_halfedge_handle(f1, h1next);

	//set half_edge attributes
	mesh.set_face_handle(h0prev, f1);
	mesh.set_face_handle(h0next, f1);

	mesh.set_next_halfedge_handle(h1prev, h0next);
	mesh.set_next_halfedge_handle(h0prev, h1next);

	mesh.set_prev_halfedge_handle(h1next, h0prev);
	mesh.set_prev_halfedge_handle(h0next, h1prev);

	//set vertex attributes
	mesh.set_halfedge_handle(v0, h1next);
	mesh.set_halfedge_handle(v1, h0next);

	
	if (debugData) {
		int count = 0;
		for (auto v : mesh.fv_range(f1)) {
			if (v == v0)
				qDebug() << "v0";
			if (v == v1)
				qDebug() << "v1";
			if (v == p0)
				qDebug() << "p0";
			if (v == p1)
				qDebug() << "p1";
			count++;
		}
		MyMesh::HalfedgeHandle temph = mesh.halfedge_handle(f1);

		if (temph == h0)
			qDebug() << "h0";

		if (temph == h1)
			qDebug() << "h1";

		if (temph == h0prev)
			qDebug() << "h0prev";

		if (temph == h0next)
			qDebug() << "h0next";

		if (temph == h1prev)
			qDebug() << "h1prev";

		if (temph == h1next)
			qDebug() << "h1next";
		qDebug() << "count:" << count;
	}


}