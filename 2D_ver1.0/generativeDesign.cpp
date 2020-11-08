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
void MyViewer::reMeshAll(int iterations) {
	//Calculate the average edge length
	double sum = 0.0;
	int n = 0;
	//Calculate avg edge length
	for (MyMesh::EdgeIter it = mesh.edges_begin(); it != mesh.edges_end(); ++it) {
		sum += mesh.calc_edge_length(*it);
		n++;
	}
	double L = sum / n;

	//Do the given iteration number of remeshing
	for (size_t i = 0; i < iterations; i++)
	{
		reMeshTagvertices(true);
		reMeshEdgeLength(L);
		reMeshVertexValences();
		updateMesh();
		reMeshVertexPositions();
	}

	resetFlags();
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
	int sum4 = 0;
	int sum3 = 0;
	for (auto f : mesh.faces())
	{
		if (faceSides(f) == 4)
			sum4++;
		if (faceSides(f) == 3)
			sum3++;

	}
	qDebug() << "Before 3side: " << sum3 << " 4 side: " << sum4;
	makePureQuad(190);
	sum4 = 0;
	sum3 = 0;
	for (auto f : mesh.faces())
	{
		if (faceSides(f) == 4)
			sum4++;
		if (faceSides(f) == 3)
			sum3++;

	}
	qDebug() << "After 3side:" << sum3 << " 4 side: " << sum4;
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

		OpenMesh::Vec3d v0p0 = (mesh.point(p0) - mesh.point(v0)).normalized();
		OpenMesh::Vec3d v0p1 = (mesh.point(p1) - mesh.point(v0)).normalized();
		OpenMesh::Vec3d v1p0 = (mesh.point(p0) - mesh.point(v1)).normalized();
		OpenMesh::Vec3d v1p1 = (mesh.point(p1) - mesh.point(v1)).normalized();

		double temp = abs(dot(p0v0, p0v1)) + abs(dot(p1v0, p1v1)) + abs(dot(v0p0, v0p1)) + abs(dot(v1p0, v1p1));
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
	deleteEdges();
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

				delete_edge(it.handle());
				mesh.garbage_collection();
				deleted = true;
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
void MyViewer::makePureQuad(int idx) {
	int idx2 = 0;
	int idx3 = 0;
	//int count = 0;
	for (MyMesh::FaceIter it = mesh.faces_begin(); it != mesh.faces_end(); ++it) {
		resetFaceFlags();
		if (faceSides(it.handle()) == 3) {
			//idx2++;
			//if (++idx2 < idx) continue;
			//mesh.data(it).tagged2 = true;
			///BFS - START
			/// ############
			std::queue<MyMesh::FaceHandle> prevFaces;
			bool foundOthertriangle = false;
			MyMesh::FaceHandle pairFace;

			prevFaces.push(it.handle());
			mesh.data(it).hasPrev = true;
			MyMesh::FaceHandle tempFh;
			while (prevFaces.size() > 0 &&!foundOthertriangle) {
				tempFh = prevFaces.front();
				prevFaces.pop();

				for (MyMesh::FaceFaceIter ff_iter = mesh.ff_iter(tempFh); ff_iter.is_valid(); ff_iter++) {
					if (!mesh.data(ff_iter).hasPrev) {

						prevFaces.push(ff_iter.handle());
						mesh.data(ff_iter).hasPrev = true;
						mesh.data(ff_iter).prevFace = tempFh;

						if (faceSides(ff_iter.handle()) == 3) {
							foundOthertriangle = true;
							//mesh.data(ff_iter).tagged2 = true;
							pairFace = ff_iter.handle();
							break;
						}
					}
				}
				
			}
			///BFS - END
			/// ############
			if (it.handle().idx() == 1388 || it.handle().idx() == 1417) {
				mesh.data(it).tagged2 = true;
				continue;

			}
			while (mesh.data(pairFace).prevFace != it.handle()) {

				mesh.data(pairFace).tagged = true;
				MyMesh::FaceHandle tempFace = mesh.data(pairFace).prevFace;
				MyMesh::EdgeHandle eprev = getCommonEdge(tempFace, pairFace);
				mesh.data(eprev).tagged = true;
				MyMesh::EdgeHandle enext = getCommonEdge(tempFace, mesh.data(tempFace).prevFace);

				MyMesh::VertexHandle vprev0, vprev1, vprev2, vmiddle0, vmiddle1, vnext;
				vprev0 = mesh.from_vertex_handle(mesh.halfedge_handle(eprev, 0));
				vprev1 = mesh.from_vertex_handle(mesh.halfedge_handle(eprev, 1));
				if (!isVertexOnFace(vprev0, mesh.data(tempFace).prevFace)) {
					MyMesh::VertexHandle temp = vprev0;
					vprev0 = vprev1;
					vprev1 = temp;
				}
				for (MyMesh::FaceVertexIter fv_iter = mesh.fv_iter(pairFace); fv_iter.is_valid(); fv_iter++) {
					if (fv_iter.handle() != vprev0 && fv_iter.handle() != vprev1)
						vprev2 = fv_iter.handle();
				}
				for (MyMesh::FaceHalfedgeIter fh_iter = mesh.fh_iter(tempFace); fh_iter.is_valid(); fh_iter++) {
					if (mesh.from_vertex_handle(fh_iter) == vprev0 && mesh.to_vertex_handle(fh_iter)   != vprev1) {
						vmiddle0 = mesh.to_vertex_handle(fh_iter);
					}
					if (mesh.to_vertex_handle(fh_iter) == vprev0 && mesh.from_vertex_handle(fh_iter) != vprev1) {
						vmiddle0 = mesh.from_vertex_handle(fh_iter);
					}
					if (mesh.from_vertex_handle(fh_iter) == vprev1 && mesh.to_vertex_handle(fh_iter) != vprev0) {
						vmiddle1 = mesh.to_vertex_handle(fh_iter);
					}
					if (mesh.to_vertex_handle(fh_iter) == vprev1 && mesh.from_vertex_handle(fh_iter) != vprev0) {
						vmiddle1 = mesh.from_vertex_handle(fh_iter);
					}
				}
				for (MyMesh::FaceHalfedgeIter fh_iter = mesh.fh_iter(mesh.data(tempFace).prevFace); fh_iter.is_valid(); fh_iter++) {
					if (mesh.from_vertex_handle(fh_iter) == vmiddle0 && mesh.to_vertex_handle(fh_iter) != vmiddle1 && mesh.to_vertex_handle(fh_iter) != vprev0) {
						vnext = mesh.to_vertex_handle(fh_iter);
					}
					if (mesh.to_vertex_handle(fh_iter) == vmiddle0 && mesh.from_vertex_handle(fh_iter) != vmiddle1 && mesh.from_vertex_handle(fh_iter) != vprev0) {
						vnext = mesh.from_vertex_handle(fh_iter);
					}
				}

				MyMesh::FaceHandle nextPrevFace = mesh.data(tempFace).prevFace;
				bool opposite = !hasCommonVertex(pairFace, nextPrevFace);
				MyMesh::FaceHandle newFace = delete_edge(eprev);
			

				mesh.data(newFace).prevFace = nextPrevFace;

				if (!opposite) {
					std::vector<MyMesh::VertexHandle> temp = { vprev2, vprev0, vnext, vmiddle0 };
					std::vector<MyMesh::VertexHandle> temp2 = { vmiddle1, vprev0, vnext, vmiddle0 };
					MyMesh::FaceHandle tempFace;
					if (getAngleSum(temp) > getAngleSum(temp2))
						tempFace = add_edge(newFace, vmiddle0, vprev2);
					else
						tempFace = add_edge(newFace, vprev0, vmiddle1);
					mesh.data(tempFace).prevFace = nextPrevFace;
					pairFace = tempFace;

				}
				else {
						
					resetFlags();
					

					tempFace = add_edge(newFace, vmiddle0, vprev1);

					mesh.data(tempFace).prevFace = nextPrevFace;
					pairFace = tempFace;
				}

				updateMesh();

			}

			MyMesh::FaceHandle tempFace = mesh.data(pairFace).prevFace;
			MyMesh::EdgeHandle eprev = getCommonEdge(tempFace, pairFace);
			MyMesh::FaceHandle temporary = delete_edge(eprev);
			//mesh.data(temporary).tagged2 = true;
			mesh.garbage_collection();
			updateMesh();
			//break;
		}

	}
}
bool MyViewer::isVertexOnFace(MyMesh::VertexHandle vh, MyMesh::FaceHandle fh) {
	for (MyMesh::FaceVertexIter fv_iter = mesh.fv_iter(fh); fv_iter.is_valid(); fv_iter++) {
		if (fv_iter.handle() == vh)
			return true;
	}
	return false;
}
MyViewer::MyMesh::EdgeHandle MyViewer::getCommonEdge(MyMesh::FaceHandle fh1, MyMesh::FaceHandle fh2) {
	for (MyMesh::FaceEdgeIter fe_iter = mesh.fe_iter(fh1); fe_iter.is_valid(); fe_iter++) {
		if (mesh.face_handle(mesh.halfedge_handle(fe_iter.handle(), 0)) == fh2 ||
			mesh.face_handle(mesh.halfedge_handle(fe_iter.handle(), 1)) == fh2) {
			return fe_iter.handle();
		}
	}
	throw 0;
}
bool MyViewer::hasCommonVertex(MyMesh::FaceHandle fh1, MyMesh::FaceHandle fh2) {
	for (MyMesh::FaceVertexIter fv_iter = mesh.fv_iter(fh1); fv_iter.is_valid(); fv_iter++) {
		if (isVertexOnFace(fv_iter.handle(), fh2))
			return true;
	}
	return false;
}
bool MyViewer::hasCommonEdge(MyMesh::VertexHandle vh1, MyMesh::VertexHandle vh2) {
	for (MyMesh::VertexVertexIter vv_iter = mesh.vv_iter(vh1); vv_iter.is_valid(); vv_iter++) {
		if (vv_iter.handle() == vh2)
			return true;
	}
	return false;
}
MyViewer::MyMesh::EdgeHandle MyViewer::getCommonEdge(MyMesh::VertexHandle vh1, MyMesh::VertexHandle vh2) {
	for (MyMesh::VertexIHalfedgeIter vih_iter = mesh.vih_iter(vh1); vih_iter.is_valid(); vih_iter++) {
		if (mesh.to_vertex_handle(vih_iter.handle()) == vh2 || mesh.from_vertex_handle(vih_iter.handle()) == vh2)
			return mesh.edge_handle(vih_iter.handle());
	}
	throw;
}
double MyViewer::getAngleSum(MyMesh::FaceHandle fh) {
	double angle = 0;
	for (MyMesh::FaceHalfedgeIter fh_iter = mesh.fh_iter(fh); fh_iter.is_valid(); fh_iter++) {
		MyMesh::Point vec1 = -(mesh.point(mesh.to_vertex_handle(fh_iter)) - mesh.point(mesh.from_vertex_handle(fh_iter))).normalized();
		MyMesh::HalfedgeHandle temph = mesh.next_halfedge_handle(fh_iter.handle());
		MyMesh::Point vec2 = (mesh.point(mesh.to_vertex_handle(temph)) - mesh.point(mesh.from_vertex_handle(temph))).normalized();

		angle += acos(dot(vec1, vec2));
	}
	return angle;
}
double MyViewer::getAngleSum(std::vector<MyMesh::VertexHandle> vertices) {
	double angle = 0;
	vertices.push_back(vertices[0]);
	vertices.push_back(vertices[1]);

	for (size_t i = 0; i < vertices.size() - 2; i++)
	{
		MyMesh::Point vec1 = -(mesh.point(vertices[i]) - mesh.point(vertices[i+1])).normalized();
		MyMesh::Point vec2 = (mesh.point(vertices[i + 1]) - mesh.point(vertices[i + 2])).normalized();
		angle += acos(dot(vec1, vec2));
	}

	return angle;
}
MyViewer::MyMesh::FaceHandle MyViewer::add_edge(MyMesh::FaceHandle fh, MyMesh::VertexHandle vh0, MyMesh::VertexHandle vh1) {
	std::vector<MyMesh::VertexHandle> fv1;
	std::vector<MyMesh::VertexHandle> fv2;

	MyMesh::HalfedgeHandle he;
	for (MyMesh::FaceHalfedgeIter fh_iter = mesh.fh_iter(fh); fh_iter.is_valid(); fh_iter++) {
		if (mesh.from_vertex_handle(fh_iter) == vh0)
			he = fh_iter.handle();
		//mesh.data(mesh.edge_handle(fh_iter.handle())).tagged = true;
	}

	MyMesh::VertexHandle vh = vh0;
	int idx = 0;
	while (vh != vh1) {
		fv1.push_back(vh);
		vh = mesh.to_vertex_handle(he);
		he = mesh.next_halfedge_handle(he);
	}
	fv1.push_back(vh);

	while (vh != vh0) {
		fv2.push_back(vh);
		vh = mesh.to_vertex_handle(he);
		he = mesh.next_halfedge_handle(he);
	}
	fv2.push_back(vh);
	mesh.delete_face(fh,false);
	MyMesh::FaceHandle newFace1 = mesh.add_face(fv1);
	MyMesh::FaceHandle newFace2 = mesh.add_face(fv2);

	if (fv1.size() == 4) {
		if (hasCommonEdge(fv1[0], fv1[2]))
			delete_edge(getCommonEdge(fv1[0], fv1[2]));
		if (hasCommonEdge(fv1[1], fv1[3]))
			delete_edge(getCommonEdge(fv1[1], fv1[3]));
		return newFace2;
	}
	if (fv2.size() == 4) {
		if (hasCommonEdge(fv2[0], fv2[2]))
			delete_edge(getCommonEdge(fv2[0], fv2[2]));
		if (hasCommonEdge(fv2[1], fv2[3]))
			delete_edge(getCommonEdge(fv2[1], fv2[3]));
		return newFace1;
	}
	return fh;
	
}
void MyViewer::resetFaceFlags() {
	for (auto f : mesh.faces()) {
		mesh.data(f).hasPrev = false;
		mesh.data(f).tagged = false;
		//mesh.data(f).tagged2 = false;
	}
}
int MyViewer::faceSides(MyMesh::FaceHandle fh) {
	int count = 0;
	for (MyMesh::FaceVertexIter fv_iter = mesh.fv_iter(fh); fv_iter.is_valid(); fv_iter++) {
		count++;
	}
	return count;
}
MyViewer::MyMesh::FaceHandle MyViewer::delete_edge(MyMesh::EdgeHandle _eh) {
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


	

	std::vector<MyMesh::HalfedgeHandle> temph;
	for (MyMesh::FaceHalfedgeIter fh_iter = mesh.fh_iter(f0); fh_iter.is_valid(); fh_iter++) {
		if (fh_iter.handle() != h0)
			temph.push_back(fh_iter.handle());
	}
	mesh.status(f0).set_deleted(true);
	mesh.status(h0).set_deleted(true);
	mesh.status(h1).set_deleted(true);
	mesh.status(_eh).set_deleted(true);
	//set face attributes
	mesh.set_halfedge_handle(f1, h1next);

	//set half_edge attributes
	for (auto h : temph) {
		mesh.set_face_handle(h, f1);
	}
	//mesh.set_face_handle(h0prev, f1);
	//mesh.set_face_handle(h0next, f1);
	//qDebug() << "ASD";
	mesh.set_next_halfedge_handle(h1prev, h0next);
	mesh.set_next_halfedge_handle(h0prev, h1next);

	mesh.set_prev_halfedge_handle(h1next, h0prev);
	mesh.set_prev_halfedge_handle(h0next, h1prev);

	//set vertex attributes
	mesh.set_halfedge_handle(v0, h1next);
	mesh.set_halfedge_handle(v1, h0next);

	return f1;
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