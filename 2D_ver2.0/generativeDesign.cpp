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
	findPointsInBoundary(20.0f);
	//Calculate the average edge length
	//double sum = 0.0;
	//int n = 0;
	////Calculate avg edge length
	//for (MyMesh::EdgeIter it = mesh.edges_begin(); it != mesh.edges_end(); ++it) {
	//	sum += mesh.calc_edge_length(*it);
	//	n++;
	//}
	//double L = (sum / n) * 2;
	//// Get an iterator over all halfedges
	//MyMesh::HalfedgeIter he_it, he_end = mesh.halfedges_end();
	//// If halfedge is boundary, lock the corresponding vertices
	//for (he_it = mesh.halfedges_begin(); he_it != he_end; ++he_it)
	//	if (mesh.is_boundary(*he_it)) {
	//		mesh.status(mesh.to_vertex_handle(*he_it)).set_locked(true);
	//		mesh.status(mesh.from_vertex_handle(*he_it)).set_locked(true);
	//	}
	//for (auto v : mesh.vertices()) {
	//	if (mesh.status(v).locked()) {
	//		int count = 0;
	//		for (MyMesh::VertexVertexIter vv_iter = mesh.vv_iter(v); vv_iter.is_valid(); vv_iter++) {

	//			if (mesh.status(*vv_iter).locked() && mesh.calc_edge_length(getCommonEdge(v, *vv_iter)) < L)
	//				count++;
	//		}
	//		if (count > 1)
	//			mesh.status(v).set_locked(false);
	//	}
	//}
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
			//if (it.handle().idx() == 1388 || it.handle().idx() == 1417) {
			//	mesh.data(it).tagged2 = true;
			//	continue;

			//}
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
bool MyViewer::hasCommonFace(MyMesh::VertexHandle vh1, MyMesh::VertexHandle vh2) {
	for (MyMesh::VertexFaceIter vf_iter = mesh.vf_iter(vh1); vf_iter.is_valid(); vf_iter++) {
		if (isVertexOnFace(vh2, vf_iter.handle()))
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

	MyMesh::FaceHandle returnHandlee = newFace1;
	if (fv1.size() == 4) {
		if (hasCommonEdge(fv1[0], fv1[2]))
			delete_edge(getCommonEdge(fv1[0], fv1[2]));
		if (hasCommonEdge(fv1[1], fv1[3]))
			delete_edge(getCommonEdge(fv1[1], fv1[3]));
		returnHandlee = newFace2;
	}
	if (fv2.size() == 4) {
		if (hasCommonEdge(fv2[0], fv2[2]))
			delete_edge(getCommonEdge(fv2[0], fv2[2]));
		if (hasCommonEdge(fv2[1], fv2[3]))
			delete_edge(getCommonEdge(fv2[1], fv2[3]));
		returnHandlee = newFace1;
	}
	return returnHandlee;
	
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

	mesh.adjust_outgoing_halfedge(v0);
	mesh.adjust_outgoing_halfedge(v1);


	return f1;
}
void MyViewer::printUnregularVertices() {
	for (auto f : mesh.faces()) {
		if (mesh.valence(f) != 4) {
			qDebug() << "BAJ VAN";
			//mesh.data(f).tagged2 = true;
		}
	}
	int countRegular = 0;
	int countUnregular = 0;
	for (auto v : mesh.vertices()) {
		if (mesh.valence(v) == 4) {
			countRegular++;
		}
		else
		{ if (!mesh.is_boundary(v))
			countUnregular++;
		}
	}
	qDebug() << "---Regular vertices: " << countRegular;
	qDebug() << "---Unregular vertices: " << countUnregular;

}
void MyViewer::quadRegularization(int iterations) {
	
	for (size_t i = 0; i < iterations; i++)
	{
		qDebug() << "First iteration";
		for (size_t i = 0; i < iterations; i++)
		{

			printUnregularVertices();
			qDebug() << "Remove2ValenceVertices";
			quadRegularizationRemove2ValenceVertices();
			printUnregularVertices();
			qDebug() << "quadRegularizationSwapping";
			quadRegularizationSwapping();
			printUnregularVertices();
			qDebug() << "quadRegularizationCollapsing";
			quadRegularizationCollapsing();
			printUnregularVertices();
			qDebug() << "quadRegularizationSplitting";
			quadRegularizationSplitting();
			printUnregularVertices();
			qDebug() << "quadRegularizationCompositons";
			quadRegularizationCompositons();
			printUnregularVertices();
			updateMesh();
		}
		update();
		qDebug() << "First transfer";
		for (size_t i = 0; i < iterations; i++)
		{
			quadRegularizationTransfer();
		}
		qDebug() << "Second iteration";
		for (size_t i = 0; i < iterations; i++)
		{

			printUnregularVertices();
			qDebug() << "Remove2ValenceVertices";
			quadRegularizationRemove2ValenceVertices();
			printUnregularVertices();
			qDebug() << "quadRegularizationSwapping";
			quadRegularizationSwapping();
			printUnregularVertices();
			qDebug() << "quadRegularizationCollapsing";
			quadRegularizationCollapsing();
			printUnregularVertices();
			qDebug() << "quadRegularizationSplitting";
			quadRegularizationSplitting();
			printUnregularVertices();
			qDebug() << "quadRegularizationCompositons";
			quadRegularizationCompositons();
			printUnregularVertices();
			updateMesh();
		}
		qDebug() << "Second transfer";
		for (size_t i = 0; i < iterations; i++)
		{
			quadRegularizationTransfer2();
		}
		
	}
	smoothQuadMesh();
	update();
}
void MyViewer::eliminate2and4ValenceBounradies() {
	for (auto v : mesh.vertices()) {
		qDebug() << "Test0";
		if (false &&mesh.is_boundary(v) && mesh.valence(v) == 4) {
			mesh.data(v).flags.tagged = true;
			MyMesh::VertexHandle helper1;
			MyMesh::VertexHandle helper2;
			MyMesh::FaceHandle helperFace;
			bool found1 = false;
			bool found2 = false;
			for (MyMesh::VertexVertexIter vv_iter = mesh.vv_iter(v); vv_iter.is_valid(); vv_iter++) {
				if (found1 && !mesh.is_boundary(*vv_iter)) {
					helper2 = vv_iter.handle();
					found2 = true;
					break;
				}
				if (!found1 && !mesh.is_boundary(*vv_iter)) {
					helper1 = vv_iter.handle();
					found1 = true;
				}
				
			}
			qDebug() << "Test1";
			if (found1 && found2) {
				found1 = false;
				for (MyMesh::VertexFaceIter vf_iter = mesh.vf_iter(v); vf_iter.is_valid(); vf_iter++) {
					if (isVertexOnFace(helper1, vf_iter.handle()) && isVertexOnFace(helper2, vf_iter.handle())) {
						found1 = true;
						helperFace = vf_iter.handle();
					}
				}
				if (found1) {
					qDebug() << "Test2";
					collapseVertices(helperFace,helper1,helper2);

					qDebug() << "Test3";
					mesh.garbage_collection();
					break;
				}

			}
			qDebug() << "Test4";
		}
		if (mesh.is_boundary(v) && mesh.valence(v) ==2) {
			MyMesh::VertexHandle helper1;
			MyMesh::VertexHandle helper2;
			MyMesh::FaceHandle helperFace;
			bool found1 = false;
			bool found2 = false;
			for (MyMesh::VertexVertexIter vv_iter = mesh.vv_iter(v); vv_iter.is_valid(); vv_iter++) {
				if (found1) {
					helper2 = vv_iter.handle();
					found2 = true;
					break;
				}
				if (!found1) {
					helper1 = vv_iter.handle();
					found1 = true;
				}

			}
			if (found1 && found2) {
				found1 = false;
				for (MyMesh::VertexFaceIter vf_iter = mesh.vf_iter(v); vf_iter.is_valid(); vf_iter++) {
					if (isVertexOnFace(helper1, vf_iter.handle()) && isVertexOnFace(helper2, vf_iter.handle())) {
						found1 = true;
						helperFace = vf_iter.handle();
					}
				}
				if (found1) {
					qDebug() << "Test2";
					collapseVertices(helperFace, helper1, helper2);

					qDebug() << "Test3";
					mesh.garbage_collection();
					break;
				}

			}

		}
	}
}
void MyViewer::quadRegularizationRemove2ValenceVertices() {
	int count = 1;

	int count2 = 0;
	while (count2 < count) {
		count = 0;
		count2 = 0;
		for (auto v : mesh.vertices())
			if (mesh.valence(v) == 2)
				count++;
		for (auto v : mesh.vertices()) {
			if (mesh.valence(v) == 2 && !mesh.is_boundary(v)) {

				MyMesh::HalfedgeHandle h = mesh.halfedge_handle(v);
				MyMesh::EdgeHandle e = mesh.edge_handle(mesh.next_halfedge_handle(h));


				MyMesh::HalfedgeHandle h2 = mesh.opposite_halfedge_handle(h);
				if (!mesh.is_boundary(e)) {
					mesh.data(e).tagged = true;
					MyMesh::VertexHandle v2 = mesh.from_vertex_handle(mesh.prev_halfedge_handle(mesh.opposite_halfedge_handle(mesh.next_halfedge_handle(h))));
					MyMesh::FaceHandle f = delete_edge(e);


					add_edge(f, v, v2);

				}
				else {
					MyMesh::EdgeHandle e = mesh.edge_handle(mesh.next_halfedge_handle(mesh.next_halfedge_handle(h)));
					if (!mesh.is_boundary(e)) {
						//mesh.data(e).tagged = true;
						MyMesh::VertexHandle v2 = mesh.to_vertex_handle(mesh.next_halfedge_handle(mesh.opposite_halfedge_handle(mesh.next_halfedge_handle(mesh.next_halfedge_handle(h)))));
						MyMesh::FaceHandle f = delete_edge(e);


						add_edge(f, v, v2);
					}
					else {
						MyMesh::EdgeHandle e = mesh.edge_handle(mesh.next_halfedge_handle(mesh.next_halfedge_handle(mesh.opposite_halfedge_handle(h))));
						if (!mesh.is_boundary(e)) {
							mesh.data(e).tagged = true;

							MyMesh::VertexHandle v2 = mesh.from_vertex_handle(mesh.prev_halfedge_handle(mesh.opposite_halfedge_handle(mesh.next_halfedge_handle(mesh.next_halfedge_handle(mesh.opposite_halfedge_handle(h))))));

							MyMesh::FaceHandle f = delete_edge(e);						

							add_edge(f, v, v2);

						}
						else {
							MyMesh::EdgeHandle e = mesh.edge_handle(mesh.prev_halfedge_handle(mesh.opposite_halfedge_handle(h)));
							if (!mesh.is_boundary(e)) {
								mesh.data(e).tagged = true;
								MyMesh::VertexHandle v2 = mesh.to_vertex_handle(mesh.next_halfedge_handle(mesh.opposite_halfedge_handle(mesh.prev_halfedge_handle(mesh.opposite_halfedge_handle(h)))));
								MyMesh::FaceHandle f = delete_edge(e);


								add_edge(f, v, v2);
							}
							else
								qDebug() << "Cant delete";
						}
					}

				}

			}
		}
		mesh.garbage_collection();
		for (auto v : mesh.vertices())
			if (mesh.valence(v) == 2)
				count2++;
	}


	for (auto v : mesh.vertices()) {
		//if (mesh.valence(v) == 2)
		//	mesh.data(v).flags.tagged = true;
		//else
			mesh.data(v).flags.tagged = false;
	}
}
void MyViewer::quadRegularizationSwapping() {
	for (MyMesh::EdgeIter it = mesh.edges_begin(); it != mesh.edges_end(); ++it) {
		MyMesh::VertexHandle v1 = mesh.from_vertex_handle(mesh.halfedge_handle(it, 0));
		MyMesh::VertexHandle v2 = mesh.from_vertex_handle(mesh.halfedge_handle(it, 1));

		if (mesh.valence(v1) == 5 && mesh.valence(v2) == 5) {
			MyMesh::VertexHandle v11 = mesh.to_vertex_handle(mesh.next_halfedge_handle(mesh.halfedge_handle(it, 0)));
			MyMesh::VertexHandle v12 = mesh.to_vertex_handle(mesh.next_halfedge_handle(mesh.halfedge_handle(it, 1)));

			MyMesh::VertexHandle v21 = mesh.from_vertex_handle(mesh.prev_halfedge_handle(mesh.halfedge_handle(it, 0)));
			MyMesh::VertexHandle v22 = mesh.from_vertex_handle(mesh.prev_halfedge_handle(mesh.halfedge_handle(it, 1)));

			if (mesh.valence(v11) == 3 && mesh.valence(v12) == 3) {
				MyMesh::FaceHandle tempf= delete_edge(it.handle());
				add_edge(tempf,v11, v12);
			}
			else if (mesh.valence(v21) == 3 && mesh.valence(v22) == 3) {
				MyMesh::FaceHandle tempf = delete_edge(it.handle());
				add_edge(tempf, v21, v22);
			}
			else if ((mesh.valence(v11) == 3 && mesh.valence(v12) == 4) || (mesh.valence(v11) == 4 && mesh.valence(v12) == 3)) {
				MyMesh::FaceHandle tempf = delete_edge(it.handle());
				add_edge(tempf, v11, v12);
			}
			else if ((mesh.valence(v21) == 3 && mesh.valence(v22) == 4) || (mesh.valence(v21) == 4 && mesh.valence(v22) == 3)) {
				MyMesh::FaceHandle tempf = delete_edge(it.handle());
				add_edge(tempf, v21, v22);
			}

		}
	}
	mesh.garbage_collection();
}
void MyViewer::collapseVertices(MyMesh::FaceHandle fh, MyMesh::VertexHandle vh1, MyMesh::VertexHandle vh2) {
	MyMesh::FaceHandle f = add_edge(fh, vh1, vh2);
	for (MyMesh::FaceHalfedgeIter fh_iter = mesh.fh_iter(fh); fh_iter.is_valid(); fh_iter++) {
		if ((mesh.from_vertex_handle(fh_iter) == vh1 && mesh.to_vertex_handle(fh_iter) == vh2) ||
			(mesh.from_vertex_handle(fh_iter) == vh2 && mesh.to_vertex_handle(fh_iter) == vh1)) {
			MyMesh::HalfedgeHandle he = fh_iter;
			if (mesh.is_boundary(mesh.from_vertex_handle(fh_iter))) {
				he = mesh.opposite_halfedge_handle(he);
			}
			mesh.collapse(he);
			//mesh.data(mesh.edge_handle(fh_iter)).tagged = true;

			break;
		}
	}

}
void MyViewer::quadRegularizationCollapsing() {
	bool collapse = true;
	while (collapse) {
		collapse = false;
		for (MyMesh::FaceIter it = mesh.faces_begin(); it != mesh.faces_end(); ++it) {
			std::vector<MyMesh::VertexHandle> vertices;
			for (MyMesh::FaceVertexIter fv_iter = mesh.fv_iter(it); fv_iter.is_valid(); fv_iter++)
			{
				vertices.push_back(fv_iter.handle());
			}
			if (mesh.valence(vertices[0]) == 5 && mesh.valence(vertices[2]) == 5 && mesh.valence(vertices[1]) == 3 && mesh.valence(vertices[3]) == 3) {
				MyMesh::FaceHandle f = add_edge(it.handle(), vertices[1], vertices[3]);
				for (MyMesh::FaceHalfedgeIter fh_iter = mesh.fh_iter(it); fh_iter.is_valid(); fh_iter++) {
					if ((mesh.from_vertex_handle(fh_iter) == vertices[1] && mesh.to_vertex_handle(fh_iter) == vertices[3]) ||
						(mesh.from_vertex_handle(fh_iter) == vertices[3] && mesh.to_vertex_handle(fh_iter) == vertices[1])) {
						MyMesh::HalfedgeHandle he = fh_iter;
						if (mesh.is_boundary(mesh.from_vertex_handle(fh_iter))) {
							he = mesh.opposite_halfedge_handle(he);
						}
						mesh.collapse(he);
						//mesh.data(mesh.edge_handle(fh_iter)).tagged = true;
						break;
					}
				}
				collapse = true;
				mesh.garbage_collection();
				break;
			}
			else if (mesh.valence(vertices[1]) == 5 && mesh.valence(vertices[3]) == 5 && mesh.valence(vertices[0]) == 3 && mesh.valence(vertices[2]) == 3) {
				collapseVertices(it.handle(), vertices[0], vertices[2]);				
				collapse = true;
				mesh.garbage_collection();
				break;
			}
			else if (mesh.valence(vertices[1]) == 5 && mesh.valence(vertices[3]) == 5 &&
				(mesh.valence(vertices[0]) == 3 && mesh.valence(vertices[2]) == 4 || mesh.valence(vertices[0]) == 4 && mesh.valence(vertices[2]) == 3)) {
				collapseVertices(it.handle(), vertices[0], vertices[2]);
				collapse = true;
				mesh.garbage_collection();
				break;
			}
			else if (mesh.valence(vertices[0]) == 5 && mesh.valence(vertices[2]) == 5 &&
				(mesh.valence(vertices[1]) == 3 && mesh.valence(vertices[3]) == 4 || mesh.valence(vertices[1]) == 4 && mesh.valence(vertices[3]) == 3)) {
				collapseVertices(it.handle(), vertices[1], vertices[3]);				
				collapse = true;
				mesh.garbage_collection();
				break;
			}
		}
	}
}
void MyViewer::quadRegularizationSplitting() {
	for (auto v : mesh.vertices()) {
		if (mesh.valence(v) > 4) {
			int lowValenceCount = 0;
			MyMesh::VertexHandle vl;
			MyMesh::VertexHandle vr;
			int edgeBetween = 0;
			for (MyMesh::VertexVertexIter vv_iter = mesh.vv_iter(v); vv_iter.is_valid(); vv_iter++) {
				if (mesh.valence(vv_iter) < 4) {
					if (lowValenceCount == 0) {
						vl = vv_iter.handle();
						edgeBetween = 0;
					}

					if (lowValenceCount == 1) {
						vr = vv_iter.handle();
						/*if (edgeBetween == 1) {
							MyMesh::VertexHandle temp = vr;
							vr = vl;
							vl = temp;
						}*/
					}
					lowValenceCount++;

				}
				edgeBetween++;
			}
			if (lowValenceCount >= 2) {
				MyMesh::Point p = mesh.point(v);
				mesh.set_point(v, (mesh.point(v) + mesh.point(vl) + mesh.point(vr)) / 3);
				MyMesh::HalfedgeHandle temph = mesh.vertex_split(p, v, vl, vr);
				delete_edge(mesh.edge_handle(temph));
				mesh.garbage_collection();

			
			}
				
		}
			
	}
}
void MyViewer::quadRegularizationCompositons() {
	qDebug() << "quadRegularizationCompositonSwapSplit";
	quadRegularizationCompositonSwapSplit();
	printUnregularVertices();
	qDebug() << "quadRegularizationCompositonSwapCollapse";
	quadRegularizationCompositonSwapCollapse();
	printUnregularVertices();
	qDebug() << "quadRegularizationCompositonSplitSplit";
	quadRegularizationCompositonSplitSplit();
	printUnregularVertices();
	qDebug() << "quadRegularizationCompositonCollapseCollapse";
	quadRegularizationCompositonCollapseCollapse();
	printUnregularVertices();
	qDebug() << "quadRegularizationCompositonSwapSwap";
	quadRegularizationCompositonSwapSwap();
	printUnregularVertices();
	qDebug() << "quadRegularizationCompositonSwapSwap2";
	quadRegularizationCompositonSwapSwap2();
	printUnregularVertices();
}
void MyViewer::quadRegularizationCompositonSwapSplit() {
	for (auto v : mesh.vertices()) {
		if (mesh.valence(v) == 5) {

			std::vector<MyMesh::VertexHandle> neighbours;
			std::vector<int> valences;
			for (MyMesh::VertexVertexIter vv_iter = mesh.vv_iter(v); vv_iter.is_valid(); vv_iter++) {
				neighbours.push_back(vv_iter.handle());
				valences.push_back(mesh.valence(vv_iter));
			}
			MyMesh::VertexHandle fourV;
			MyMesh::VertexHandle threeV;
			MyMesh::VertexHandle middle;
			bool found = true;
			if (valences[0] == 4 && valences[2] == 3) {
				fourV = neighbours[0];
				middle = neighbours[1];
				threeV = neighbours[2];
			}
			else 	if (valences[0] == 3 && valences[2] == 4) {
				fourV = neighbours[2];
				middle = neighbours[1];
				threeV = neighbours[0];
			}
			else if (valences[1] == 4 && valences[3] == 3) {
				fourV = neighbours[1];
				middle = neighbours[2];
				threeV = neighbours[3];
			}
			else if (valences[1] == 3 && valences[3] == 4) {
				fourV = neighbours[3];
				middle = neighbours[2];
				threeV = neighbours[1];
			}
			else if (valences[2] == 4 && valences[4] == 3) {
				fourV = neighbours[2];
				middle = neighbours[3];
				threeV = neighbours[4];
			}
			else if (valences[2] == 3 && valences[4] == 4) {
				fourV = neighbours[4];
				middle = neighbours[3];
				threeV = neighbours[2];
			}
			else if (valences[0] == 4 && valences[3] == 3) {
				fourV = neighbours[0];
				middle = neighbours[4];
				threeV = neighbours[3];
			}
			else if (valences[0] == 3 && valences[3] == 4) {
				fourV = neighbours[3];
				middle = neighbours[4];
				threeV = neighbours[0];
			}
			else if (valences[1] == 4 && valences[4] == 3) {
				fourV = neighbours[1];
				middle = neighbours[0];
				threeV = neighbours[4];
			}
			else if (valences[1] == 3 && valences[4] == 4) {
				fourV = neighbours[4];
				middle = neighbours[0];
				threeV = neighbours[1];
			}
			else
				found = false;
			if (found) {
				MyMesh::EdgeHandle e = getCommonEdge(fourV, v);
				MyMesh::FaceHandle f = mesh.face_handle(mesh.halfedge_handle(e, 0));
				if ( !isVertexOnFace(middle, f) ) {
					if (mesh.is_boundary(e))
						continue;
					f = mesh.face_handle(mesh.halfedge_handle(e, 1));
				}
				MyMesh::VertexHandle A;
				qDebug() << mesh.valence(f);
				for (MyMesh::FaceVertexIter fv_iter = mesh.fv_iter(f); fv_iter.is_valid(); fv_iter++) {
					if (fv_iter.handle() != middle && fv_iter.handle() != v && fv_iter.handle() != fourV) {
						A = fv_iter.handle();
						break;
					}
				}

				MyMesh::VertexHandle B;
				MyMesh::EdgeHandle e2 = getCommonEdge(fourV, A);
				MyMesh::FaceHandle f2 = mesh.face_handle(mesh.halfedge_handle(e2, 0));

				if (f2 == f) {
					if (mesh.is_boundary(e2))
						continue;
					f2 = mesh.face_handle(mesh.halfedge_handle(e2, 1));
				}

				for (MyMesh::VertexVertexIter vv_iter = mesh.vv_iter(A); vv_iter.is_valid(); vv_iter++) {
					if (isVertexOnFace(vv_iter.handle(), f2) && vv_iter.handle() != fourV) {
						B = vv_iter.handle();
						break;
					}
				}
				

				if ((mesh.valence(A) == 5 && mesh.valence(B) <= 4) || (mesh.valence(A) >= 4 && mesh.valence(B) == 3)) {
					MyMesh::FaceHandle temp = delete_edge(e2);
					add_edge(temp, B, v);
					MyMesh::Point p = mesh.point(v);
					mesh.set_point(v, (mesh.point(v) + mesh.point(fourV) + mesh.point(threeV)) / 3);
					MyMesh::HalfedgeHandle temph = mesh.vertex_split(p, v, fourV, threeV);
					delete_edge(mesh.edge_handle(temph));
					mesh.garbage_collection();

					break;
				}


			}


		}
	}
}
void MyViewer::quadRegularizationCompositonSwapCollapse() {
	int testcount = 0;
	for (auto v : mesh.vertices()) {
		if (mesh.valence(v) == 5) {
			
			quadRegularizationCompositonSwapCollapseHelper(v);

		}
	}
	qDebug() << "Muveletek: " << testcount;
	printUnregularVertices();
}
void MyViewer::quadRegularizationCompositonSwapCollapseHelper(MyMesh::VertexHandle v, bool otherNeighbour) {
	std::vector<MyMesh::VertexHandle> neighbours;
	std::vector<int> valences;
	for (MyMesh::VertexVertexIter vv_iter = mesh.vv_iter(v); vv_iter.is_valid(); vv_iter++) {
		neighbours.push_back(vv_iter.handle());
		valences.push_back(mesh.valence(vv_iter));
	}
	MyMesh::VertexHandle fourV;
	MyMesh::VertexHandle threeV;
	MyMesh::VertexHandle middle;
	bool found = true;
	if (valences[0] == 4 && valences[1] == 3 && !otherNeighbour) {
		fourV = neighbours[0];
		threeV = neighbours[1];
	}
	else if (valences[1] == 4 && valences[0] == 3 && !otherNeighbour) {
		fourV = neighbours[1];
		threeV = neighbours[0];
	}
	else if (valences[1] == 4 && valences[2] == 3 && !otherNeighbour) {
		fourV = neighbours[1];
		threeV = neighbours[2];
	}
	else if (valences[2] == 4 && valences[1] == 3) {
		fourV = neighbours[2];
		threeV = neighbours[1];
	}
	else if (valences[2] == 4 && valences[3] == 3 && !otherNeighbour) {
		fourV = neighbours[2];
		threeV = neighbours[3];
	}
	else if (valences[3] == 4 && valences[2] == 3) {
		fourV = neighbours[3];
		threeV = neighbours[2];
	}
	else if (valences[3] == 4 && valences[4] == 3) {
		fourV = neighbours[3];
		threeV = neighbours[4];
	}
	else if (valences[4] == 4 && valences[3] == 3) {
		fourV = neighbours[4];
		threeV = neighbours[3];
	}
	else if (valences[0] == 4 && valences[4] == 3) {
		fourV = neighbours[0];
		threeV = neighbours[4];
	}
	else if (valences[4] == 4 && valences[0] == 3) {
		fourV = neighbours[4];
		threeV = neighbours[0];
	}
	else {
		found = false;
	}
	if (!found)
		return;
	MyMesh::EdgeHandle e = getCommonEdge(v, fourV);
	if (mesh.is_boundary(e))
		return;
	MyMesh::FaceHandle f1 = mesh.face_handle(mesh.halfedge_handle(e, 0));
	MyMesh::FaceHandle f2 = mesh.face_handle(mesh.halfedge_handle(e, 1));
	MyMesh::VertexHandle fourV2;
	MyMesh::VertexHandle A;
	MyMesh::VertexHandle B;
	if (!isVertexOnFace(threeV, f1)) {
		MyMesh::FaceHandle temp = f1;
		f1 = f2;
		f2 = temp;
	}
	for (MyMesh::VertexVertexIter vv_iter = mesh.vv_iter(fourV); vv_iter.is_valid(); vv_iter++) {
		if (isVertexOnFace(vv_iter.handle(), f2) && vv_iter.handle() != v)
			fourV2 = vv_iter.handle();
		if (!isVertexOnFace(vv_iter.handle(), f2) && !isVertexOnFace(vv_iter.handle(), f1))
			A = vv_iter.handle();
	}

	MyMesh::EdgeHandle e2 = getCommonEdge(fourV, fourV2);
	MyMesh::FaceHandle f3 = mesh.face_handle(mesh.halfedge_handle(e2, 0));
	if (!isVertexOnFace(A, f3)) {
		if (mesh.is_boundary(e2))
			return;
		f3 = mesh.face_handle(mesh.halfedge_handle(e2, 1));

	}

	for (MyMesh::FaceVertexIter fv_iter = mesh.fv_iter(f3); fv_iter.is_valid(); fv_iter++) {
		if (fv_iter.handle() != A && fv_iter.handle() != fourV && fv_iter.handle() != fourV2)
			B = fv_iter.handle();
	}

	if ((mesh.valence(B) == 3 && mesh.valence(A) >= 4) || (mesh.valence(B) <= 4 && mesh.valence(A) >= 5)) {
		mesh.data(A).flags.tagged = true;
		mesh.data(fourV2).flags.tagged = true;

		MyMesh::FaceHandle temp = delete_edge(getCommonEdge(v, fourV));
		add_edge(temp, fourV2, threeV);
		mesh.garbage_collection();

		//collapse
		MyMesh::FaceHandle tempf = add_edge(f3, B, fourV);
		mesh.data(tempf).tagged2 = true;
		for (MyMesh::FaceHalfedgeIter fh_iter = mesh.fh_iter(tempf); fh_iter.is_valid(); fh_iter++) {
			if ((mesh.from_vertex_handle(fh_iter) == B && mesh.to_vertex_handle(fh_iter) == fourV) ||
				(mesh.from_vertex_handle(fh_iter) == fourV && mesh.to_vertex_handle(fh_iter) == B)) {
				MyMesh::HalfedgeHandle he = fh_iter.handle();
				if (mesh.is_boundary(mesh.from_vertex_handle(fh_iter))) {
					he = mesh.opposite_halfedge_handle(he);
				}

				mesh.data(mesh.edge_handle(he)).tagged = true;
				mesh.collapse(he);
				mesh.garbage_collection();
				//mesh.data(mesh.edge_handle(fh_iter)).tagged = true;
				break;
			}
		}
	}
	else if (!otherNeighbour)
		quadRegularizationCompositonSwapCollapseHelper(v, true);
}
void MyViewer::quadRegularizationCompositonSplitSplit() {
	for (auto v : mesh.vertices()) {
		if (mesh.valence(v) == 5) {
			std::vector<MyMesh::VertexHandle> neighbours;
			MyMesh::VertexHandle threeV;
			bool found = false;
			int idx = 0;
			int i = 0;
			for (MyMesh::VertexVertexIter vv_iter = mesh.vv_iter(v); vv_iter.is_valid(); vv_iter++) {
				if (mesh.valence(vv_iter.handle()) == 3) {
					threeV = vv_iter.handle();
					MyMesh::EdgeHandle eh = getCommonEdge(v, threeV);
					idx = i;
					found = true;
				}
				neighbours.push_back(vv_iter.handle());
				i++;
			}
			if (!found)
				continue;
			MyMesh::VertexHandle A1;
			MyMesh::VertexHandle A2;
			MyMesh::VertexHandle helper1;
			MyMesh::VertexHandle helper2;
			if (idx == 0) {
				A1 = neighbours[3];
				helper1 = neighbours[2];
				A2 = neighbours[2];
				helper2 = neighbours[3];
			} else 	if (idx == 1) {
				A1 = neighbours[3];
				helper1 = neighbours[4];
				A2 = neighbours[4];
				helper2 = neighbours[0];
			}
			else if (idx == 2) {
				A1 = neighbours[0];
				helper1 = neighbours[4];
				A2 = neighbours[4];
				helper2 = neighbours[0];
			}
			else if (idx == 3) {
				A1 = neighbours[0];
				helper1 = neighbours[1];
				A2 = neighbours[1];
				helper2 = neighbours[0];
			}
			else if (idx == 4) {
				A1 = neighbours[1];
				helper1 = neighbours[0];
				A2 = neighbours[2];
				helper2 = neighbours[1];
			}
			MyMesh::FaceHandle fHelper1;
			MyMesh::FaceHandle fHelper2;
			for (MyMesh::VertexFaceIter vf_iter = mesh.vf_iter(A1); vf_iter.is_valid(); vf_iter++) {
				if (isVertexOnFace(v, vf_iter.handle()) && isVertexOnFace(helper1, vf_iter.handle()))
					fHelper1 = vf_iter.handle();
			}
			for (MyMesh::VertexFaceIter vf_iter = mesh.vf_iter(A2); vf_iter.is_valid(); vf_iter++) {
				if (isVertexOnFace(v, vf_iter.handle()) && isVertexOnFace(helper2, vf_iter.handle()))
					fHelper2 = vf_iter.handle();
			}
			MyMesh::FaceHandle fHelper12;
			MyMesh::FaceHandle fHelper22;
			for (MyMesh::FaceFaceIter ff_iter = mesh.ff_iter(fHelper1); ff_iter.is_valid(); ff_iter++) {
				if (!isVertexOnFace(v, ff_iter.handle()) && !isVertexOnFace(helper1, ff_iter.handle())) {
					fHelper12 = ff_iter.handle();
					break;
				}
			}
			for (MyMesh::FaceFaceIter ff_iter = mesh.ff_iter(fHelper2); ff_iter.is_valid(); ff_iter++) {
				if (!isVertexOnFace(v, ff_iter.handle()) && !isVertexOnFace(helper2, ff_iter.handle())) {
					fHelper22 = ff_iter.handle();
					break;
				}
			}
			MyMesh::VertexHandle B1;
			MyMesh::VertexHandle B2;
			for (MyMesh::VertexVertexIter vv_iter = mesh.vv_iter(A1); vv_iter.is_valid(); vv_iter++) {
				if (isVertexOnFace(vv_iter.handle(), fHelper12) && !isVertexOnFace(vv_iter.handle(), fHelper1)) {
					B1 = vv_iter.handle();
					break;
				}
			}
			for (MyMesh::VertexVertexIter vv_iter = mesh.vv_iter(A2); vv_iter.is_valid(); vv_iter++) {
				if (isVertexOnFace(vv_iter.handle(), fHelper22) && !isVertexOnFace(vv_iter.handle(), fHelper1)) {
					B2 = vv_iter.handle();
					break;
				}
			}
			if ((mesh.valence(A1) == 5 && mesh.valence(B1) <= 4) || (mesh.valence(A1) == 4 && mesh.valence(B1) == 3)) {
				MyMesh::Point p = mesh.point(v);
				mesh.set_point(v, (mesh.point(v) + mesh.point(threeV) + mesh.point(A1)) / 3);
				MyMesh::HalfedgeHandle temph = mesh.vertex_split(p, v, threeV, A1);
				delete_edge(mesh.edge_handle(temph));
				mesh.garbage_collection();
			}
			else if ((mesh.valence(A2) == 5 && mesh.valence(B2) <= 4) || (mesh.valence(A2) == 4 && mesh.valence(B2) == 3)) {
				MyMesh::Point p = mesh.point(v);
				mesh.set_point(v, (mesh.point(v) + mesh.point(threeV) + mesh.point(A1)) / 3);
				MyMesh::HalfedgeHandle temph = mesh.vertex_split(p, v, A1, threeV);
				delete_edge(mesh.edge_handle(temph));
				mesh.garbage_collection();
			}
		}
	}
}
void MyViewer::quadRegularizationCompositonCollapseCollapse() {
	for (auto v : mesh.vertices()) {
		if (mesh.valence(v) == 5) {
			MyMesh::VertexHandle threeV;
			bool found = false;
			for (MyMesh::VertexVertexIter vv_iter = mesh.vv_iter(v); vv_iter.is_valid(); vv_iter++) {
				if (mesh.valence(vv_iter) == 3) {
					threeV = vv_iter.handle();
					found = true;
					break;
				}
			}
			if (!found)
				continue;
			found = false;
			MyMesh::VertexHandle fourV;
			MyMesh::VertexHandle fourV2;
			for (MyMesh::VertexVertexIter vv_iter = mesh.vv_iter(v); vv_iter.is_valid(); vv_iter++) {
				for (MyMesh::VertexVertexIter vv_iter2 = mesh.vv_iter(threeV); vv_iter2.is_valid(); vv_iter2++) {
					if (vv_iter2.handle() != v && mesh.valence(vv_iter2) == 4 && hasCommonEdge(vv_iter.handle(), vv_iter2) && mesh.valence(vv_iter.handle()) == 4) {
						fourV = vv_iter.handle();
						found = true;
						fourV2 = vv_iter2.handle();
						break;
					}
				}
				if (found)
					break;
			}
			if (!found)
				continue;

			found = false;
			MyMesh::VertexHandle A;
			MyMesh::VertexHandle B;
			for (MyMesh::VertexVertexIter vv_iter = mesh.vv_iter(fourV); vv_iter.is_valid(); vv_iter++) {
				for (MyMesh::VertexVertexIter vv_iter2 = mesh.vv_iter(fourV2); vv_iter2.is_valid(); vv_iter2++) {
					if (vv_iter.handle() != v && vv_iter.handle() != fourV2 && vv_iter2.handle() != fourV && vv_iter2.handle() != threeV &&
						hasCommonEdge(vv_iter.handle(), vv_iter2.handle())) {
						A = vv_iter.handle();
						B = vv_iter2.handle();
						found =true;
						break;
					}
				}
				if (found)
					break;
			}
			if (!found)
				continue;
			if ((mesh.valence(A) == 3 && mesh.valence(B) >= 4) || (mesh.valence(A) <= 4 && mesh.valence(B) >= 5)) {
				for (MyMesh::VertexFaceIter vf_iter = mesh.vf_iter(fourV); vf_iter.is_valid(); vf_iter++) {
					if (isVertexOnFace(v, vf_iter.handle()) && isVertexOnFace(threeV, vf_iter.handle())) {
						collapseVertices(vf_iter.handle(), fourV, threeV);
						break;
					}
				}
				for (MyMesh::VertexFaceIter vf_iter = mesh.vf_iter(A); vf_iter.is_valid(); vf_iter++) {
					if (isVertexOnFace(B, vf_iter.handle()) && isVertexOnFace(fourV2, vf_iter.handle())) {
						collapseVertices(vf_iter.handle(), A, fourV2);
						break;
					}
				}
				mesh.garbage_collection();
			}

		}
	}
}
void MyViewer::quadRegularizationCompositonSwapSwap() {
	for (auto f : mesh.faces()) {
		MyMesh::VertexHandle fiveV;
		MyMesh::VertexHandle threeV;
		MyMesh::VertexHandle fourV1;
		MyMesh::VertexHandle fourV2;
		std::vector<MyMesh::VertexHandle> vertices;
		std::vector<int> valences;
		for (MyMesh::FaceVertexIter fv_iter = mesh.fv_iter(f); fv_iter.is_valid(); fv_iter++) {
			vertices.push_back(fv_iter.handle());
			valences.push_back(mesh.valence(fv_iter.handle()));
		}
		bool found1 = false;
		bool found2 = false;

		if (valences[0] == 3 && valences[2] == 5) {
			fiveV = vertices[2];
			threeV = vertices[0];
			if (valences[1] == 4)
			{
				found1 = true;
				fourV1 = vertices[1];
			}
			if (valences[3] == 4)
			{
				found2 = true;
				fourV2 = vertices[3];
			}
		} else if (valences[2] == 3 && valences[0] == 5) {
			fiveV = vertices[0];
			threeV = vertices[2];
			if (valences[1] == 4)
			{
				found1 = true;
				fourV1 = vertices[1];
			}
			if (valences[3] == 4)
			{
				found2 = true;
				fourV2 = vertices[3];
			}
		}
		else if (valences[1] == 3 && valences[3] == 5) {
			fiveV = vertices[3];
			threeV = vertices[1];
			if (valences[0] == 4)
			{
				found1 = true;
				fourV1 = vertices[0];
			}
			if (valences[2] == 4)
			{
				found2 = true;
				fourV2 = vertices[2];
			}
		} else if (valences[3] == 3 && valences[1] == 5) {
			fiveV = vertices[1];
			threeV = vertices[3];
			if (valences[0] == 4)
			{
				found1 = true;
				fourV1 = vertices[0];
			}
			if (valences[2] == 4)
			{
				found2 = true;
				fourV2 = vertices[2];
			}
		}
		MyMesh::FaceHandle f1;
		MyMesh::FaceHandle f2;
		if (found1 && mesh.is_boundary(getCommonEdge(fiveV, fourV1)))
			found1 = false;
		if (found2 && mesh.is_boundary(getCommonEdge(fiveV, fourV2)))
			found2 = false;
		for (MyMesh::FaceFaceIter ff_iter = mesh.ff_iter(f); ff_iter.is_valid(); ff_iter++) {
			if (found1 && isVertexOnFace(fiveV, ff_iter.handle()) && isVertexOnFace(fourV1, ff_iter.handle())) {
				f1 = ff_iter.handle();
			}
			if (found2 && isVertexOnFace(fiveV, ff_iter.handle()) && isVertexOnFace(fourV2, ff_iter.handle())) {
				f2 = ff_iter.handle();
			}
		}
		if (found1) {
			MyMesh::VertexHandle fourV3;
			MyMesh::VertexHandle A;
			bool found12 = false;
			for (MyMesh::FaceVertexIter fv_iter = mesh.fv_iter(f1); fv_iter.is_valid(); fv_iter++) {
				if (fv_iter.handle() != fourV1 && fv_iter.handle() != fiveV && hasCommonEdge(fiveV, fv_iter.handle()) && mesh.valence(fv_iter.handle()) == 4) {

					fourV3 = fv_iter.handle();
					found12 = true;
				}
				if (fv_iter.handle() != fourV1 && fv_iter.handle() != fiveV && hasCommonEdge(fourV1, fv_iter.handle()))
					A = fv_iter.handle();
			}
			if (found12) {
				MyMesh::VertexHandle B;
				found12 = false;
				for (MyMesh::VertexVertexIter vv_iter = mesh.vv_iter(fourV3); vv_iter.is_valid(); vv_iter++) {
					if (vv_iter.handle() != fiveV && hasCommonFace(A, vv_iter.handle())) {
						B = vv_iter.handle();
						found12 = true;
						break;
					}
				}
				if (found12 && ((mesh.valence(A) == 5 && mesh.valence(B) <= 4) ||(mesh.valence(A) >=4 && mesh.valence(B) == 3) )) {
					MyMesh::EdgeHandle e = getCommonEdge(fiveV, fourV1);
					MyMesh::FaceHandle f = delete_edge(e);
					add_edge(f, threeV, fourV3);

					e = getCommonEdge(fourV3, A);
					f = delete_edge(e);
					add_edge(f, fourV1, B);
					mesh.garbage_collection();
					

				}
			}
		}

		else if (found2) {
			MyMesh::VertexHandle fourV3;
			MyMesh::VertexHandle A;
			bool found22 = false;
			for (MyMesh::FaceVertexIter fv_iter = mesh.fv_iter(f2); fv_iter.is_valid(); fv_iter++) {
				if (fv_iter.handle() != fourV2 && fv_iter.handle() != fiveV && hasCommonEdge(fiveV, fv_iter.handle()) && mesh.valence(fv_iter.handle()) == 4) {

					fourV3 = fv_iter.handle();
					found22 = true;
				}
				if (fv_iter.handle() != fourV2 && fv_iter.handle() != fiveV && hasCommonEdge(fourV2, fv_iter.handle()))
					A = fv_iter.handle();
			}
			if (found22) {
				MyMesh::VertexHandle B;
				found22 = false;
				for (MyMesh::VertexVertexIter vv_iter = mesh.vv_iter(fourV3); vv_iter.is_valid(); vv_iter++) {
					if (vv_iter.handle() != fiveV && hasCommonFace(A, vv_iter.handle())) {
						B = vv_iter.handle();
						found22 = true;
						break;
					}
				}
				if (found22 && ((mesh.valence(A) == 5 && mesh.valence(B) <= 4) || (mesh.valence(A) >= 4 && mesh.valence(B) == 3))) {
					MyMesh::EdgeHandle e = getCommonEdge(fiveV, fourV2);
					MyMesh::FaceHandle f = delete_edge(e);
					add_edge(f, threeV, fourV3);

					e = getCommonEdge(fourV3, A);
					f = delete_edge(e);
					add_edge(f, fourV2, B);

					mesh.garbage_collection();
					
				}
			}
		}
	}
}
void MyViewer::quadRegularizationCompositonSwapSwap2() {
	for (auto v : mesh.vertices()) {
		if (mesh.valence(v) == 5) {
			std::vector<MyMesh::VertexHandle> neighbours;
			std::vector<int> valences;
			for (MyMesh::VertexVertexIter vv_iter = mesh.vv_iter(v); vv_iter.is_valid(); vv_iter++) {
				valences.push_back(mesh.valence(vv_iter.handle()));
				neighbours.push_back(vv_iter.handle());
			}
			MyMesh::VertexHandle threeV;
			MyMesh::VertexHandle fourV1;
			MyMesh::VertexHandle fourV2;
			bool found1 = false;
			bool found2 = false;
			if (valences[0] == 3 && !mesh.is_boundary(neighbours[0])) {
				threeV = neighbours[0];
				if (valences[1] == 4) {
					fourV1 = neighbours[1];
					found1 = true;
				}
				if (valences[4] == 4) {
					fourV2 = neighbours[4];
					found2 = true;
				}
			} else 	if (valences[1] == 3 && !mesh.is_boundary(neighbours[1])) {
				threeV = neighbours[1];
				if (valences[0] == 4) {
					fourV1 = neighbours[0];
					found1 = true;
				}
				if (valences[2] == 4) {
					fourV2 = neighbours[2];
					found2 = true;
				}
			} else if (valences[2] == 3 && !mesh.is_boundary(neighbours[2])) {
				threeV = neighbours[2];
				if (valences[1] == 4) {
					fourV1 = neighbours[1];
					found1 = true;
				}
				if (valences[3] == 4) {
					fourV2 = neighbours[3];
					found2 = true;
				}
			} else if (valences[3] == 3 && !mesh.is_boundary(neighbours[3])) {
				threeV = neighbours[3];
				if (valences[2] == 4) {
					fourV1 = neighbours[2];
					found1 = true;
				}
				if (valences[4] == 4 ) {
					fourV2 = neighbours[4];
					found2 = true;
				}
			} else if (valences[4] == 3 && !mesh.is_boundary(neighbours[4])) {
				threeV = neighbours[4];
				if (valences[3] == 4) {
					fourV1 = neighbours[3];
					found1 = true;
				}
				if (valences[0] == 4) {
					fourV2 = neighbours[0];
					found2 = true;
				}
			}
			if (found1) {
				MyMesh::VertexHandle fourV3;
				bool found3 = false;
				for (MyMesh::VertexVertexIter vv_iter = mesh.vv_iter(fourV1); vv_iter.is_valid(); vv_iter++) {
					if (vv_iter.handle() != v && hasCommonFace(vv_iter.handle(), v) && !hasCommonFace(vv_iter.handle(), threeV) && mesh.valence(vv_iter.handle()) == 4) {
						found3 = true;
						fourV3 = vv_iter.handle();
						break;
					}
				}
				if (found3) {
					MyMesh::VertexHandle A;
					found3 = false;
					for (MyMesh::VertexVertexIter vv_iter = mesh.vv_iter(fourV3); vv_iter.is_valid(); vv_iter++) {
						if (vv_iter.handle() != fourV1 && hasCommonFace(vv_iter.handle(), fourV1) && !hasCommonFace(vv_iter.handle(), v)) {
							found3 = true;
							A = vv_iter.handle();
							break;
						}
					}
					if (found3) {
						MyMesh::VertexHandle B;
						found3 = false;
						for (MyMesh::VertexVertexIter vv_iter = mesh.vv_iter(A); vv_iter.is_valid(); vv_iter++) {
							if (vv_iter.handle() != fourV3 && hasCommonFace(vv_iter.handle(), fourV3) && !hasCommonFace(vv_iter.handle(), fourV1)) {
								found3 = true;
								B = vv_iter.handle();
								break;
							}
						}
						if (found3 && ((mesh.valence(A) >= 5 && mesh.valence(B) <= 4) || (mesh.valence(A) >= 4 && mesh.valence(B) == 3))) {
							MyMesh::EdgeHandle e = getCommonEdge(v, fourV1);
							MyMesh::FaceHandle f = delete_edge(e);
							add_edge(f, threeV, fourV3);

							e = getCommonEdge(fourV3, A);
							f = delete_edge(e);
							add_edge(f, fourV1, B);
							mesh.garbage_collection();
						}

					}
				}
			}
			if (found2) {
				MyMesh::VertexHandle fourV3;
				bool found3 = false;
				for (MyMesh::VertexVertexIter vv_iter = mesh.vv_iter(fourV2); vv_iter.is_valid(); vv_iter++) {
					if (vv_iter.handle() != v && hasCommonFace(vv_iter.handle(), v) && !hasCommonFace(vv_iter.handle(), threeV) && mesh.valence(vv_iter.handle()) == 4) {
						found3 = true;
						fourV3 = vv_iter.handle();
						break;
					}
				}
				if (found3) {
					MyMesh::VertexHandle A;
					found3 = false;
					for (MyMesh::VertexVertexIter vv_iter = mesh.vv_iter(fourV3); vv_iter.is_valid(); vv_iter++) {
						if (vv_iter.handle() != fourV2 && hasCommonFace(vv_iter.handle(), fourV2) && !hasCommonFace(vv_iter.handle(), v)) {
							found3 = true;
							A = vv_iter.handle();
							break;
						}
					}
					if (found3) {
						MyMesh::VertexHandle B;
						found3 = false;
						for (MyMesh::VertexVertexIter vv_iter = mesh.vv_iter(A); vv_iter.is_valid(); vv_iter++) {
							if (vv_iter.handle() != fourV3 && hasCommonFace(vv_iter.handle(), fourV3) && !hasCommonFace(vv_iter.handle(), fourV2)) {
								found3 = true;
								B = vv_iter.handle();
								break;
							}
						}
						if (found3 && ((mesh.valence(A) >= 5 && mesh.valence(B) <= 4) || (mesh.valence(A) >= 4 && mesh.valence(B) == 3))) {
							MyMesh::EdgeHandle e = getCommonEdge(v, fourV2);
							MyMesh::FaceHandle f = delete_edge(e);
							add_edge(f, threeV, fourV3);

							e = getCommonEdge(fourV3, A);
							f = delete_edge(e);
							add_edge(f, fourV2, B);
							mesh.garbage_collection();
						}
					}
				}
			}
		}
	}
}
void MyViewer::quadRegularizationTransfer() {
	for (auto v : mesh.vertices()) {
		if (mesh.valence(v) >= 5) {
			bool found = false;
			MyMesh::VertexHandle threeV;
			for (MyMesh::VertexFaceIter vf_iter = mesh.vf_iter(v); vf_iter.is_valid(); vf_iter++) {
				found = false;
				for (MyMesh::FaceVertexIter fv_iter = mesh.fv_iter(vf_iter); fv_iter.is_valid(); fv_iter++) {
					if (!mesh.is_boundary(fv_iter.handle()) && mesh.valence(fv_iter.handle()) == 3 && fv_iter.handle() != v && !hasCommonEdge(fv_iter,v))
					{
						found = true;
						threeV = fv_iter.handle();
						break;
					}
				}
				if (found)
					break;
			}
			if (found) {
				MyMesh::VertexHandle helper;
				MyMesh::VertexHandle helper2;
				found = false;
				bool found2 = false;
				for (MyMesh::VertexVertexIter vv_iter = mesh.vv_iter(threeV); vv_iter.is_valid(); vv_iter++) {
					if (hasCommonFace(vv_iter.handle(), v))
					{
						helper = vv_iter.handle();
						found = true;
						break;
					}
				}
				for (MyMesh::VertexVertexIter vv_iter = mesh.vv_iter(v); vv_iter.is_valid(); vv_iter++) {
					if (hasCommonFace(vv_iter.handle(), helper) && !hasCommonFace(vv_iter.handle(), threeV))
					{
						helper2 = vv_iter.handle();
						found2 = false;
						break;
					}
				}
				if (!found || !found2)
					continue;
				if (mesh.is_boundary(helper2) || (mesh.is_boundary(helper2) && mesh.valence(helper) == 3))
					continue;
				MyMesh::EdgeHandle e = getCommonEdge(helper, v);
				MyMesh::FaceHandle temp = delete_edge(e);
				add_edge(temp, helper2, threeV);
				mesh.data(helper).flags.tagged = true;
				mesh.data(v).flags.tagged = true;
				mesh.garbage_collection();
			}
		}
	}
}
void MyViewer::quadRegularizationTransfer2() {
	for (auto v : mesh.vertices()) {
		//if (mesh.valence(v) == 4 || mesh.is_boundary(v))
		//	mesh.data(v).flags.tagged = false;
		//else
		//	mesh.data(v).flags.tagged = true;

		if (mesh.valence(v) >= 5) {
			MyMesh::VertexHandle threeV;
			bool found = false;
			for (MyMesh::VertexVertexIter vv_iter = mesh.vv_iter(v); vv_iter.is_valid(); vv_iter++) {
				if (mesh.valence(vv_iter.handle()) == 3 && !mesh.is_boundary(vv_iter.handle()))
				{
					threeV = vv_iter.handle();
					//mesh.data(vv_iter.handle()).flags.tagged = true;
					//mesh.data(v).flags.tagged = true;
					found = true;
					break;
				}
			}
			if (!found)
				continue;
			qDebug() << "TEST0";

			MyMesh::VertexHandle helper;
			MyMesh::VertexHandle helper2;
			found = false;
			bool found2 = false;
			for (MyMesh::VertexVertexIter vv_iter = mesh.vv_iter(v); vv_iter.is_valid(); vv_iter++) {
				if (vv_iter.handle() != threeV && hasCommonFace(vv_iter.handle(), threeV))
				{
					helper = vv_iter.handle();
					found = true;
					break;
				}
			}
			for (MyMesh::VertexVertexIter vv_iter = mesh.vv_iter(helper); vv_iter.is_valid(); vv_iter++) {
				if (vv_iter.handle() != v && hasCommonFace(vv_iter.handle(), v) &&!hasCommonFace(vv_iter.handle(), threeV))
				{
					helper2 = vv_iter.handle();
					found2 = true;
					break;
				}
			}
			qDebug() << "TEST1";
			if (!found || !found2)
				continue;
			MyMesh::EdgeHandle e = getCommonEdge(helper, v);
			MyMesh::FaceHandle temp = delete_edge(e);
			add_edge(temp, helper2, threeV);
			mesh.garbage_collection();
			qDebug() << "TEST2";


		}
	}
}
MyViewer::MyMesh::VertexHandle MyViewer::findWaytoNearestUnregular(MyMesh::VertexHandle vh, MyMesh::VertexHandle notThis) {
	for (auto v : mesh.vertices()) {
		mesh.data(v).Bfs.hasPrev = false;
	}
	///BFS - START
	/// ############
	std::queue<MyMesh::VertexHandle> prevVertices;
	bool foundOtherVertex = false;
	MyMesh::VertexHandle pairVertex;
	prevVertices.push(vh);
	mesh.data(vh).Bfs.hasPrev = true;
	MyMesh::VertexHandle tempVh;
	while (prevVertices.size() > 0 && !foundOtherVertex) {
		tempVh = prevVertices.front();
		prevVertices.pop();

		for (MyMesh::VertexVertexIter vv_iter = mesh.vv_iter(tempVh); vv_iter.is_valid(); vv_iter++) {
			if (!mesh.data(vv_iter).Bfs.hasPrev) {

				prevVertices.push(vv_iter.handle());
				mesh.data(vv_iter).Bfs.hasPrev = true;
				mesh.data(vv_iter).Bfs.prevVertex = tempVh;

				if (mesh.valence(vv_iter.handle()) != 4 && !mesh.is_boundary(vv_iter) && notThis != vv_iter.handle()) {
					foundOtherVertex = true;
					pairVertex = vv_iter.handle();
					break;
				}
			}
		}

	}
	while (pairVertex != vh) {
		if (mesh.data(pairVertex).Bfs.prevVertex != vh)
			pairVertex = mesh.data(pairVertex).Bfs.prevVertex;
		else
			return pairVertex;
	}
	return pairVertex;
	///BFS - END
	/// ############
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

void MyViewer::smoothQuadMesh(int iterations) {
	for (size_t i = 0; i < iterations; i++)
	{
		for (auto v : mesh.vertices()) {
			if (mesh.is_boundary(v)) {
			//	mesh.data(v).flags.tagged = true;
				mesh.data(v).newPos = mesh.point(v);
			}
			//else
			//	mesh.data(v).flags.tagged = false;

		}
		OpenMesh::Smoother::JacobiLaplaceSmootherT<MyMesh> smoother(mesh);
		smoother.initialize(OpenMesh::Smoother::SmootherT<MyMesh>::Tangential_and_Normal, // or: Tangential_and_Normal
			OpenMesh::Smoother::SmootherT<MyMesh>::C1);
		for (size_t i = 1; i <= 10; ++i) {
			smoother.smooth(10);
		}

		updateMesh(false);
		for (auto v : mesh.vertices()) {
			if (mesh.is_boundary(v))
				mesh.set_point(v, mesh.data(v).newPos);
			mesh.data(v).flags.tagged = true;
		}
		reMeshVertexPositions();
		for (auto v : mesh.vertices()) {
			mesh.data(v).flags.tagged = false;
		}
		updateMesh(false);
	}
}


//Ver2.0

void MyViewer::createQuadPartitioning(double maxLength) {
	int taggedNum = findPointsInBoundary(maxLength);
	createQuadsFromTagged(taggedNum);
}
int MyViewer::findPointsInBoundary(double maxLength) {
	mesh.reset_status();
	double maxAngle = -1.0f;
	MyMesh::VertexHandle maxangled = mesh.vertices_begin();
	std::vector<MyMesh::VertexHandle> orderedBoundaryVertices;
	for (auto v : mesh.vertices()) {
		mesh.data(v).flags.tagged = false;
		mesh.data(v).flags.temporary_tagged = false;
		if (mesh.is_boundary(v)) {
			MyMesh::VertexHandle neighbour0;
			MyMesh::VertexHandle neighbour1;
			bool found1 = false;
			for (MyMesh::VVIter vv_iter = mesh.vv_iter(v); vv_iter.is_valid(); vv_iter++) {
				if (mesh.is_boundary(vv_iter) && mesh.is_boundary(getCommonEdge(v,vv_iter.handle()))) {
					if (!found1)
					{
						neighbour0 = vv_iter.handle();
						found1 = true;
					}
					else {
						neighbour1 = vv_iter.handle();
						break;
					}
				}
			}
			OpenMesh::Vec3d vNeighbour0 = (mesh.point(v) - mesh.point(neighbour0)).normalized();
			OpenMesh::Vec3d vNeighbour1 = (mesh.point(v) - mesh.point(neighbour1)).normalized();

			double dotProduct = dot(vNeighbour0, vNeighbour1);
			mesh.data(v).angle = dotProduct;
			/*if (dotProduct > -0.75f) {
				mesh.data(v).flags.tagged = true;
			}*/
			//if (dotProduct > maxAngle) {
			//	maxAngle = dotProduct;
			//	mesh.data(maxangled).flags.tagged = false;
			//	mesh.data(v).flags.tagged = true;
			//	maxangled = v;
			//}
			orderedBoundaryVertices.push_back(v);
		}
	}
	MyViewer::Comparer comparer = MyViewer::Comparer(mesh);
	std::sort(orderedBoundaryVertices.begin(), orderedBoundaryVertices.end(), comparer);
	int counter = 0;
	std::vector<MyMesh::VertexHandle> taggedV;
	mesh.data(*orderedBoundaryVertices.begin()).flags.tagged = true;
	taggedV.push_back(*orderedBoundaryVertices.begin());

	for (auto v : orderedBoundaryVertices)
	{
		bool good = true;
		/*for (auto v2 : taggedV)
		{
			if ((mesh.point(v2) - mesh.point(v)).length() < 0.8f * maxLength) {
				good = false;
				break;
			}
		}*/

		MyMesh::EdgeHandle e1;
		MyMesh::EdgeHandle e2;
		bool found1 = false;
		for (MyMesh::VEIter ve_iter = mesh.ve_iter(v); ve_iter.is_valid(); ve_iter++) {
			if (mesh.is_boundary(ve_iter)) {
				if (found1) {
					e2 = ve_iter.handle();
				}
				else {
					e1 = ve_iter.handle();
					found1 = true;
				}
			}
		}
		MyMesh::VertexHandle currentHandle = v;
		double sumLength = 0;
		while (!mesh.data(currentHandle).flags.tagged) {
			MyMesh::HalfedgeHandle h = mesh.halfedge_handle(e1, 0);
			if (mesh.from_vertex_handle(h) != currentHandle) {
				h = mesh.halfedge_handle(e1, 1);
			}
			currentHandle = mesh.to_vertex_handle(h);
			sumLength += (mesh.point(mesh.to_vertex_handle(h)) - mesh.point(mesh.from_vertex_handle(h))).length();
			for (MyMesh::VEIter ve_iter = mesh.ve_iter(currentHandle); ve_iter.is_valid(); ve_iter++) {
				if (mesh.is_boundary(ve_iter) && ve_iter.handle() != e1) {
					e1 = ve_iter.handle();
					break;
				}
			}
		}
		double sumLength2 = 0;
		currentHandle = v;
		while (!mesh.data(currentHandle).flags.tagged) {
			MyMesh::HalfedgeHandle h = mesh.halfedge_handle(e2, 0);
			if (mesh.from_vertex_handle(h) != currentHandle) {
				h = mesh.halfedge_handle(e2, 1);
			}
			currentHandle = mesh.to_vertex_handle(h);
			sumLength2 += (mesh.point(mesh.to_vertex_handle(h)) - mesh.point(mesh.from_vertex_handle(h))).length();
			for (MyMesh::VEIter ve_iter = mesh.ve_iter(currentHandle); ve_iter.is_valid(); ve_iter++) {
				if (mesh.is_boundary(ve_iter) && ve_iter.handle() != e2) {
					e2 = ve_iter.handle();
					break;
				}
			}
		}
		if (sumLength < 0.8f * maxLength || sumLength2 < 0.8f * maxLength)
			good = false;
		//qDebug() << sumLength << "\n";
		//qDebug() << sumLength2 <<"\n";

		if (good) {
			mesh.data(v).flags.tagged = true;
			mesh.status(v).set_locked(true);
			taggedV.push_back(v);
		}
	}
	if (taggedV.size() % 2 != 0) {
		return findPointsInBoundary(maxLength * 0.9f);
	} else
		return taggedV.size();
}

bool MyViewer::Comparer::operator()(const MyMesh::VertexHandle& v1, const MyMesh::VertexHandle& v2) {
	return (mesh.data(v1).angle > mesh.data(v2).angle);
}
void MyViewer::createQuadsFromTagged(int taggedNum) {
	std::vector<MyMesh::VertexHandle> tagged;
	for (auto v: mesh.vertices())
	{
		if (mesh.data(v).flags.tagged == true)
		{
			tagged.push_back(v);
			mesh.data(v).flags.temporary_tagged = true;
			break;
		}
	}
	while (tagged.size() < taggedNum) {
		MyMesh::EdgeHandle e1;
		for (MyMesh::VEIter ve_iter = mesh.ve_iter(tagged[tagged.size()-1]); ve_iter.is_valid(); ve_iter++) {
			MyMesh::HalfedgeHandle h = mesh.halfedge_handle(ve_iter,0);
			if (mesh.is_boundary(ve_iter) &&  !(mesh.data(mesh.from_vertex_handle(h)).flags.temporary_tagged && mesh.data(mesh.to_vertex_handle(h)).flags.temporary_tagged)){
					e1 = ve_iter.handle();
					break;
			}
		}
		MyMesh::VertexHandle currentHandle = tagged[tagged.size() - 1];
		bool first = true;
		while (!mesh.data(currentHandle).flags.tagged || first) {
			MyMesh::HalfedgeHandle h = mesh.halfedge_handle(e1, 0);
			if (mesh.from_vertex_handle(h) != currentHandle) {
				h = mesh.halfedge_handle(e1, 1);
			}
			currentHandle = mesh.to_vertex_handle(h);
			for (MyMesh::VEIter ve_iter = mesh.ve_iter(currentHandle); ve_iter.is_valid(); ve_iter++) {
				if (mesh.is_boundary(ve_iter) && ve_iter.handle() != e1) {
					e1 = ve_iter.handle();
					break;
				}
			}
			first = false;
		}
		tagged.push_back(currentHandle);
		mesh.data(currentHandle).flags.temporary_tagged = true;
	}
	//int count = 3;
	while (tagged.size() > 0 /*&& count > 0*/) {
		qDebug() << "#######################xxx\n";
		int bestIdx = 0;
		double minDot = 100000;
		for (size_t i = 0; i < tagged.size(); i++)
		{
			int idx0 = i;
			int idx1 = (i + 1) % tagged.size();
			int idx2 = (i + 2) % tagged.size();
			int idx3 = (i + 3) % tagged.size();

			MyMesh::VertexHandle v0 = tagged[idx0];
			MyMesh::VertexHandle p0 = tagged[idx1];
			MyMesh::VertexHandle v1 = tagged[idx2];
			MyMesh::VertexHandle p1 = tagged[idx3];

			OpenMesh::Vec3d p0v0 = (mesh.point(v0) - mesh.point(p0)).normalized();
			OpenMesh::Vec3d p0v1 = (mesh.point(v1) - mesh.point(p0)).normalized();
			OpenMesh::Vec3d p1v0 = (mesh.point(v0) - mesh.point(p1)).normalized();
			OpenMesh::Vec3d p1v1 = (mesh.point(v1) - mesh.point(p1)).normalized();

			OpenMesh::Vec3d v0p0 = (mesh.point(p0) - mesh.point(v0)).normalized();
			OpenMesh::Vec3d v0p1 = (mesh.point(p1) - mesh.point(v0)).normalized();
			OpenMesh::Vec3d v1p0 = (mesh.point(p0) - mesh.point(v1)).normalized();
			OpenMesh::Vec3d v1p1 = (mesh.point(p1) - mesh.point(v1)).normalized();
			bool convex = true;
			if (isLeft(mesh.point(v0), mesh.point(p0), mesh.point(v1)) != isLeft(mesh.point(v0), mesh.point(p0), mesh.point(p1)))
				convex = false;
			if (isLeft(mesh.point(p0), mesh.point(v1), mesh.point(v0)) != isLeft(mesh.point(p0), mesh.point(v1), mesh.point(p1)))
				convex = false;
			if (isLeft(mesh.point(v1), mesh.point(p1), mesh.point(v0)) != isLeft(mesh.point(v1), mesh.point(p1), mesh.point(p0)))
				convex = false;
			if (isLeft(mesh.point(p1), mesh.point(v0), mesh.point(v1)) != isLeft(mesh.point(p1), mesh.point(v0), mesh.point(p0)))
				convex = false;
			if (!convex)
				continue;
			double temp = abs(dot(p0v0, p0v1)) + abs(dot(p1v0, p1v1)) + abs(dot(v0p0, v0p1)) + abs(dot(v1p0, v1p1));
			qDebug() << "idx: " << idx0 << " dot: " << temp << "\n";
			if (temp < minDot) {
				bestIdx = idx0;
				minDot = temp;
			}
		}
		qDebug() << "bestidx: " << bestIdx << " dot: " << minDot << "\n";
		int idx0 = bestIdx;
		int idx1 = (bestIdx + 1) % tagged.size();
		int idx2 = (bestIdx + 2) % tagged.size();
		int idx3 = (bestIdx + 3) % tagged.size();

		MyMesh::VertexHandle vhandle[4];
		MyMesh::Point offset = MyMesh::Point(0, 0, 0.3);
		vhandle[0] = quadPartition.add_vertex(mesh.point(tagged[idx0]) + offset);
		vhandle[1] = quadPartition.add_vertex(mesh.point(tagged[idx1]) + offset);
		vhandle[2] = quadPartition.add_vertex(mesh.point(tagged[idx2]) + offset);
		vhandle[3] = quadPartition.add_vertex(mesh.point(tagged[idx3]) + offset);
		std::vector<MyMesh::VertexHandle>  face_vhandles;
		face_vhandles.clear();
		face_vhandles.push_back(vhandle[0]);
		face_vhandles.push_back(vhandle[1]);
		face_vhandles.push_back(vhandle[2]);
		face_vhandles.push_back(vhandle[3]);
		quadPartition.add_face(face_vhandles);

		std::vector<MyMesh::VertexHandle> taggedTemp;
		for (size_t i = 0; i < tagged.size(); i++)
		{
			if (i != idx1 && i != idx2) {
				taggedTemp.push_back(tagged[i]);
			}
		}
		tagged = taggedTemp;
		//count--;
	}
}
bool MyViewer::isLeft(MyMesh::Point a, MyMesh::Point b, MyMesh::Point c) {
	return ((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])) > 0;
}