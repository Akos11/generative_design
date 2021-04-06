#include "MyViewer.h"
#include <OpenMesh/Tools/Smoother/JacobiLaplaceSmootherT.hh>
#include <iostream>
#include <gsl/gsl_linalg.h>
#include <iostream>
#include <fstream>

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
/// For every quad perform a flip if we get better valences that way
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
	double angleThreshold = (3.14159 / 2.0) * 1.5f;
	initiateBoundaryConstraints(angleThreshold);
}

void MyViewer::initiateBoundaryConstraints(double angleThreshold) {
	for (auto v : mesh.vertices()) {
		if (mesh.is_boundary(v)) {
			MyMesh::Point neighbour1;
			MyMesh::Point neighbour2;
			bool found = false;
			for (MyMesh::VertexVertexIter vv_iter = mesh.vv_iter(v); vv_iter.is_valid(); vv_iter++) {
				if (mesh.is_boundary(vv_iter.handle()) && mesh.is_boundary(getCommonEdge(v,vv_iter.handle()))) {
					if (!found) {
						neighbour1 = mesh.point(vv_iter.handle());
						found = true;
					}
					else {
						neighbour2 = mesh.point(vv_iter.handle());
						break;
					}
				}
			}
			OpenMesh::Vec3d vNeighbour0 = (neighbour1 - mesh.point(v)).normalized();
			OpenMesh::Vec3d vNeighbour1 = (neighbour2 - mesh.point(v)).normalized();

			double dotProduct = dot(vNeighbour0, vNeighbour1);
			dotProduct = fmaxf(-1, dotProduct);
			dotProduct = fminf(1, dotProduct);
			double alpha = acos(dotProduct);
			//printf("%f %f\n", dotProduct, alpha);
			double randomNumber = ((double)rand()/ RAND_MAX )*2- 1;
			double randomNumber2 = ((double)rand() / RAND_MAX) * 2 - 1;
			mesh.data(v).u = OpenMesh::Vec3d(0, 0,0.0);
			OpenMesh::Vec3d temp = vNeighbour0 + vNeighbour1;
			double wtemp[3];
			MyMesh::VertexHandle v_arraytemp[3];
			Vector tempMiddle = neighbour1 + 0.5 * (neighbour2 - neighbour1);
			MyMesh::FaceHandle tempF2 = getFace(tempMiddle, wtemp, v_arraytemp);
			
			if (alpha > angleThreshold) {
				if (abs(alpha - 3.14159) < 0.001)
					temp = vNeighbour0;
				temp = temp.normalized();
			} 
			else {
				temp = temp.normalized();
				OpenMesh::Vec3d temp2 = OpenMesh::Vec3d(temp[1],-temp[0],0);
				temp = temp + temp2;
				temp = temp.normalized();
				if (wtemp[0] < -0.5) {
					double wtemp[3];
					MyMesh::VertexHandle v_arraytemp[3];
					MyMesh::FaceHandle tempF = getFace(mesh.point(v), wtemp, v_arraytemp);
					singularities.push_back(Singularity(mesh.point(v), tempF));
					corners.push_back(Corner(mesh.point(v), std::vector<SeparatricePart>(), true, v, v));
				}
				else {
					corners.push_back(Corner(mesh.point(v), std::vector<SeparatricePart>(), true, v, v,true));
					//temp = vNeighbour0;
					//temp = temp.normalized();
				}
			}
			if (temp[0] < 0.0)
			{
				temp[0] = temp[0] * -1;
				temp[1] = temp[1] * -1;
			}
			if (temp[1] < 0.0)
			{
				double t = temp[0];
				temp[0] = temp[1] * -1;
				temp[1] = t;
			}
			double theta = acos(temp[0]);
			theta = theta * 4.0;
			mesh.data(v).u = OpenMesh::Vec3d(cos(theta), sin(theta), 0.0);
			
		}
	}
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
void MyViewer::initVFunction(int iterations) {
	clock_t start = clock() / (CLOCKS_PER_SEC / 1000);
	MyMesh::VertexHandle* omega = new MyMesh::VertexHandle[mesh.n_vertices()];

	std::vector<int>* neighbours = new std::vector<int>[mesh.n_vertices()];

	//V
	double* vFunc = new double[mesh.n_vertices()];
	double* gradVFunc = new double[mesh.n_vertices()];
	//U
	double* uFunc = new double[mesh.n_vertices()];
	double* gradUFunc = new double[mesh.n_vertices()];
	//U0
	double* u0Func = new double[mesh.n_vertices()];
	double* gradU0Func = new double[mesh.n_vertices()];

	int* matrixIdx = new int[mesh.n_vertices()];

	int rowNum = 0;
	for (auto v : mesh.vertices()) {
		if (!mesh.is_boundary(v))
			rowNum++;
	}
	double* A_data = new double[rowNum * rowNum * 4];
	double* b_data = new double[rowNum * 2];
	for (size_t Ui = 0; Ui < 2; Ui++)
	{
		int idx = 0;
		int mIdx = 0;
		for (auto v : mesh.vertices()) {
			omega[idx] = v;
			vFunc[idx] = 0;
			gradVFunc[idx] = 0;
			if (mesh.is_boundary(v)) {
				u0Func[idx] = mesh.data(v).u[Ui];
				matrixIdx[idx] = 0;
			}
			else {
				u0Func[idx] = 0;
				matrixIdx[idx] = mIdx++;
			}
			uFunc[idx] = 0;
			gradU0Func[idx] = 0;
			gradUFunc[idx] = 0;
			idx++;
		}
		buildNeighbours(omega, neighbours);
		calculateGrad(omega, neighbours, u0Func, gradU0Func);

		int rowIdx = 0;
		for (int i = 0; i < mesh.n_vertices(); i++) {
			start = clock() / (CLOCKS_PER_SEC / 1000);
			if (mesh.is_boundary(omega[i]))
				continue;
			for (int j = 0; j < mesh.n_vertices(); j++) {
				vFunc[j] = 0;
				gradVFunc[j] = 0;
				A_data[(rowIdx*2 + Ui) * rowNum*2 + matrixIdx[j]*2 +Ui] = 0;
			}
			vFunc[i] = 1;
			calculateGrad(omega, neighbours, vFunc, gradVFunc);
			double sumRight = 0;

			for (int j = 0; j < mesh.n_vertices(); j++) {
				if (gradVFunc[j] != 0) {
					sumRight += gradU0Func[j] * gradVFunc[j];
				}
				double L = calculateL(omega, neighbours, j);
				for (auto n : neighbours[j]) {
					if (!mesh.is_boundary(omega[n])) {
						A_data[(rowIdx * 2 + Ui) * rowNum*2 + matrixIdx[n] * 2 + Ui] += (1 / (mesh.point(omega[n]) - mesh.point(omega[j])).length()) * gradVFunc[j];
					}
				}
				if (!mesh.is_boundary(omega[j]))
					A_data[(rowIdx * 2 + Ui) * rowNum*2 + matrixIdx[j] * 2 + Ui] += -1 * gradVFunc[j] * L;
			}
			b_data[(rowIdx * 2 + Ui)] = -sumRight;
			rowIdx++;

		}
	}
	auto A = gsl_matrix_view_array(A_data, rowNum*2, rowNum*2);
	auto b = gsl_vector_view_array(b_data, rowNum*2);
	gsl_vector* x = gsl_vector_alloc(rowNum*2);
	gsl_permutation * perm = gsl_permutation_alloc(rowNum*2);

	// Decompose A into the LU form:
	int signum;       // Sign of the permutation
	qDebug() << "LU decomp";
	gsl_linalg_LU_decomp(&A.matrix, perm, &signum);
	// Solve the linear system
	qDebug() << "LU solve";
	gsl_linalg_LU_solve(&A.matrix, perm, &b.vector, x);

	double* A_data2 = new double[rowNum * 4 *rowNum * 4];
	double* b_data2 = new double[rowNum * 4];
	gsl_vector* x2 = gsl_vector_alloc(rowNum * 4);
	gsl_permutation* perm2 = gsl_permutation_alloc(rowNum * 4);

	for (size_t iter = 0; iter < iterations; iter++)
	{
		for (size_t i = 0; i < rowNum * 4; i++)
		{
			for (size_t j = 0; j < rowNum * 4; j++)
			{
				if (i < rowNum * 2 && j < rowNum * 2) {
					A_data2[i * rowNum * 4 + j] = A_data[i * rowNum * 2 + j];
					b_data2[i] = b_data[i];
				}
				else if (i < rowNum * 2 && j >= rowNum * 2) {
					if (i == j - rowNum * 2) {
						A_data2[i * rowNum * 4 + j] = iter == 0 ? x->data[i] : x2->data[i];
					}
					else
						A_data2[i * rowNum * 4 + j] = 0;
				}
				else if (i >= rowNum * 2 && j < rowNum * 2) {
					if (i - rowNum * 2 == j)
						A_data2[i * rowNum * 4 + j] = iter == 0 ? x->data[j] : x2->data[j];
					else
						A_data2[i * rowNum * 4 + j] = 0;
					b_data2[i] = 1;
				}
				else if (i >= rowNum * 2 && j >= rowNum * 2) {
					A_data2[i * rowNum * 4 + j] = 0;
				}
			}
		}
		auto A2 = gsl_matrix_view_array(A_data2, rowNum * 4, rowNum * 4);
		auto b2 = gsl_vector_view_array(b_data2, rowNum * 4);


		// Decompose A into the LU form:
		int signum2;       // Sign of the permutation
		qDebug() << "LU decomp";
		gsl_linalg_LU_decomp(&A2.matrix, perm2, &signum2);
		// Solve the linear system
		qDebug() << "LU solve";
		gsl_linalg_LU_solve(&A2.matrix, perm2, &b2.vector, x2);
	}
	

	for (size_t Ui = 0; Ui < 2; Ui++)
	{
		for (size_t i = 0; i < mesh.n_vertices(); i++)
		{
			if (!mesh.is_boundary(omega[i])) {
				mesh.data(omega[i]).u[Ui] = x2->data[matrixIdx[i] * 2 + Ui];
			}
		}
	}
	gsl_vector_free(x2);
	gsl_permutation_free(perm2);
	gsl_vector_free(x);
	gsl_permutation_free(perm);

	delete[] omega;
	delete[] vFunc;
	delete[] gradVFunc;
	delete[] uFunc;
	delete[] gradUFunc;
	delete[] u0Func;
	delete[] gradU0Func;
	//delete[] uCoefficient;
	delete[] A_data;
	delete[] b_data;
	delete[] A_data2;
	delete[] b_data2;
	for (auto v : mesh.vertices())
	{
		mesh.data(v).u = mesh.data(v).u.normalize();
	}
}

void MyViewer::calculateGrad(MyMesh::VertexHandle* omega, std::vector<int>* neighbours, double* func, double* gradFunc) {

	for (size_t i = 0; i < mesh.n_vertices(); i++)
	{
		double grad = 0;
		for (auto n : neighbours[i]) {
			grad += (func[n] - func[i]) / (mesh.point(omega[n]) - mesh.point(omega[i])).length();
		}
		gradFunc[i] = grad;
	}
}
std::vector<int> MyViewer::getNeighbours(MyMesh::VertexHandle* omega, int idx) {
	std::vector<int> neighbours;
	for (size_t j = 0; j < mesh.n_vertices(); j++)
	{
		if (idx != j && hasCommonEdge(omega[idx], omega[j])) {
			neighbours.push_back(j);
		}
	}
	return neighbours;
}
double MyViewer::calculateL(MyMesh::VertexHandle* omega, std::vector<int>* neighbours, int idx) {
	double L = 0;
	for (auto n : neighbours[idx]) {
		L += 1 / ((mesh.point(omega[n]) - mesh.point(omega[idx])).length());
	}
	return L;
}
void MyViewer::buildNeighbours(MyMesh::VertexHandle* omega, std::vector<int>* neighbours) {
	for (size_t i = 0; i < mesh.n_vertices(); i++)
	{
		neighbours[i] = getNeighbours(omega, i);
	}
}
void MyViewer::findSingularities() {

	double* A_data = new double[3 * 3];
	double* b_data = new double[3];

	gsl_vector* x = gsl_vector_alloc(3);
	gsl_permutation* perm = gsl_permutation_alloc(3);
	for (auto f : mesh.faces()) {
		mesh.data(f).hasSingularity = false;
		MyMesh::VertexHandle v1;
		MyMesh::VertexHandle v2;
		MyMesh::VertexHandle v3;
		int i = 0;
		for (auto v : mesh.fv_range(f))
		{
			switch (i)
			{
			case 0 :
				v1 = v;
				i++;
				break;
			case 1:
				v2 = v;
				i++;
				break;
			case 2:
				v3 = v;
				i++;
				break;
			}
		}
		if (i < 3)
			throw;
		//Row 0
		A_data[0] = 1;
		A_data[1] = 1;
		A_data[2] = 1;
		b_data[0] = 1;

		//Row 1
		A_data[3] = mesh.data(v1).u[0];
		A_data[4] = mesh.data(v2).u[0];
		A_data[5] = mesh.data(v3).u[0];
		b_data[1] = 0;

		//Row 2
		A_data[6] = mesh.data(v1).u[1];
		A_data[7] = mesh.data(v2).u[1];
		A_data[8] = mesh.data(v3).u[1];
		b_data[2] = 0;
		auto A = gsl_matrix_view_array(A_data, 3, 3);
		auto b = gsl_vector_view_array(b_data, 3);


		// Decompose A into the LU form:
		int signum;       // Sign of the permutation
		gsl_linalg_LU_decomp(&A.matrix, perm, &signum);
		double det = gsl_linalg_LU_det(&A.matrix, signum);
		//qDebug() << det;
		if (fabs(det) > 0.00001) {
			gsl_linalg_LU_solve(&A.matrix, perm, &b.vector, x);
			if (x->data[0] > 0 && x->data[1] > 0 && x->data[2] > 0) {
				//qDebug() << "w1: " << x->data[0] << "w2: " << x->data[1] << "w3: " << x->data[2];
				//qDebug() << "x: " << mesh.data(v1).u[0] * x->data[0] + mesh.data(v2).u[0] * x->data[1] + mesh.data(v3).u[0] * x->data[2];
				//qDebug() << "y: " << mesh.data(v1).u[1] * x->data[0] + mesh.data(v2).u[1] * x->data[1] + mesh.data(v3).u[1] * x->data[2];
				//mesh.data(f).tagged2 = true;
				MyMesh::Point p = mesh.point(v1) * x->data[0] + mesh.point(v2) * x->data[1] + mesh.point(v3) * x->data[2];
				//qDebug() << "sing x: " << p[0] << " y: " << p[1];
				//for (double yd = -0.2; yd <= 0.2; yd+=0.1)
				//{
				//	double w_array[3];
				//	MyMesh::VertexHandle v_array[3];
				//	Vector temp = p;
				//	temp[1] += yd;
				//	MyMesh::FaceHandle f = getFace(temp, w_array, v_array);
				//	Vector Utemp = Vector(0, 0, 0);
				//	for (size_t i = 0; i < 3; i++)
				//	{
				//		Utemp += w_array[i] * mesh.data(v_array[i]).u;
				//	}
				//	qDebug() << "yd : " << yd << " U: x:" << Utemp[0] << " y: " << Utemp[1];
				//}
				singularities.push_back(Singularity(p,f));
				mesh.data(f).hasSingularity = true;
			}
		}

		

	}
	gsl_vector_free(x);
	gsl_permutation_free(perm);
	delete[] A_data;
	delete[] b_data;
}

void MyViewer::findSeparatrices() {

	for (auto sin : singularities) {
		MyMesh::VertexHandle vertices[6];
		MyMesh::VertexHandle v1;
		MyMesh::VertexHandle v2;
		MyMesh::VertexHandle v3;
		int i = 0;
		for (auto v : mesh.fv_range(sin.f))
		{
			switch (i)
			{
			case 0:
				v1 = v;
				i++;
				break;
			case 1:
				v2 = v;
				i++;
				break;
			case 2:
				v3 = v;
				i++;
				break;
			}
		}
		if (i < 3)
			throw;
		vertices[0] = v1; vertices[1] = v2; vertices[2] = v1; vertices[3] = v3; vertices[4] = v2; vertices[5] = v3;
		for (size_t vi = 0; vi < 6; vi+=2)
		{
			double minDot[4];
			minDot[0] = 999999999;
			minDot[1] = 999999999;
			minDot[2] = 999999999;
			minDot[3] = 999999999;

			double minW[4];
			minW[0] = -1;
			minW[1] = -1;
			minW[2] = -1;
			minW[3] = -1;

			Vector minCross[4];
			for (double w = 0; w < 1.0; w+=0.001)
			{
				MyMesh::Point P = mesh.point(vertices[vi]) * w + (1-w) * mesh.point(vertices[vi+1]);
				Vector UP = (P - sin.pos).normalize();
				Vector U = mesh.data(vertices[vi]).u * w + (1 - w) * mesh.data(vertices[vi + 1]).u;
				U = U.normalize();
				Vector cross[4];
				cross[0] = (calculateCrossFromU(U)[0]).normalize();
				cross[1] = (calculateCrossFromU(U)[1]).normalize();
				cross[2] = cross[0] * -1;
				cross[3] = cross[1] * -1;

				for (size_t i = 0; i < 4; i++)
				{
					if ((1 - dot(cross[i], UP)) < minDot[i]) {
						minDot[i] = (1 - dot(cross[i], UP));
						minW[i] = w;
						minCross[i] = cross[i];
					}
				}
				{
				
				//MyMesh::Point P = mesh.point(vertices[vi]) * w + (1-w) * mesh.point(vertices[vi+1]);
				//Vector UP = (P - sin.pos).normalize();
				//Vector U = mesh.data(vertices[vi]).u * w + (1 - w) * mesh.data(vertices[vi + 1]).u;
				//U = U.normalize();
				//Vector cross1[4];
				//cross1[0] = (calculateCrossFromU(mesh.data(vertices[vi]).u)[0]).normalize();
				//cross1[1] = (calculateCrossFromU(mesh.data(vertices[vi]).u)[1]).normalize();
				//cross1[2] = cross1[0] * -1;
				//cross1[3] = cross1[1] * -1;
				//Vector cross2[4];
				//cross2[0] = (calculateCrossFromU(mesh.data(vertices[vi + 1]).u)[0]).normalize();
				//cross2[1] = (calculateCrossFromU(mesh.data(vertices[vi + 1]).u)[1]).normalize();
				//cross2[2] = cross2[0] * -1;
				//cross2[3] = cross2[1] * -1;
				//Vector cross[4];
				//cross[0] = (w * cross1[0] + (1 - w) * cross2[0]).normalize();
				//cross[1] = (w * cross1[1] + (1 - w) * cross2[1]).normalize();
				//cross[2] = (w * cross1[2] + (1 - w) * cross2[2]).normalize();
				//cross[3] = (w * cross1[3] + (1 - w) * cross2[3]).normalize();

				//for (size_t i = 0; i < 4; i++)
				//{
				//	if ((1 - dot(cross[i], UP)) < minDot[i]) {
				//		minDot[i] = (1 - dot(cross[i], UP));
				//		minW[i] = w;
				//		minCross[i] = cross[i];
				//	}
				//}
				
				}
			}
			for (size_t i = 0; i < 4; i++)
			{
				if (minDot[i] < 0.00001) {
					std::vector<Vector> separatrice;
					separatrice.push_back(sin.pos);
					MyMesh::Point P = mesh.point(vertices[vi]) * minW[i] + (1 - minW[i]) * mesh.point(vertices[vi + 1]);
					separatrice.push_back(P);
					//separatrice.push_back(P + minCross[i] * 5);

					buildSeparatrices(&separatrice, minCross[i], vertices[vi], vertices[vi + 1], sin.f);
					separatrices.push_back(separatrice);
				}
			}
		}
		/*for (size_t vi = 0; vi < 6; vi+=2)
		{
			Vector CFS1[2];
			Vector CFS2[2];
			CFS1[0] = calculateCrossFromU(mesh.data(vertices[vi]).u)[0].normalize();
			CFS1[1] = calculateCrossFromU(mesh.data(vertices[vi]).u)[1].normalize();
			CFS2[0] = calculateCrossFromU(mesh.data(vertices[vi +1]).u)[0].normalize();
			CFS2[1] = calculateCrossFromU(mesh.data(vertices[vi +1]).u)[1].normalize();
			for (size_t i = 0; i < 2; i++)
			{
				double a = mesh.point(vertices[vi +1])[0] - sin.pos[0];
				double b = CFS2[i][0];
				double c = mesh.point(vertices[vi])[0] - mesh.point(vertices[vi +1])[0];
				double d = CFS1[i][0] - CFS2[i][0];

				double e = mesh.point(vertices[vi])[1] - mesh.point(vertices[vi +1])[1];
				double f = CFS1[i][1] - CFS2[i][1];
				double g = mesh.point(vertices[vi +1])[1] - sin.pos[1];
				double h = CFS2[i][1];

				double A = a * e - g * c;
				double B = h * c + g * d - b * e - f * a;
				double C = -h * d + f * b;
				double det = B * B - 4 * A * C;
				if (fabs(det) < 0.001) {
					double k1 = -B / (2 * A);
					//qDebug() << "k1: " << k1;
					double w = (-k1 * a + b) / (k1 * c - d);
					//qDebug() << "w: " << w;
					if (w >= 0 && w <= 1) {
						std::vector<Vector> separatrice;
						separatrice.push_back(sin.pos);
						Vector temp = mesh.point(vertices[vi]) * w + mesh.point(vertices[vi +1]) * (1 - w);
						Vector dirtemp = CFS1[i] * w + CFS2[i] * (1 - w);
						dirtemp = dirtemp.normalize();
						if (k1 < 0)
							dirtemp = dirtemp.normalize()*-1;
						separatrice.push_back(temp);
						buildSeparatrices(&separatrice, dirtemp, vertices[vi], vertices[vi + 1], sin.f);
						separatrices.push_back(separatrice);
					}
				}
				else if (det > 0) {
					double k1 = (-B + sqrt(det)) / (2 * A);
					double k2 = (-B - sqrt(det)) / (2 * A);
					//qDebug() << "k1: " << k1 << "k2: " << k2;
					double w1 = (-k1 * a + b) / (k1 * c - d);
					double w2 = (-k2 * a + b) / (k2 * c - d);
					//qDebug() << "w1: " << w1 << " w2: " << w2;
					if (w1 >= 0 && w1 <= 1) {
						std::vector<Vector> separatrice;
						separatrice.push_back(sin.pos);
						Vector temp = mesh.point(vertices[vi]) * w1 + mesh.point(vertices[vi +1]) * (1 - w1);
						Vector dirtemp = CFS1[i] * w1 + CFS2[i] * (1 - w1);
						dirtemp = dirtemp.normalize();
						if (k1 < 0)
							dirtemp = dirtemp.normalize()*-1;
						separatrice.push_back(temp);
						//Vector U = w1 * mesh.data(vertices[vi]).u + (1 - w1) * mesh.data(vertices[vi + 1]).u;
					//	qDebug() << "dirTemp x: " << dirtemp[0] << " y: " << dirtemp[1];
						buildSeparatrices(&separatrice, dirtemp, vertices[vi], vertices[vi+1],sin.f);
						separatrices.push_back(separatrice);
					}
					if (w2 >= 0 && w2 <= 1) {
						std::vector<Vector> separatrice;
						separatrice.push_back(sin.pos);
						Vector temp = mesh.point(vertices[vi]) * w2 + mesh.point(vertices[vi +1]) * (1 - w2);
						Vector dirtemp = CFS1[i] * w2 + CFS2[i] * (1 - w2);
						dirtemp = dirtemp.normalize();
						if (k2 < 0)
							dirtemp = dirtemp.normalize()*-1;
						separatrice.push_back(temp);
						//Vector U = w2 * mesh.data(vertices[vi]).u + (1 - w2) * mesh.data(vertices[vi + 1]).u;

						//qDebug() << "dirTemp x: " << dirtemp[0] << " y: " << dirtemp[1];
						buildSeparatrices(&separatrice, dirtemp, vertices[vi], vertices[vi + 1], sin.f);
						separatrices.push_back(separatrice);
					}
				}
			}
		}
		*/
	}
}
void MyViewer::findSeparatrices2() {
	double r = 0.05;
	for (auto sing : singularities) {
		int countSep = 0;
		for (double alpha = 0; alpha < 3.14159*2; alpha+=0.002)
		{
			Vector P = sing.pos + Vector(cos(alpha) * r, sin(alpha) * r,0);
			double w[3];
			MyMesh::VertexHandle v_array[3];
			MyMesh::FaceHandle f = getFace(P, w, v_array);
			Vector SP = P - sing.pos;
			SP = SP.normalize();
			Vector U = Vector(0, 0, 0);
			for (size_t i = 0; i < 3; i++)
			{
				U += w[i] * mesh.data(v_array[i]).u;
			}
			U = U.normalize();
			Vector cross[4];
			cross[0] = (calculateCrossFromU(U)[0]).normalize();
			cross[1] = (calculateCrossFromU(U)[1]).normalize();
			cross[2] = cross[0] * -1;
			cross[3] = cross[1] * -1;

			double minDif = 99999;
			int minIdx = 0;
			for (size_t j = 0; j < 4; j++)
			{
				double temp = dot(cross[j], SP);
				if (1 - temp < minDif) {
					minDif = 1 - temp;
					minIdx = j;
				}
			}
			if (minDif < 0.00001) {
				std::vector<Vector> separatrice;
				separatrice.push_back(sing.pos);
				separatrice.push_back(P);
				followStreamLine(P, cross[minIdx], &separatrice, separatrices.size(),0.1,2000);

				separatrices.push_back(separatrice);
				countSep++;
				alpha += 0.1;
			}

		}
		qDebug() << countSep;
	}
	//eliminateDuplicateSeparatrices(0.5);
}
void MyViewer::eliminateDuplicateSeparatrices(double limit) {
	struct singularityPair {
		int separatriceIdx;
		int singularityIdx1;
		int singularityIdx2;
		int separatriceT;
		int separatriceIdx2;
		int separatriceT2;
		bool added;

	};
	int idx = 0;
	double stepSize = 0.1;
	std::vector<singularityPair> singularityPairs;
	for (auto sep : separatrices) {
		int startI = (limit / stepSize) + 1;
		int startingSing = 0;
		for (auto sing : singularities) {
			if (getDistance(sing.pos, sep[0]) < 0.001) {
				break;
			}
			startingSing++;
		}
		int endingSing = -1;
		int separatriceT = -1;
		for (size_t i = 0; i < sep.size(); i++)
		{
			for (size_t j = 0; j < singularities.size(); j++)
			{
				if (j != startingSing && getDistance(singularities[j].pos,sep[i]) < limit ) {
					endingSing = j;
					separatriceT = i;
					singularityPair temp;
					temp.separatriceIdx = idx;
					temp.separatriceIdx2 = -1;
					temp.singularityIdx1 = startingSing;
					temp.singularityIdx2 = endingSing;
					temp.separatriceT = separatriceT;
					temp.separatriceT2 = -1;
					temp.added = false;
					singularityPairs.push_back(temp);
					//break;
				}
			}
			if (endingSing != -1)
				break;
		}

		idx++;
	}
	//for (auto singPair : singularityPairs) {
	//	qDebug() << "sep Idx: " << singPair.separatriceIdx << " sep T: " << singPair.separatriceT << " sing1: " << singPair.singularityIdx1 << " sing2: " << singPair.singularityIdx2;
	//}
	for (size_t i = 0; i < singularityPairs.size(); i++)
	{
		for (size_t j = 0; j < singularityPairs.size(); j++)
		{
			if (singularityPairs[i].singularityIdx1 == singularityPairs[j].singularityIdx2 &&
				singularityPairs[i].singularityIdx2 == singularityPairs[j].singularityIdx1 && !singularityPairs[i].added)
			{
				singularityPairs[i].separatriceIdx2 = singularityPairs[j].separatriceIdx;
				singularityPairs[i].separatriceT2 = singularityPairs[j].separatriceT;
				singularityPairs[i].added = true;
				singularityPairs[j].added = true;
			}
		}
	}
//	qDebug() << "####################";
	std::vector<int> separatricesToDelete;
	std::vector<int> cornersToDelete;
	for (size_t i = 0; i < singularityPairs.size(); i++)
	{
		if (singularityPairs[i].separatriceIdx2 != -1) {
			for (size_t cornerIdx = 0; cornerIdx < corners.size(); cornerIdx++)
			{
				if (corners[cornerIdx].separatrices.size() > 0 && corners[cornerIdx].separatrices[0].separatriceIdx == singularityPairs[i].separatriceIdx &&
					corners[cornerIdx].separatrices[0].fromI >= singularityPairs[i].separatriceT) {
					cornersToDelete.push_back(cornerIdx);
				}
			}
			/*qDebug() << "sep Idx1: " << singularityPairs[i].separatriceIdx << " sep T: " << singularityPairs[i].separatriceT << " sing1: " << singularityPairs[i].singularityIdx1 << " sing2: " << singularityPairs[i].singularityIdx2
				<< "sep Idx2: " << singularityPairs[i].separatriceIdx2<< " sep T2: " << singularityPairs[i].separatriceT2;*/
			auto startIt = separatrices[singularityPairs[i].separatriceIdx].begin() + singularityPairs[i].separatriceT + 1;
			separatrices[singularityPairs[i].separatriceIdx].erase(startIt, separatrices[singularityPairs[i].separatriceIdx].end());
			separatrices[singularityPairs[i].separatriceIdx][singularityPairs[i].separatriceT] = singularities[singularityPairs[i].singularityIdx2].pos;
			separatricesToDelete.push_back(singularityPairs[i].separatriceIdx2);
			//qDebug() << "Sep edit: " << singularityPairs[i].separatriceIdx;
			if (singularityPairs[i].separatriceT == 0 || separatrices[singularityPairs[i].separatriceIdx].size() == 0)
				separatricesToDelete.push_back(singularityPairs[i].separatriceIdx);		
		}
	}
	std::sort(separatricesToDelete.begin(), separatricesToDelete.end(), std::greater<int>());
	separatricesToDelete.erase(std::unique(separatricesToDelete.begin(), separatricesToDelete.end()), separatricesToDelete.end());
	for (auto d : separatricesToDelete) {
		separatrices.erase(separatrices.begin() + d);
		for (size_t cornerIdx = 0; cornerIdx < corners.size(); cornerIdx++)
		{
			std::vector<int> toDelete;
			int tempIdx = 0;
			for (int tempIdx = 0; tempIdx < corners[cornerIdx].separatrices.size(); tempIdx++) {
				if (corners[cornerIdx].separatrices[tempIdx].separatriceIdx == d) {
					toDelete.push_back(tempIdx);
				}
				if (corners[cornerIdx].separatrices[tempIdx].separatriceIdx > d)
					corners[cornerIdx].separatrices[tempIdx].separatriceIdx--;
			}
			std::sort(toDelete.begin(), toDelete.end(), std::greater<int>());
			for (auto cd : toDelete) {
				corners[cornerIdx].separatrices.erase(corners[cornerIdx].separatrices.begin() + cd);
			}
		}
	}
	std::sort(cornersToDelete.begin(), cornersToDelete.end(), std::greater<int>());
	cornersToDelete.erase(std::unique(cornersToDelete.begin(), cornersToDelete.end()), cornersToDelete.end());
	for (auto cd : cornersToDelete) {
		corners.erase(corners.begin() + cd);
	}
	std::vector<int> singularitySepNum(singularities.size());
	for (size_t i = 0; i < singularities.size(); i++)
	{
		singularitySepNum[i] = 0;
	}
	int sepIdx = 0;
	for (auto sep : separatrices) {
		for (size_t i = 0; i < singularities.size(); i++)
		{
			if (getDistance(singularities[i].pos, sep[0]) < 0.0001 || getDistance(singularities[i].pos, sep[sep.size()-1]) < 0.0001) {
				singularitySepNum[i]++;
			}
		}
		sepIdx++;
	}
	std::vector<int> singularitiesToDelete;
	for (size_t i = 0; i < singularities.size(); i++)
	{
		//qDebug() << "i: " << i << " SepNum: " << singularitySepNum[i] << "pos: (" << singularities[i].pos[0] << ", " << singularities[i].pos[1] <<")";
		if (singularitySepNum[i] == 0) {
			singularitiesToDelete.push_back(i);
		}
	}
	std::sort(singularitiesToDelete.begin(), singularitiesToDelete.end(), std::greater<int>());
	for (auto d : singularitiesToDelete) {
		singularities.erase(singularities.begin() + d);
	}
	refreshSeparatricePartsAtFaces();
}
void MyViewer::refreshSeparatricePartsAtFaces() {
	for (auto f : mesh.faces()) {
		mesh.data(f).separaticeParts.clear();
	}
	double wtemp[3];
	MyMesh::VertexHandle v_arraytemp[3];
	int idx = 0;
	for (auto sep : separatrices) {
		MyMesh::FaceHandle prevFace = getFace(sep[0], wtemp, v_arraytemp);
		int startIdx = 0;
		for (int i = 1; i < sep.size(); i++) {
			MyMesh::FaceHandle tempFace = getFace(sep[i], wtemp, v_arraytemp);
			if (tempFace != prevFace || i == sep.size() - 1) {
				mesh.data(prevFace).separaticeParts.push_back(SeparatricePart(idx,startIdx,i < sep.size() -1 ? i+1 : i));
				if (tempFace != prevFace && i == sep.size() - 1) {
					mesh.data(tempFace).separaticeParts.push_back(SeparatricePart(idx, i-1, i));
				}
				prevFace = tempFace;
				startIdx = i - 1;
			}
		}
		idx++;
	}
}
void MyViewer::buildSeparatrices(std::vector<Vector>* separatice, Vector dir,MyMesh::VertexHandle v1, MyMesh::VertexHandle v2, MyMesh::FaceHandle f) {
	MyMesh::VertexHandle v1_ = v1;
	MyMesh::VertexHandle v2_ = v2;
	Vector dir_ = dir;
	MyMesh::FaceHandle f_ = f;
	for (size_t i = 0; i < 1000; i++)
	{


		MyMesh::HalfedgeHandle h = mesh.halfedge_handle(getCommonEdge(v1_, v2_), 0);
		if (mesh.face_handle(h) == f_)
			h = mesh.halfedge_handle(getCommonEdge(v1_, v2_), 1);
		MyMesh::FaceHandle nextF = mesh.face_handle(h);
		f_ = nextF;
		MyMesh::VertexHandle v3;
		for (auto vertex : mesh.fv_range(nextF)) {
			if (vertex != v1_ && vertex != v2_)
				v3 = vertex;
		}
		Vector tempXin;
		Vector tempDin;
		MyMesh::VertexHandle tempV1;
		MyMesh::VertexHandle tempV2;
		findXiNext((*separatice)[separatice->size() - 1], dir_, v1_, v2_, v3, &tempXin, &tempDin, &tempV1, &tempV2);

		Vector din = ((tempDin.normalize() + dir_.normalize()) / 2).normalize();
		Vector Xin;
		findXiNext((*separatice)[separatice->size() - 1], din, v1_, v2_, v3, &Xin, &dir_, &v1_, &v2_);

		separatice->push_back(Xin);
		if (mesh.is_boundary(v1_) && mesh.is_boundary(v2_))
			break;
		//if (mesh.data(f_).hasSingularity)
		//	break;
	}
}
std::vector<MyViewer::Vector> MyViewer::calculateCrossFromU(Vector u) {
	std::vector<MyViewer::Vector> returnV;
	double atang = atan2(u[1], u[0]);
	if (atang < 0)
		atang = 3.14159 * 2 + atang;
	double theta = (atang) / 4.0;
	returnV.push_back(Vector(cos(theta), sin(theta), 0.0));
	returnV.push_back(Vector(cos(theta + 3.14159 / 2), sin(theta + 3.14159 / 2), 0.0) );
	return returnV;
}
MyViewer::Vector MyViewer::findClosestCrossVector(Vector di, Vector Xi, MyMesh::VertexHandle v1, MyMesh::VertexHandle v2) {
	double w;
	if ((mesh.point(v1)[0] - mesh.point(v2)[0]) != 0)
		w = (Xi[0] - mesh.point(v2)[0]) / (mesh.point(v1)[0] - mesh.point(v2)[0]);
	else
		w = (Xi[1] - mesh.point(v2)[1]) / (mesh.point(v1)[1] - mesh.point(v2)[1]);
	Vector U = w * mesh.data(v1).u + (1 - w) * mesh.data(v2).u;
	U = U.normalize();
	Vector cross[4];
	cross[0] = (calculateCrossFromU(U)[0]).normalize();
	cross[1] = (calculateCrossFromU(U)[1]).normalize();
	cross[2] = cross[0] * -1;
	cross[3] = cross[1] * -1;

	//qDebug() << "cross1 x: " << cross[0][0] << " y: " << cross[0][1];
	//qDebug() << "cross2 x: " << cross[1][0] << " y: " << cross[1][1];
	//qDebug() << "cross3 x: " << cross[2][0] << " y: " << cross[2][1];
	//qDebug() << "cross4 x: " << cross[3][0] << " y: " << cross[3][1];
	//qDebug() << "#########################";
	di = di.normalize();
	double minDif = 99999;
	int minIdx = 0;
	for (size_t i = 0; i < 4; i++)
	{
		double temp = dot(cross[i], di);
		if (1 - temp < minDif) {
			minDif = 1 - temp;
			minIdx = i;
		}
	}
	Vector Vdi = cross[minIdx];
	return Vdi;
}
void MyViewer::findXiNext(Vector Xi, Vector di, MyMesh::VertexHandle v1, MyMesh::VertexHandle v2, MyMesh::VertexHandle v3, Vector* Xin, Vector* din, MyMesh::VertexHandle* newV1, MyMesh::VertexHandle* newV2) {
	Vector Vdi = findClosestCrossVector(di, Xi, v1, v2);

	//separatice->push_back((*separatice)[separatice->size() - 1] + Vdi * 5);
	double a = (mesh.point(v1)[0] - mesh.point(v3)[0]) / Vdi[0];
	double b = (mesh.point(v3)[0] - Xi[0]) / Vdi[0];
	double w = (mesh.point(v3)[1] - Xi[1] - b * Vdi[1]) / (Vdi[1] * a - mesh.point(v1)[1] + mesh.point(v3)[1]);
	if (w >= 0 && w <= 1) {
		*Xin = (mesh.point(v1) * w + (1 - w) * mesh.point(v3));
		*din = findClosestCrossVector(di, *Xin, v1, v3);
		*newV1 = v1;
		*newV2 = v3;
	}
	else {
		a = (mesh.point(v2)[0] - mesh.point(v3)[0]) / Vdi[0];
		b = (mesh.point(v3)[0] - Xi[0]) / Vdi[0];
		double w = (mesh.point(v3)[1] - Xi[1] - b * Vdi[1]) / (Vdi[1] * a - mesh.point(v2)[1] + mesh.point(v3)[1]);
		if (w >= 0 && w <= 1) {
			*Xin = (mesh.point(v2) * w + (1 - w) * mesh.point(v3));
			*din = findClosestCrossVector(di, *Xin, v2, v3);
			*newV1 = v2;
			*newV2 = v3;
		}
	}
}
void MyViewer::buildStreamLine(Vector v) {
	double w[3];
	MyMesh::VertexHandle v_array[3];
	MyMesh::FaceHandle f = getFace(v, w, v_array);

	//singularities.push_back(Singularity(v, f));
	double step = 0.5;
	Vector U = Vector(0, 0, 0);
	for (size_t i = 0; i < 3; i++)
	{
		U += w[i] * mesh.data(v_array[i]).u;
	}
	U = U.normalize();
	Vector cross[4];
	cross[0] = (calculateCrossFromU(U)[0]).normalize();
	cross[1] = (calculateCrossFromU(U)[1]).normalize();
	cross[2] = cross[0] * -1;
	cross[3] = cross[1] * -1;
	for (int i = 0; i < 4; i++) {
		std::vector<Vector> separatrice;
		separatrice.push_back(v);

		Vector prevDir = cross[i];
		Vector vtemp = v + cross[i] * step;
		separatrice.push_back(vtemp);
		followStreamLine(vtemp, cross[i], &separatrice,-1);
		double wtemp[3];
		MyMesh::VertexHandle v_arraytemp[3];

		separatrices2.push_back(separatrice);
	}

}
MyViewer::MyMesh::FaceHandle MyViewer::getFace(Vector v, double* w, MyMesh::VertexHandle* v_array) {
	double* A_data = new double[3 * 3];
	double* b_data = new double[3];

	gsl_vector* x = gsl_vector_alloc(3);
	gsl_permutation* perm = gsl_permutation_alloc(3);
	for (auto f : mesh.faces()) {
		mesh.data(f).hasSingularity = false;
		MyMesh::VertexHandle v1;
		MyMesh::VertexHandle v2;
		MyMesh::VertexHandle v3;
		int i = 0;
		for (auto v : mesh.fv_range(f))
		{
			switch (i)
			{
			case 0:
				v1 = v;
				i++;
				break;
			case 1:
				v2 = v;
				i++;
				break;
			case 2:
				v3 = v;
				i++;
				break;
			}
		}
		if (i < 3)
			throw;
		//Row 0
		A_data[0] = 1;
		A_data[1] = 1;
		A_data[2] = 1;
		b_data[0] = 1;

		//Row 1
		A_data[3] = mesh.point(v1)[0];
		A_data[4] = mesh.point(v2)[0];
		A_data[5] = mesh.point(v3)[0];
		b_data[1] = v[0];

		//Row 2
		A_data[6] = mesh.point(v1)[1];
		A_data[7] = mesh.point(v2)[1];
		A_data[8] = mesh.point(v3)[1];
		b_data[2] = v[1];
		auto A = gsl_matrix_view_array(A_data, 3, 3);
		auto b = gsl_vector_view_array(b_data, 3);


		// Decompose A into the LU form:
		int signum;       // Sign of the permutation
		gsl_linalg_LU_decomp(&A.matrix, perm, &signum);
		double det = gsl_linalg_LU_det(&A.matrix, signum);
		//qDebug() << det;
		if (fabs(det) > 0.00001) {
			gsl_linalg_LU_solve(&A.matrix, perm, &b.vector, x);

			if (x->data[0] > -0.0001 && x->data[1] > -0.0001 && x->data[2] > -0.0001) {
				w[0] = x->data[0];
				w[1] = x->data[1];
				w[2] = x->data[2];
				v_array[0] = v1;
				v_array[1] = v2;
				v_array[2] = v3;
				gsl_vector_free(x);
				gsl_permutation_free(perm);
				delete[] A_data;
				delete[] b_data;
				return f;
			}
		}



	}
	gsl_vector_free(x);
	gsl_permutation_free(perm);
	delete[] A_data;
	delete[] b_data;
	w[0] = -1;
	w[1] = -1;
	w[2] = -1;
	return mesh.faces().begin();
}

void MyViewer::followStreamLine(Vector v, Vector prevDir, std::vector<Vector>* streamline, int separatriceIdx, double step, double iterations) {
	double wtemp[3];
	MyMesh::VertexHandle v_arraytemp[3];
	Vector vtemp = v;
	MyMesh::FaceHandle prevFace;
	if (streamline->size() > 0)
		prevFace = getFace((*streamline)[0], wtemp, v_arraytemp);
	else
		prevFace = getFace(v, wtemp, v_arraytemp);
	int from = 0;
	for (size_t j = 0; j < iterations; j++)
	{
		MyMesh::FaceHandle ftemp = getFace(vtemp, wtemp, v_arraytemp);
		if (ftemp != prevFace) {
			if (from > 0)
				from--;
			mesh.data(prevFace).separaticeParts.push_back(SeparatricePart(separatriceIdx, from, streamline->size()));
			prevFace = ftemp;
			from = streamline->size()-1;
		}
		if (wtemp[0] < -0.5) {
			MyMesh::VertexHandle v1Boundary;
			MyMesh::VertexHandle v2Boundary;
			MyMesh::VertexHandle v3Boundary;
			bool v1BoundaryFound = false;
			bool v2BoundaryFound = false;
			bool v3BoundaryFound = false;
			for (auto v : v_arraytemp) {
				if (mesh.is_boundary(v)) {
					if (v1BoundaryFound && v2BoundaryFound) {
						v3BoundaryFound = true;
						v3Boundary = v;
					}
					else {
						if (v1BoundaryFound && !v2BoundaryFound) {
							v2Boundary = v;
							v2BoundaryFound = true;
						}
						else {
							v1Boundary = v;
							v1BoundaryFound = true;
						}
					}
				}
			}
			Vector P;
			if (v3BoundaryFound) {
				if (doIntersect(mesh.point(v1Boundary), mesh.point(v2Boundary), (*streamline)[streamline->size() - 2], (*streamline)[streamline->size() - 1], &P)) {
					std::vector<SeparatricePart> separatrices;
					separatrices.push_back(SeparatricePart(separatriceIdx, streamline->size() - 2, streamline->size() - 1, 1));
					corners.push_back(Corner(P, separatrices, true, v1Boundary, v2Boundary));
				} else if (doIntersect(mesh.point(v1Boundary), mesh.point(v3Boundary), (*streamline)[streamline->size() - 2], (*streamline)[streamline->size() - 1], &P)) {
					std::vector<SeparatricePart> separatrices;
					separatrices.push_back(SeparatricePart(separatriceIdx, streamline->size() - 2, streamline->size() - 1, 1));
					corners.push_back(Corner(P, separatrices, true, v1Boundary, v3Boundary));
				}
				else if (doIntersect(mesh.point(v2Boundary), mesh.point(v3Boundary), (*streamline)[streamline->size() - 2], (*streamline)[streamline->size() - 1], &P)) {
					std::vector<SeparatricePart> separatrices;
					separatrices.push_back(SeparatricePart(separatriceIdx, streamline->size() - 2, streamline->size() - 1, 1));
					corners.push_back(Corner(P, separatrices, true, v2Boundary, v3Boundary));
				}
			}
			else if (v2BoundaryFound) {
				if (doIntersect(mesh.point(v1Boundary),mesh.point(v2Boundary), (*streamline)[streamline->size() - 2], (*streamline)[streamline->size() - 1], &P)) {
					std::vector<SeparatricePart> separatrices;
					separatrices.push_back(SeparatricePart(separatriceIdx, streamline->size() - 2, streamline->size() - 1,1));
					corners.push_back(Corner(P, separatrices, true, v1Boundary, v2Boundary));
				}
			}
			else {
				//qDebug() << "StreamLine: (" << (*streamline)[streamline->size() - 1][0] << ", " << (*streamline)[streamline->size() - 1][1] << ") Boundary: " << mesh.point(v1Boundary)[0] << ", " << mesh.point(v1Boundary)[1] << ")";
				if (getDistance(mesh.point(v1Boundary), (*streamline)[streamline->size() - 1]) < 0.1) {
					std::vector<SeparatricePart> separatrices;
					separatrices.push_back(SeparatricePart(separatriceIdx, streamline->size() - 2, streamline->size() - 1,1));
					corners.push_back(Corner(mesh.point(v1Boundary), separatrices, true, v1Boundary, v1Boundary));
				}
			}
			break;
		}

		Vector Utemp = Vector(0, 0, 0);
		for (size_t i = 0; i < 3; i++)
		{
			Utemp += wtemp[i] * mesh.data(v_arraytemp[i]).u;
		}
		Utemp = Utemp.normalize();
		Vector crosstemp[4];
		crosstemp[0] = (calculateCrossFromU(Utemp)[0]).normalize();
		crosstemp[1] = (calculateCrossFromU(Utemp)[1]).normalize();
		crosstemp[2] = crosstemp[0] * -1;
		crosstemp[3] = crosstemp[1] * -1;
		prevDir = prevDir.normalize();
		double minDif = 99999;
		int minIdx = 0;
		for (size_t j = 0; j < 4; j++)
		{
			double temp = dot(crosstemp[j], prevDir);
			if (1 - temp < minDif) {
				minDif = 1 - temp;
				minIdx = j;
			}
		}
		vtemp = vtemp + crosstemp[minIdx] * step;
		streamline->push_back(vtemp);
		prevDir = crosstemp[minIdx];
	}

}
void MyViewer::findPartitionCorners() {
	for (auto f : mesh.faces()) {
		if (mesh.data(f).separaticeParts.size() > 1) {
			//mesh.data(f).tagged2 = true;
			for (size_t i = 0; i < mesh.data(f).separaticeParts.size(); i++)
			{
				int separatriceIdx = mesh.data(f).separaticeParts[i].separatriceIdx;
				for (size_t j = 0; j < mesh.data(f).separaticeParts.size(); j++)
				{
					int separatriceIdx2 = mesh.data(f).separaticeParts[j].separatriceIdx;
					if (i == j || (separatriceIdx == separatriceIdx2 && !(mesh.data(f).separaticeParts[i].fromI > mesh.data(f).separaticeParts[j].toI 
																		|| mesh.data(f).separaticeParts[j].fromI > mesh.data(f).separaticeParts[i].toI)))
						continue;
					for (size_t s1Idx = mesh.data(f).separaticeParts[i].fromI; s1Idx < mesh.data(f).separaticeParts[i].toI; s1Idx++)
					{
						for (size_t s2Idx = mesh.data(f).separaticeParts[j].fromI; s2Idx < mesh.data(f).separaticeParts[j].toI; s2Idx++)
						{
							Vector P;
							if (doIntersect(separatrices[separatriceIdx][s1Idx], separatrices[separatriceIdx][s1Idx + 1], separatrices[separatriceIdx2][s2Idx], separatrices[separatriceIdx2][s2Idx+1], &P)) {	
								std::vector<SeparatricePart> separatricesTemp;
								double t1 = (separatrices[separatriceIdx][s1Idx][0] - P[0]) / (separatrices[separatriceIdx][s1Idx][0] - separatrices[separatriceIdx][s1Idx + 1][0]);
								double t2 = (separatrices[separatriceIdx2][s2Idx][0] - P[0]) / (separatrices[separatriceIdx2][s2Idx][0] - separatrices[separatriceIdx2][s2Idx + 1][0]);
								separatricesTemp.push_back(SeparatricePart(separatriceIdx, s1Idx, s1Idx + 1,t1));
								separatricesTemp.push_back(SeparatricePart(separatriceIdx2, s2Idx, s2Idx + 1, t2));
								corners.push_back(Corner(P, separatricesTemp, false));
							}
						}
					}
				}
			}
		}
	}
	eliminateDoubleCorners();
}
bool MyViewer::doIntersect(Vector p1, Vector p2, Vector p3, Vector p4, Vector* P) {
	double t = ((p1[0] - p3[0]) * (p3[1] - p4[1]) - (p1[1] - p3[1]) * (p3[0] - p4[0]))
		/ ((p1[0] - p2[0]) * (p3[1] - p4[1]) - (p1[1] - p2[1]) * (p3[0] - p4[0]));
	double u = ((p2[0] - p1[0]) * (p1[1] - p3[1]) - (p2[1] - p1[1]) * (p1[0] - p3[0]))
		/ ((p1[0] - p2[0]) * (p3[1] - p4[1]) - (p1[1] - p2[1]) * (p3[0] - p4[0]));

	//qDebug() << "u: " << u << "t: " << t;
	if (t >= 0 && t <= 1 && u >= 0 && u <= 1) {
		(*P)[0] = p1[0] + t * (p2[0] - p1[0]);
		(*P)[1] = p1[1] + t * (p2[1] - p1[1]);
		(*P)[2] = 0;
		return true;
	}
	return false;
}
void MyViewer::eliminateDoubleCorners() {
	std::vector<Corner> tempCorners;
	for (int i = 0; i < corners.size(); i++) {
		for (auto c2 : corners) {
			if (getDistance(corners[i].pos, c2.pos) < 0.001){
				bool sameSp = false;
				for (auto sp2 : c2.separatrices)
				{
					for (auto sp : corners[i].separatrices)
					{
						if (sp == sp2)
						{
							sameSp = true;
							break;
						}
					}
					if (!sameSp) {
						corners[i].separatrices.push_back(sp2);
					}
				}
			}
		}

	}
	for (auto c : corners) {
		bool same = false;
		for (auto c2 : tempCorners) {
			if (getDistance(c.pos, c2.pos) < 0.001) {
				same = true;
				break;
			}
		}
		if (!same) {
			tempCorners.push_back(c);
		}
	}
	corners = tempCorners;
	std::vector<int> cornersToDelete;
	for (size_t cornerIdx = 0; cornerIdx < corners.size(); cornerIdx++)
	{
		if (corners[cornerIdx].separatrices.size() == 0 && !corners[cornerIdx].noSeparatrice)
			cornersToDelete.push_back(cornerIdx);
	}
	std::sort(cornersToDelete.begin(), cornersToDelete.end(), std::greater<int>());
	cornersToDelete.erase(std::unique(cornersToDelete.begin(), cornersToDelete.end()), cornersToDelete.end());
	for (auto cd : cornersToDelete)
		corners.erase(corners.begin() + cd);
	
	/*qDebug() << "After ----------------------------------------------";
	int idx = 0;
	for (auto c : corners) {
		qDebug() <<idx << " - pos: (" << c.pos[0] << ", " << c.pos[1] << ", " << c.pos[2]
			<< "), Boundary: " << c.boundary << ", BoundaryV1: " << c.boundaryV1.idx() << ", BoundaryV2: " << c.boundaryV2.idx();
		for (auto sp : c.separatrices) {
			qDebug() << "\t" << "sp index: " << sp.separatriceIdx << "fromI: " << sp.fromI << "toI: " << sp.toI;
		}
		idx++;
	}*/
}
double MyViewer::getDistance(Vector p1, Vector p2) {
	return (p1 - p2).length();
}

void MyViewer::detectRegions() {
	addCornersToEdgesAndVertices();
	std::vector<std::vector<int>> neighbours = getCornerNeighbours();
	int idx = 0;
	//for (auto c : corners) {
	//	qDebug() << idx << " - pos: (" << c.pos[0] << ", " << c.pos[1] << ", " << c.pos[2]
	//		<< "), Boundary: " << c.boundary << ", BoundaryV1: " << c.boundaryV1.idx() << ", BoundaryV2: " << c.boundaryV2.idx();
	//	for (auto sp : c.separatrices) {
	//		qDebug() << "\t" << "sp index: " << sp.separatriceIdx << "fromI: " << sp.fromI << "toI: " << sp.toI << " t: " << sp.t;
	//	}
	//	for (auto n : neighbours[idx])
	//		qDebug() << "\t n: " << n;
	//	idx++;
	//}
	idx = 0;
	int count = 0;
	for (auto c : corners) {
		//qDebug() << idx;
		for (auto n : neighbours[idx]) {
			std::vector<Vector> points;
			std::vector<int> pointsCorner;

			points.push_back(c.pos);
			pointsCorner.push_back(idx);
			points.push_back(corners[n].pos);
			pointsCorner.push_back(n);

			//qDebug() << "Pushed Back" << idx;
			//qDebug() << "Pushed Back" << n;
			int prevCorner = n;
			int prevCorner2 = idx;
			bool returned = true;
			int debugCount = 0;
			while (true) {
				double maxDot = -9999999;
				int maxCorner = -1;
				for (auto n2 : neighbours[prevCorner]) {
					if (n2 == prevCorner2)
						continue;
					if (!isLeft(points[points.size() - 2], points[points.size() - 1], corners[n2].pos)) {
						//qDebug() << "Right" << n2;
						Vector v1 = (points[points.size() - 1] - points[points.size() - 2]).normalize();
						Vector v2 = (points[points.size() - 1]- corners[n2].pos ).normalize();
						//qDebug() << "Dot: " << dot(v1, v2) <<"  MaxDot: " <<maxDot;
						if (dot(v1, v2) > maxDot) {
						//	qDebug() << " New MaxDot0: " << dot(v1, v2);
							maxDot = dot(v1, v2);
							maxCorner = n2;
						//	qDebug() << " New MaxDot1: " << maxDot;

						}
					}
				}
				if (maxCorner == -1)
				{
					returned = false;
					break;
				}
				else if (maxCorner != idx) {
					pointsCorner.push_back(maxCorner);
					points.push_back(corners[maxCorner].pos);
					prevCorner2 = prevCorner;
					prevCorner = maxCorner;
				}
				else
					break;
				if (debugCount++ > 10) {
					qDebug() << "NOT GOOD";
					return;
					break;
				}
			}
			if (returned) {
				bool alreadyIn = false;
				for (auto region : regions)
				{
					int same = 0;
					for (auto p : region.corners) {
						for (int i = 0; i < points.size(); i++) {
							if (pointsCorner[i] == p) {
								same++;
							}
						}
						if (same == points.size()) {
							alreadyIn = true;
							break;
						}
					}
				}
				if (!alreadyIn) {
					regions.push_back(Region(pointsCorner));
				}
				count++;
			}
			//qDebug() << "##############################";
		}
		idx++;
	}
	//updateRegionSides();
}

std::vector<std::vector<int>> MyViewer::getCornerNeighbours() {
	std::vector<std::vector<int>> neighbours;
	std::vector<int> boundaryCornerChain = buildBoundaryCornerChain();
	int idx = 0;
	for (auto c : corners) {
		std::vector<int> tempNeighbours;
		for (auto sp : c.separatrices) {
			int minDist1 = 9999999;
			double minDistT1 = 9999;
			int minCornerIdx1 = -1;
			int minDist2 = 9999999;
			double minDistT2 = 9999;
			int minCornerIdx2 = -1;
			int idx2 = -1;
			for (auto c2 : corners) {
				idx2++;
				if (getDistance(c.pos, c2.pos) < 0.01)
					continue;
				for (auto sp2 : c2.separatrices) {
					if (sp.separatriceIdx == sp2.separatriceIdx) {

						int tempIdx = sp.fromI - sp2.fromI;
						double tempT = sp.t - sp2.t;
						if (tempIdx > 0) {
							if (tempIdx < minDist1) {
								minDist1 = tempIdx;
								minDistT1 =  1 - sp2.t;
								minCornerIdx1 = idx2;
							}
							else if (tempIdx == minDist1) {
								if (1 - sp2.t < minDistT1) {
									minDist1 = tempIdx;
									minDistT1 = 1 - sp2.t;
									minCornerIdx1 = idx2;
								}
							}
						}
						else if (tempIdx < 0 ){
							tempIdx = -1 * tempIdx;
							if (tempIdx < minDist2) {
								minDist2 = tempIdx;
								minDistT2 = sp2.t;
								minCornerIdx2 = idx2;
							}
							else if (tempIdx == minDist2) {
								if (sp2.t < minDistT2) {
									minDist2 = tempIdx;
									minDistT2 = sp2.t;
									minCornerIdx2 = idx2;
								}
							}
						}
						else {
							if (tempT > 0) {
								if ((minDist1 == 0 && tempT < minDistT1) || minDist1 > 0) {
									minDist1 = 0;
									minDistT1 = tempT;
									minCornerIdx1 = idx2;
								}
							}
							else {
								if ((minDist2 == 0 && sp2.t < minDistT2) || minDist2 > 0) {
									minDist2 = 0;
									minDistT2 = sp2.t;
									minCornerIdx2 = idx2;
								}
							}
						}

					}
				}
			}
			if (minCornerIdx1 > -1)
				tempNeighbours.push_back(minCornerIdx1);
			if (minCornerIdx2 > -1)
				tempNeighbours.push_back(minCornerIdx2);
		}
		if (c.boundary) {
			std::vector<int> neighboursBorder = getBoundaryNeighbourCorner(idx, boundaryCornerChain);
			tempNeighbours.push_back(neighboursBorder[0]);
			tempNeighbours.push_back(neighboursBorder[1]);
		}
		neighbours.push_back(tempNeighbours);
		idx++;
	}
	return neighbours;
}
void MyViewer::addCornersToEdgesAndVertices() {
	int idx = 0;
	for (auto c : corners) {
		if (c.boundary) {
			if (c.boundaryV1 != c.boundaryV2) {
				mesh.data(c.boundaryV1).corners.push_back(idx);
				mesh.data(c.boundaryV2).corners.push_back(idx);
				MyMesh::EdgeHandle eh = getCommonEdge(c.boundaryV1, c.boundaryV2);
				MyMesh::VertexHandle v = mesh.from_vertex_handle(mesh.halfedge_handle(eh, 0));
				double dist = getDistance(c.pos, mesh.point(v));
				if (mesh.data(eh).corners.size() > 0) {
					bool added = false;
					for (size_t i = 0; i < mesh.data(eh).corners.size(); i++)
					{
						if (getDistance(corners[mesh.data(eh).corners[i]].pos, mesh.point(v)) > dist) {
							auto itPos = mesh.data(eh).corners.begin() + i;
							mesh.data(eh).corners.insert(itPos, idx);
							added = true;
							break;
						}
					}
					if (!added) {
						mesh.data(eh).corners.push_back(idx);
					}
				}
				else {
					mesh.data(eh).corners.push_back(idx);
				}
			}
			else {
				mesh.data(c.boundaryV1).corners.push_back(idx);
				for (MyMesh::VertexEdgeIter ve_iter = mesh.ve_iter(c.boundaryV1); ve_iter.is_valid(); ve_iter++) {
					if (mesh.is_boundary(ve_iter.handle())) {
						MyMesh::VertexHandle v = mesh.from_vertex_handle(mesh.halfedge_handle(ve_iter, 0));
						double dist = getDistance(c.pos, mesh.point(v));
						if (mesh.data(ve_iter.handle()).corners.size() > 0) {
							bool added = false;
							for (size_t i = 0; i < mesh.data(ve_iter.handle()).corners.size(); i++)
							{
								if (getDistance(corners[mesh.data(ve_iter.handle()).corners[i]].pos, mesh.point(v)) > dist) {
									auto itPos = mesh.data(ve_iter.handle()).corners.begin() + i;
									mesh.data(ve_iter.handle()).corners.insert(itPos, idx);
									added = true;
									break;
								}
							}
							if (!added) {
								mesh.data(ve_iter.handle()).corners.push_back(idx);
							}
						}
						else {
							mesh.data(ve_iter.handle()).corners.push_back(idx);
						}
					}
				}
			}
		}
		idx++;
	}
}
MyViewer::MyMesh::VertexHandle MyViewer::getNextBoundary(MyMesh::VertexHandle v, MyMesh::VertexHandle vPrev) {
	for (MyMesh::VertexEdgeIter ve_iter = mesh.ve_iter(v); ve_iter.is_valid(); ve_iter++) {
		if (mesh.is_boundary(ve_iter.handle())) {
			MyMesh::HalfedgeHandle  he = mesh.halfedge_handle(ve_iter.handle(), 0);
			if (mesh.from_vertex_handle(he) == v && mesh.to_vertex_handle(he) != vPrev)
				return mesh.to_vertex_handle(he);
			if (mesh.to_vertex_handle(he) == v && mesh.from_vertex_handle(he) != vPrev)
				return mesh.from_vertex_handle(he);
		}
	}
	throw;
}

int MyViewer::findNextBoundaryCorner(MyMesh::VertexHandle boundaryV1, MyMesh::VertexHandle boundaryV2) {
	MyMesh::VertexHandle prevV = boundaryV1;
	MyMesh::VertexHandle v = getNextBoundary(boundaryV1, boundaryV2);
	while (mesh.data(v).corners.size() == 0) {
		MyMesh::VertexHandle temp = v;
		v = getNextBoundary(v, prevV);
		prevV = temp;
	}
	if (mesh.data(v).corners.size() > 0) {
		double minDist = 9999;
		int minCornerIdx;
		for (auto c2 : mesh.data(v).corners)
		{
			if (getDistance(mesh.point(v), corners[c2].pos) < minDist) {
				minDist = getDistance(mesh.point(v), corners[c2].pos);
				minCornerIdx = c2;
			}
		}
		return minCornerIdx;
	}
	return -1;
}
bool MyViewer::isLeft(Vector a, Vector b, Vector c) {
	return ((b[0]- a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])) > 0;
}

std::vector<int> MyViewer::buildBoundaryCornerChain() {
	std::vector<int> chain;
	for (auto v: mesh.vertices())
	{
		if (mesh.is_boundary(v)) {
			MyMesh::VertexHandle v1 = v;
			MyMesh::VertexHandle v2 = getNextBoundary(v, v);
			int count = 0;
			while (v2 != v) {
				MyMesh::EdgeHandle e = getCommonEdge(v1, v2);
				if (mesh.data(e).corners.size() > 0) {
					if (getDistance(corners[mesh.data(e).corners[0]].pos, mesh.point(v1)) < getDistance(corners[mesh.data(e).corners[mesh.data(e).corners.size() - 1]].pos, mesh.point(v1))) {
						for (size_t i = 0; i < mesh.data(e).corners.size(); i++)
						{
							if (!(chain.size() > 0 && chain[chain.size() - 1] == mesh.data(e).corners[i])) {
								chain.push_back(mesh.data(e).corners[i]);
							}
						}
					}
					else {
						for (int i = mesh.data(e).corners.size() - 1; i >= 0; i--)
						{
							if (!(chain.size() > 0 && chain[chain.size() - 1] == mesh.data(e).corners[i])) {
								chain.push_back(mesh.data(e).corners[i]);

							}
						}
					}
				}

				MyMesh::VertexHandle tempVh = v2;
				v2 = getNextBoundary(v2, v1);
				v1 = tempVh;

			}
			return chain;
		}
	}
	return chain;
}
std::vector<int> MyViewer::getBoundaryNeighbourCorner(int corner, std::vector<int> boundaryCornerChain) {
	int cornerIdx = 0;
	for (size_t i = 0; i < boundaryCornerChain.size(); i++)
	{
		if (boundaryCornerChain[i] == corner) {
			cornerIdx = i;
			break;
		}
	}
	int prevIdx = (cornerIdx > 0) ? cornerIdx- 1 : boundaryCornerChain.size()-1;
	int nextIdx = (cornerIdx < boundaryCornerChain.size()-1)? cornerIdx + 1 : 0;
	std::vector<int> returnVec;
	returnVec.push_back(boundaryCornerChain[prevIdx]);
	returnVec.push_back(boundaryCornerChain[nextIdx]);
	return returnVec;
}
void MyViewer::updateRegionSides() {
	for (size_t i = 0; i < regions.size(); i++)
	{
		std::vector<int> tempIDxs = {3,0, 0,1, 1,2, 2,3, 3,0 };
		for (size_t j = 0; j < tempIDxs.size(); j+=2)
		{
			int j1 = tempIDxs[j];
			int j2 = tempIDxs[j+1];

		}
	}
}
void MyViewer::collapseRegions() {
	std::vector<regionPair> regionPairs = getRegionPairs();
	
	double epsilon = 0.1;
	double minAlpha = 99999999999;
	double minSmallerArea = 99999999999;
	int minIdx = -1;
	int idx = 0;
	for (auto rp : regionPairs) {
		
		idx++;
	}
}
std::vector<MyViewer::regionPair> MyViewer::getRegionPairs() {
	std::vector<regionPair> regionPairs;
	int idx = 0;
	for (auto region : regions) {
		int idx2 = 0;
		//qDebug() << idx << "-region";
		//for (auto c : region.corners) {
		//	qDebug() << "\tcorner: " << c << " pos: (" << corners[c].pos[0] << ", " << corners[c].pos[1] << ")";
		//}
		for (auto region2 : regions) {

			if (idx >= idx2) {
				idx2++;
				continue;
			}
			std::vector<int> tempIDxs = { 0,1, 1,2, 2,3, 3,0 };
			int tempi = -1;
			int tempi2 = -1;
			int tempj = -1;
			int tempj2 = -1;
			for (size_t i = 0; i < tempIDxs.size(); i += 2)
			{
				int i1 = tempIDxs[i];
				int i2 = tempIDxs[i + 1];
				for (size_t j = 0; j < tempIDxs.size(); j += 2)
				{
					int j1 = tempIDxs[j];
					int j2 = tempIDxs[j + 1];
					if (region.corners[i1] == region2.corners[j1] && region.corners[i2] == region2.corners[j2]) {
						tempi = i1;
						tempi2 = i2;
						tempj = j1;
						tempj2 = j2;
					}
					if (region.corners[i1] == region2.corners[j2] && region.corners[i2] == region2.corners[j1]) {
						tempi = i1;
						tempi2 = i2;
						tempj = j2;
						tempj2 = j1;
					}
				}
			}
			if (tempi != -1) {
				if (!corners[region.corners[tempi]].boundary && corners[region.corners[tempi]].separatrices.size() <= 3 && !corners[region.corners[tempi2]].boundary && corners[region.corners[tempi2]].separatrices.size() <= 3)
					continue;
				Vector i1V = getSideVector(region.corners, tempi, tempi2);
				Vector i2V = getSideVector(region.corners, tempi2, tempi);
				Vector j1V = getSideVector(region2.corners, tempj, tempj2);
				Vector j2V = getSideVector(region2.corners, tempj2, tempj);

				float ratio = ((i1V.length() + i2V.length()) / 2) / (corners[region.corners[tempi]].pos - corners[region.corners[tempi2]].pos).length();
				ratio += ((j1V.length() + j2V.length()) / 2) / (corners[region2.corners[tempj]].pos - corners[region2.corners[tempj2]].pos).length();

				i1V = i1V.normalize();
				i2V = i2V.normalize();
				j1V = j1V.normalize();
				j2V = j2V.normalize();

				//qDebug() << "i1V: (" << i1V[0] << ", " <<  i1V[1] << ") i2V: (" << i2V[0] << ", " << i2V[1] << ")";
				//qDebug() << "j1V: (" << j1V[0] << ", " <<  j1V[1] << ") j2V: (" << j2V[0] << ", " << j2V[1] << ")";
				//qDebug() << "--------------------------------";
				//qDebug() << idx << "-region";
				//for (auto c : region.corners) {
				//	qDebug() << "\tcorner: " << c << " pos: (" << corners[c].pos[0] << ", " << corners[c].pos[1] << ")";
				//}
				//qDebug() << idx2 << "-region";
				//for (auto c : region2.corners) {
				//	qDebug() << "\tcorner: " << c << " pos: (" << corners[c].pos[0] << ", " << corners[c].pos[1] << ")";
				//}
				//qDebug() << "--------------------------------";
				double area1 = calculateArea(corners[region.corners[0]].pos, corners[region.corners[1]].pos, corners[region.corners[2]].pos, corners[region.corners[3]].pos);
				double area2 = calculateArea(corners[region2.corners[0]].pos, corners[region2.corners[1]].pos, corners[region2.corners[2]].pos, corners[region2.corners[3]].pos);
				float alpha = abs(-1 - dot(i1V, j1V)) + abs(-1 - dot(i2V, j2V));
				regionPairs.push_back(regionPair(idx, idx2, alpha, ratio, area1 < area2 ? area1 : area2));
			}
			idx2++;
		}
		idx++;
	}
	for (auto rp : regionPairs) {
		qDebug() << "regionIdx1: " << rp.regionIdx1 << "regionIdx2: " << rp.regionIdx2 << "alpha: " << rp.alpha << " ratio: " << rp.ratio << " area: " << rp.smallerArea;
	}
	return regionPairs;
}
MyViewer::Vector MyViewer::getSideVector(std::vector<int> tmpcorners, int i1, int i2) {
	if (i2 > i1) {
		if (i2 != 3 || i1 != 0) {
			int tempi3 = i1 - 1;
			tempi3 = tempi3 < 0 ? 3 : tempi3;
			return corners[tmpcorners[tempi3]].pos - corners[tmpcorners[i1]].pos;
		}
		else {
			int tempi3 = 1;
			return corners[tmpcorners[tempi3]].pos - corners[tmpcorners[i1]].pos;
		}
	}
	else {
		if (i1 != 3 || i2 != 0) {
			int tempi3 = i1 + 1;
			tempi3 = tempi3 > 3 ? 0 : tempi3;
			return corners[tmpcorners[tempi3]].pos - corners[tmpcorners[i1]].pos;
		}
		else {
			int tempi3 = 2;
			return corners[tmpcorners[tempi3]].pos - corners[tmpcorners[i1]].pos;
		}
	}
}
double MyViewer::calculateArea(Vector A, Vector B, Vector C, Vector D) {
	double a = (A - B).length();
	double b = (B - C).length();
	double c = (C - D).length();
	double d = (D - A).length();
	double s = (a + b + c + d) / 2;

	double alpha = acos(dot((B-A).normalize(), (D-A).normalize()));
	double gamma = acos(dot((B-C).normalize(), (D-C).normalize()));

	double T = sqrt((s - a) * (s - b) * (s - c) * (s - d) - a * b * c * d * cos((alpha + gamma) / 2) * cos((alpha + gamma) / 2));
	return T;
}