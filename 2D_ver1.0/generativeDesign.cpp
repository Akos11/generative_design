#include "MyViewer.h"
#include <OpenMesh/Tools/Smoother/JacobiLaplaceSmootherT.hh>

//Generative Design
void MyViewer::calculateIncidence() {
	//mesh.reset_status();
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

void MyViewer::reMeshOrganicBoundaries(double L, int iterations) {
	
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

void MyViewer::reMeshEdgeLength(double L) {
	{
		//Split
		//bool hasSplit = true;
		//while (hasSplit) {
		//	hasSplit = false;
		//	for (MyMesh::EdgeIter it = mesh.edges_begin(); it != mesh.edges_end(); ++it) {
		//		if (mesh.calc_edge_length(*it) > 4.0 / 5.0 * L) {
		//			const MyViewer::MyTraits::Point& newPoint = 
		//				(mesh.point(mesh.from_vertex_handle(mesh.halfedge_handle(*it, 0))) 
		//					+ 
		//				 mesh.point(mesh.to_vertex_handle(mesh.halfedge_handle(*it, 0))))
		//				/2;

		//			mesh.split(*it, newPoint);
		//			hasSplit = true;
		//			break;
		//		}
		//	}
		//	mesh.garbage_collection();
		//}

		////Collapse
		//bool hasCollapsed = true;
		//while (hasCollapsed /*&& mesh.n_vertices() > 269*/) {
		//	hasCollapsed = false;
		//	for (MyMesh::HalfedgeIter it = mesh.halfedges_begin(); it != mesh.halfedges_end(); ++it) {

		//		if (mesh.calc_edge_length(*it) < 4.0 / 3.0 * L && !mesh.status(mesh.edge_handle(*it)).deleted()) {
		//			//qDebug() << " Before Vertex: " << mesh.to_vertex_handle(*it).idx();
		//			//for (auto h : mesh.vih_range(mesh.to_vertex_handle(*it))) {
		//			//	qDebug() << "FROM vertex: " << mesh.from_vertex_handle(h).idx() << " To Vertex: " << mesh.to_vertex_handle(h).idx();
		//			//	qDebug() << h.idx();
		//			//}
		//			if (mesh.to_vertex_handle(*it).idx() == 146) {
		//				MyMesh::HalfedgeHandle  h = *it;
		//				MyMesh::HalfedgeHandle  hn = mesh.next_halfedge_handle(h);
		//				MyMesh::HalfedgeHandle  hp = mesh.prev_halfedge_handle(h);

		//				MyMesh::HalfedgeHandle  o = mesh.opposite_halfedge_handle(h);
		//				MyMesh::HalfedgeHandle  on = mesh.next_halfedge_handle(o);
		//				MyMesh::HalfedgeHandle  op = mesh.prev_halfedge_handle(o);

		//				MyMesh::FaceHandle      fh = mesh.face_handle(h);
		//				MyMesh::FaceHandle      fo = mesh.face_handle(o);

		//				MyMesh::VertexHandle    vh = mesh.to_vertex_handle(h);
		//				MyMesh::VertexHandle    vo = mesh.to_vertex_handle(o);

		//				qDebug() << "h " << h.idx() << "hn " << hn.idx() << "hp " << hp.idx() << "vh " << vh.idx() << "vo " << vo.idx();
		//				for (auto h : mesh.vih_range(mesh.to_vertex_handle(*it))) {
		//					qDebug() << "FROM vertex: " << mesh.from_vertex_handle(h).idx() << " To Vertex: " << mesh.to_vertex_handle(h).idx();
		//					qDebug() << h.idx();
		//				}
		//				qDebug() << "--------------------------------------TEST: " <<mesh.n_vertices();
		//			}
		//			MyMesh::VertexHandle fromVh = mesh.from_vertex_handle(*it);
		//			MyMesh::VertexHandle toVh = mesh.to_vertex_handle(*it);
		//			MyMesh::Point newPos = (mesh.point(fromVh) + mesh.point(toVh)) / 2.0f;
		//			mesh.set_point(fromVh, newPos);
		//			mesh.set_point(toVh, newPos);
		//			mesh.collapse(*it);
		//			hasCollapsed = true;
		//			mesh.garbage_collection();
		//			for (auto h : mesh.vih_range(mesh.to_vertex_handle(*it))) {
		//				qDebug() << "FROM vertex: " << mesh.from_vertex_handle(h).idx() << " To Vertex: " << mesh.to_vertex_handle(h).idx();
		//				qDebug() << h.idx();
		//			}
		//			//qDebug() << " Vertex: " << mesh.to_vertex_handle(*it).idx();
		//			//for (auto h : mesh.vih_range(mesh.to_vertex_handle(*it))) {
		//			//	qDebug() << "FROM vertex: " << mesh.from_vertex_handle(h).idx() << " To Vertex: " << mesh.to_vertex_handle(h).idx();
		//			//	qDebug() << h.idx();
		//			//}
		//			qDebug() << mesh.n_vertices();

		//			
		//			
		//			break;
		//		}
		//		
		//	}
		//	

		//	mesh.garbage_collection();
		//}
		////mesh.garbage_collection();
		//for (auto v : mesh.vertices()) {
		//	qDebug() << " Vertex: " << v.idx();
		//	for (auto h : mesh.vih_range(v)) {
		//		qDebug() << "FROM vertex: " << mesh.from_vertex_handle(h).idx() << " To Vertex: " << mesh.to_vertex_handle(h).idx();
		//		qDebug() << h.idx();
		//	}
		//}
	}
	bool split = true;
	bool collapse = true;
	while (split) {
		split = false;
		for (MyMesh::EdgeIter it = mesh.edges_begin(); it != mesh.edges_end(); ++it) {
			if (!checkIfEdgeTagged(it.handle()))
				continue;

			double length = mesh.calc_edge_length(*it);

			if (length > 4.0 / 3.0 * L) {
				//qDebug() << "split";
				const MyViewer::MyTraits::Point& newPoint =
					(mesh.point(mesh.from_vertex_handle(mesh.halfedge_handle(*it, 0)))
						+
						mesh.point(mesh.to_vertex_handle(mesh.halfedge_handle(*it, 0))))
					/ 2;
				MyMesh::VertexHandle tempVh = mesh.split_copy(*it, newPoint);
				mesh.status(tempVh).set_tagged(true);
				split = true;
				mesh.garbage_collection();
				break;
			}
		}
	}
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

void MyViewer::reMeshTagvertices() {
	mesh.reset_status();
	//Tag vertices on the boundary of the 2 regions
	for (auto v : mesh.vertices()) {
		if (mesh.data(v).I.size() > 0) {
			for (MyMesh::VertexVertexIter vv_iter = mesh.vv_iter(v); vv_iter.is_valid(); vv_iter++) {
				if (mesh.data(*vv_iter).I.size() == 0) {
					mesh.status(v).set_tagged(true);
					break;
				}
			}

		}

	}
	//Tag the 2 ring neighbours of the previously tagged vertices
	for (size_t i = 0; i < 2; i++)
	{
		for (auto v : mesh.vertices()) {
			if (mesh.status(v).tagged()) {
				for (MyMesh::VertexVertexIter vv_iter = mesh.vv_iter(v); vv_iter.is_valid(); vv_iter++) {
					mesh.status(*vv_iter).set_tagged2(true);
				}
			}
		}
		for (auto v : mesh.vertices()) {
			if (mesh.status(v).tagged2()) {
				mesh.status(v).set_tagged2(false);
				mesh.status(v).set_tagged(true);
			}
		}
	}
}
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
				if (squaredDifferenceFlip < squaredDifference) {
					mesh.flip(*it);
					flipped = true;
					break;
				}
			}
		}
	}
}
void MyViewer::reMeshVertexPositions() {
	for (auto vh : mesh.vertices())
	{
		if (!mesh.status(vh).tagged() || mesh.is_boundary(vh))
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
	return (mesh.status(mesh.from_vertex_handle(hh)).tagged() && mesh.status(mesh.to_vertex_handle(hh)).tagged());
}
void MyViewer::reMeshSmoothing(int iterations) {
	for (size_t i = 0; i < iterations; i++)
	{
		reMeshSmoothingIteration();
	}
}
void MyViewer::reMeshSmoothingIteration() {
	for (auto v : mesh.vertices()) {
		mesh.status(v).set_tagged2(false);
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
				mesh.status(v).set_tagged2(true);
			}
			
		}
		nextSet = false;
		prevSet = false;

	}
	for (auto v : mesh.vertices()) {
		if (mesh.status(v).tagged2()) {
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