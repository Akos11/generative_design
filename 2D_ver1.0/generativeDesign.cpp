#include "MyViewer.h"

//Generative Design
void MyViewer::calculateIncidence() {
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
	for (size_t i = 0; i < iterations; i++)
	{
		reMeshEdgeLength(L);
		reMeshVertexValences();
		reMeshVertexPositions();
	}
	updateMesh();
	update();
}

void MyViewer::reMeshEdgeLength(double L) {
    mesh.reset_status();

	//Split
	bool hasSplit = true;
	while (hasSplit) {
		hasSplit = false;
		for (MyMesh::EdgeIter it = mesh.edges_begin(); it != mesh.edges_end(); ++it) {
			if (mesh.calc_edge_length(*it) > 4.0 / 5.0 * L) {
				const MyViewer::MyTraits::Point& newPoint = 
					(mesh.point(mesh.from_vertex_handle(mesh.halfedge_handle(*it, 0))) 
						+ 
					 mesh.point(mesh.to_vertex_handle(mesh.halfedge_handle(*it, 0))))
					/2;

				mesh.split(*it, newPoint);
				hasSplit = true;
				//break;
			}
		}
		mesh.garbage_collection();
	}


	qDebug() << 2;

	//Collapse
	bool hasCollapsed = true;
	while (hasCollapsed) {
		hasCollapsed = false;
		for (MyMesh::HalfedgeIter it = mesh.halfedges_begin(); it != mesh.halfedges_end(); ++it) {
			if (mesh.calc_edge_length(*it) < 4.0 / 3.0 * L) {
				mesh.collapse(*it);
				hasCollapsed = true;
				break;
			}
		}
		mesh.garbage_collection();
	}
	
}

void MyViewer::reMeshVertexValences() {
}
void MyViewer::reMeshVertexPositions() {
}

