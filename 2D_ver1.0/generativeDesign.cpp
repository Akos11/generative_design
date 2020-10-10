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

void MyViewer::reMeshOrganicBoundaries() {

}