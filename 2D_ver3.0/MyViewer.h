// -*- mode: c++ -*-
#pragma once
#include <string>

#include <QGLViewer/qglviewer.h>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include "ConstraintSolid.h"
#include "CircleSolid.h"
#include <OpenMesh/Tools/Subdivider/Uniform/CatmullClarkT.hh>
#include <queue>
#include <OpenMesh/Tools/Decimater/Observer.hh>
#include <OpenMesh/Tools/Decimater/DecimaterT.hh>
#include <OpenMesh/Tools/Decimater/ModQuadricT.hh>
#include <OpenMesh/Tools/Decimater/ModRoundnessT.hh>


using qglviewer::Vec;

class MyViewer : public QGLViewer {
	Q_OBJECT

public:
	explicit MyViewer(QWidget* parent);
	virtual ~MyViewer();

	inline double getCutoffRatio() const;
	inline void setCutoffRatio(double ratio);
	inline double getMeanMin() const;
	inline void setMeanMin(double min);
	inline double getMeanMax() const;
	inline void setMeanMax(double max);
	inline const double* getSlicingDir() const;
	inline void setSlicingDir(double x, double y, double z);
	inline double getSlicingScaling() const;
	inline void setSlicingScaling(double scaling);
	bool openMesh(const std::string& filename, bool update_view = true);
	bool openBezier(const std::string& filename, bool update_view = true);
	bool openGenerative(const std::string& filename, bool update_view = true);
	bool saveBezier(const std::string& filename);
signals:
	void startComputation(QString message);
	void midComputation(int percent);
	void endComputation();

protected:
	virtual void init() override;
	virtual void draw() override;
	virtual void drawWithNames() override;
	virtual void postSelection(const QPoint& p) override;
	virtual void keyPressEvent(QKeyEvent* e) override;
	virtual void mouseMoveEvent(QMouseEvent* e) override;
	virtual QString helpString() const override;

private:
	struct Flags {
		bool tagged;
		bool temporary_tagged;
		bool temporary_tagged2;
		bool locked;
	};
	struct Bfs {
		bool hasPrev;
		OpenMesh::VertexHandle prevVertex;
	};
	struct SeparatricePart {
		int separatriceIdx;
		int fromI;
		int toI;
		SeparatricePart(int separatriceIdx_, int fromI_, int toI_) : separatriceIdx(separatriceIdx_),fromI(fromI_), toI(toI_) {}
	};

	struct MyTraits : public OpenMesh::DefaultTraits {
		using Point = OpenMesh::Vec3d; // the default would be Vec3f
		using Normal = OpenMesh::Vec3d;
		VertexTraits{
		  double mean;              // approximated mean curvature
		  std::vector<int> I;
		  OpenMesh::Vec3d newPos;
		  Flags flags;
		  Bfs Bfs;
		  bool edgeTagged;
		  double angle;
		  int idx;
		  OpenMesh::Vec3d u;		//directionality
		};
		EdgeTraits{
			double squarness;
			bool tagged;
		};
		FaceTraits{
			bool tagged;
			bool tagged2;
			bool hasPrev;
			OpenMesh::FaceHandle prevFace;
			bool hasSingularity;
			std::vector<SeparatricePart> separaticeParts;
		};
	};
	using MyMesh = OpenMesh::TriMesh_ArrayKernelT<MyTraits>;
	using Vector = OpenMesh::VectorT<double, 3>;
	class generativeDesign;
	struct Corner {
		Vector pos;
		SeparatricePart separatrice1;
		SeparatricePart separatrice2;
		bool boundary;
		MyMesh::VertexHandle boundaryV1;
		MyMesh::VertexHandle boundaryV2;
		Corner(Vector pos_, SeparatricePart separatrice1_, SeparatricePart separatrice2_, bool boundary_, MyMesh::VertexHandle boundaryV1_ = MyMesh::VertexHandle(), MyMesh::VertexHandle boundaryV2_ = MyMesh::VertexHandle())
			: pos(pos_), separatrice1(separatrice1_), separatrice2(separatrice2_),boundaryV1(boundaryV1_), boundaryV2(boundaryV2_) {}
	};
	std::vector<Corner> corners;
	// Mesh
	void updateMesh(bool update_mean_range = true);
	void updateVertexNormals();
	void updateVertexNormalsConstraint();
	void updateVertexNormalsPartition();
	void localSystem(const Vector& normal, Vector& u, Vector& v);
	double voronoiWeight(MyMesh::HalfedgeHandle in_he);
	void updateMeanMinMax();
	void updateMeanCurvature(bool update_min_max = true);

	// Bezier
	static void bernsteinAll(size_t n, double u, std::vector<double>& coeff);
	void generateMesh();

	// Visualization
	void setupCamera();
	Vec meanMapColor(double d) const;
	void drawControlNet() const;
	void drawAxes() const;
	void drawAxesWithNames() const;
	static Vec intersectLines(const Vec& ap, const Vec& ad, const Vec& bp, const Vec& bd);

	// Other
	void fairMesh();

	//generative
		void resetFlags();
		void resetEdgeProps();
		void calculateIncidence();
		//Remesh
		void decimate();
		void reMeshTagvertices(bool all = false);
		void reMeshAll(int iterations);
		void reMeshOrganicBoundaries(double L, int iterations);
		void reMeshEdgeLength(double L);
		void reMeshVertexValences();
		void reMeshVertexPositions(bool remeshUnlockedBoundaries = false);
		void reMeshSmoothing(int iterations);
		void reMeshSmoothingIteration();
		bool checkIfEdgeTagged(MyMesh::EdgeHandle eh);
		bool checkIfEdgeTagged(MyMesh::HalfedgeHandle hh);
		bool checkIfBetweenRegions(MyMesh::VertexHandle vh);
		void removeIncidentVertices();
		OpenMesh::Subdivider::Uniform::CatmullClarkT<MyMesh> catmul_clark_subdivider;
		void catmullClark();
		void partition();
			int counter;
		bool hasCommonEdge(MyMesh::VertexHandle vh1, MyMesh::VertexHandle vh2);
		MyMesh::EdgeHandle getCommonEdge(MyMesh::VertexHandle vh1, MyMesh::VertexHandle vh2);
		//PDE based approach

		void initiatePDE();
		void initiateBoundaryConstraints(double angleThreshold);
		void drawPDEu();
		void drawPDEcross();
		void initVFunction(int iterations);
		void buildNeighbours(MyMesh::VertexHandle* omega, std::vector<int>* neighbours);
		void calculateGrad(MyMesh::VertexHandle* omega, std::vector<int>* neighbours, double* func, double* gradFunc);
		std::vector<int> getNeighbours(MyMesh::VertexHandle* omega, int idx);
		double calculateL(MyMesh::VertexHandle* omega, std::vector<int>* neighbours, int idx);

		void findSingularities();
		void findSeparatrices();
		void findSeparatrices2();
		void followStreamLine(Vector v, Vector prevDir, std::vector<Vector>* streamline,int separatriceIdx, double step = 0.5, double iterations = 2000);
		void buildSeparatrices(std::vector<Vector>* separatice,Vector dir, MyMesh::VertexHandle v1, MyMesh::VertexHandle v2, MyMesh::FaceHandle f);
		void drawSeparatrices();
		void findXiNext(Vector Xi, Vector di, MyMesh::VertexHandle v1, MyMesh::VertexHandle v2, MyMesh::VertexHandle v3, Vector* Xin, Vector* din, MyMesh::VertexHandle* newV1, MyMesh::VertexHandle* newV2);
		std::vector<Vector> calculateCrossFromU(Vector u);
		Vector findClosestCrossVector(Vector di, Vector Xi, MyMesh::VertexHandle v1, MyMesh::VertexHandle v2);
		void buildStreamLine(Vector v);
		MyMesh::FaceHandle getFace(Vector v, double* w, MyMesh::VertexHandle* v_array);
		void findPartitionCorners();
		bool doIntersect(Vector x1, Vector x2, Vector x3, Vector x4, Vector* P);
	bool is_collapse_ok2(MyMesh::HalfedgeHandle v0v1);
	//////////////////////
	// Member variables //
	//////////////////////

	enum class ModelType { NONE, MESH, BEZIER_SURFACE } model_type;

	// Mesh
	MyMesh mesh;

	MyMesh keepInContsraint;
	MyMesh keepOutConstraint;

	struct Singularity {
		Vector pos;
		MyMesh::FaceHandle f;
		Singularity(Vector pos_, MyMesh::FaceHandle f_) : pos(pos_), f(f_) {}
	};
	std::vector<Singularity> singularities;
	std::vector<std::vector<Vector>> separatrices;
	std::vector<std::vector<Vector>> separatrices2;
	OpenMesh::PolyMesh_ArrayKernelT<MyTraits> quadPartition;
	// Bezier
	size_t degree[2];
	std::vector<Vec> control_points;

	//Constraint solids
	std::vector<ConstraintSolid*> constraint_solids;

	// Visualization
	double mean_min, mean_max, cutoff_ratio;
	bool show_control_points, show_solid, show_wireframe, show_constraints, show_partitioning, show_PDEu, show_PDEcross, show_singularities, show_separatrices;
	double boundaryRemeshL;
	enum class Visualization { PLAIN, MEAN, SLICING, ISOPHOTES } visualization;
	GLuint isophote_texture, environment_texture, current_isophote_texture, slicing_texture;
	Vector slicing_dir;
	double slicing_scaling;
	int selected_vertex;
	struct ModificationAxes {
		bool shown;
		float size;
		int selected_axis;
		Vec position, grabbed_pos, original_pos;
	} axes;
	std::string last_filename;
};

#include "MyViewer.hpp"
