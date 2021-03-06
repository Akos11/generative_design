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
		};
	};
	using MyMesh = OpenMesh::TriMesh_ArrayKernelT<MyTraits>;
	using Vector = OpenMesh::VectorT<double, 3>;
	class generativeDesign;

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
		void quadrangulate();
			void makeEvenTriangles();
			void makeQuadDominant();
				void calculateSquarness();
				bool canDelete(MyMesh::EdgeHandle eh);
				void deleteEdges();
				MyMesh::FaceHandle delete_edge(MyMesh::EdgeHandle _eh);
				MyMesh::FaceHandle add_edge(MyMesh::FaceHandle fh, MyMesh::VertexHandle vh0, MyMesh::VertexHandle vh1);
			void makePureQuad(int idx);
				int faceSides(MyMesh::FaceHandle fh);
				void resetFaceFlags();
				bool isVertexOnFace(MyMesh::VertexHandle vh, MyMesh::FaceHandle fh);
				bool hasCommonVertex(MyMesh::FaceHandle fh1, MyMesh::FaceHandle fh2);
				bool hasCommonEdge(MyMesh::VertexHandle vh1, MyMesh::VertexHandle vh2);
				bool hasCommonFace(MyMesh::VertexHandle vh1, MyMesh::VertexHandle vh2);
				MyMesh::EdgeHandle getCommonEdge(MyMesh::VertexHandle vh1, MyMesh::VertexHandle vh2);
				MyMesh::EdgeHandle getCommonEdge(MyMesh::FaceHandle fh1, MyMesh::FaceHandle fh2);
				double getAngleSum(MyMesh::FaceHandle fh);
				double getAngleSum(std::vector<MyMesh::VertexHandle> vertices);
			void smoothQuadMesh(int iterations = 5);
		void quadRegularization(int iterations = 5);
			void eliminate2and4ValenceBounradies();
			void printUnregularVertices();
			void quadRegularizationRemove2ValenceVertices();
			void quadRegularizationSwapping();
			void quadRegularizationCollapsing();
			void collapseVertices(MyMesh::FaceHandle fh, MyMesh::VertexHandle vh1, MyMesh::VertexHandle vh2);
			void quadRegularizationSplitting();
			void quadRegularizationCompositons();
				void quadRegularizationCompositonSwapSplit();
				void quadRegularizationCompositonSwapCollapse();
				void quadRegularizationCompositonSwapCollapseHelper(MyMesh::VertexHandle v, bool otherNeighbour = false);
				void quadRegularizationCompositonSplitSplit();
				void quadRegularizationCompositonCollapseCollapse();
				void quadRegularizationCompositonSwapSwap();
				void quadRegularizationCompositonSwapSwap2();
			void quadRegularizationTransfer();
			void quadRegularizationTransfer2();
				MyMesh::VertexHandle findWaytoNearestUnregular(MyMesh::VertexHandle vh, MyMesh::VertexHandle notThis);
		void partition();
			int counter;

	//generative ver 2.0
			void createQuadPartitioning(double maxLength);
			void createQuadPartitioningParametrization(double _maxLength, double angleWeight);
			int findPointsInBoundary(double maxLength);
			struct Comparer
			{
				MyMesh mesh;
				Comparer(MyMesh _mesh)
				{
					mesh = _mesh;
				}
				inline bool operator() (const MyMesh::VertexHandle& v1, const MyMesh::VertexHandle& v2);
			
			};
			void createQuadsFromTagged(int taggedNum);
			bool isLeft(MyMesh::Point a, MyMesh::Point b, MyMesh::Point c);
			double angleInBoundary(MyMesh::VertexHandle v);
			void smoothQuadMesh2(int iterations = 5);
			void reMeshVertexPositions2(bool remeshUnlockedBoundaries = false);
			void drawPartitioning();
	//void collapse2(MyMesh::HalfedgeHandle _hh);
	//void collapse_edge2(MyMesh::HalfedgeHandle _hh);
	//void collapse_loop2(MyMesh::HalfedgeHandle _hh);


	bool is_collapse_ok2(MyMesh::HalfedgeHandle v0v1);
	//////////////////////
	// Member variables //
	//////////////////////

	enum class ModelType { NONE, MESH, BEZIER_SURFACE } model_type;

	// Mesh
	MyMesh mesh;

	MyMesh keepInContsraint;
	MyMesh keepOutConstraint;

	OpenMesh::PolyMesh_ArrayKernelT<MyTraits> quadPartition;
	// Bezier
	size_t degree[2];
	std::vector<Vec> control_points;

	//Constraint solids
	std::vector<ConstraintSolid*> constraint_solids;

	// Visualization
	double mean_min, mean_max, cutoff_ratio;
	bool show_control_points, show_solid, show_wireframe, show_constraints, show_partitioning;
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
