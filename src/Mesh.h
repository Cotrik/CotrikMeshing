/*
 * Mesh.h
 *
 *  Created on: Dec 26, 2014
 *      Author: cotrik
 */

#ifndef __MESH_H__
#define __MESH_H__

#include <vector>
#include <map>
#include <algorithm>

#include "GeometricStruct.h"
#include "GeoUtil.h"
#include <vtkCellType.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkSmartPointer.h>
#include <vtkDataSetMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkCellIterator.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkPolyData.h>

#include "Patch.h"

#include "vtkSmartPointer.h"

#include "vtkPolyData.h"
#include "vtkPolyDataMapper.h"
#include "vtkActor.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkTimerLog.h"
#include "vtkCamera.h"
#include "vtkProperty.h"
#include "vtkInteractorStyleTrackballCamera.h"
#include "vtkPolyDataNormals.h"
#include "vtkRendererCollection.h"
#include "vtkPolyDataCollection.h"
#include "vtkObjectFactory.h"
#include "vtkIdList.h"
#include "vtkPolyDataReader.h"
#include "vtkXMLPolyDataReader.h"
#include "vtkXMLPolyDataWriter.h"
#include "vtkNew.h"
#include "vtkPointData.h"
#include "vtkContourWidget.h"
#include "vtkOrientedGlyphContourRepresentation.h"
#include "vtkPolygonalSurfacePointPlacer.h"

bool IsVertexInsideTetrahedron(const Vertex& vertex, const Cell& tet, const std::vector<Vertex>& vecVertex);

struct FaceInfo
{
	FACE_TYPE faceType;
};

struct FacePair
{
	FacePair(const Face& _sortedFace, const Face& _originFace)
	: sortedFace(_sortedFace), originFace(_originFace)
	{}
	Face sortedFace;
	Face originFace;
	bool operator < (const FacePair& right) const
	{
		return sortedFace < right.sortedFace;
	}
	bool operator == (const FacePair& right) const
	{
		return sortedFace == right.sortedFace;
	}
};

class SFace : public Face
{
public:
	SFace(){};
	SFace(const Face& face)
	{
		resize(face.size());
		std::copy(face.begin(), face.end(), begin());
	}
	bool operator < (const Face& right) const
	{
		if (size() == 3)
		{
			if (at(0) < right.at(0))
			{
				return true;
			}
			else if (at(0) == right.at(0))
			{
				if (at(1) < right.at(1))
				{
					return true;
				}
				else if (at(1) == right.at(1))
				{
					if (at(2) < right.at(2))
					{
						return true;
					}
				}
			}
		}

		return false;
	}
	bool operator == (const Face& right) const
	{
		if (size() == 4)
			return ((at(0) == right.at(0) && at(1) == right.at(3) && at(2) == right.at(2) && at(3) == right.at(1))
			     || (at(0) == right.at(0) && at(1) == right.at(1) && at(2) == right.at(2) && at(3) == right.at(3))

			     /*|| (at(0) == right.at(0) && at(1) == right.at(3) && at(2) == right.at(2)) && at(3) == right.at(1)
			     || (at(0) == right.at(3) && at(1) == right.at(2) && at(2) == right.at(1)) && at(3) == right.at(0)
			     || (at(0) == right.at(2) && at(1) == right.at(1) && at(2) == right.at(0)) && at(3) == right.at(3)
                 || (at(0) == right.at(1) && at(1) == right.at(0) && at(2) == right.at(3)) && at(3) == right.at(2)

			     || (at(0) == right.at(0) && at(1) == right.at(3) && at(2) == right.at(2)) && at(3) == right.at(1)
			     || (at(0) == right.at(3) && at(1) == right.at(2) && at(2) == right.at(1)) && at(3) == right.at(0)
			     || (at(0) == right.at(2) && at(1) == right.at(1) && at(2) == right.at(0)) && at(3) == right.at(3)
                 || (at(0) == right.at(1) && at(1) == right.at(0) && at(2) == right.at(3)) && at(3) == right.at(2)

			     || (at(1) == right.at(0) && at(2) == right.at(3) && at(3) == right.at(2)) && at(0) == right.at(1)
			     || (at(1) == right.at(3) && at(2) == right.at(2) && at(3) == right.at(1)) && at(0) == right.at(0)
			     || (at(1) == right.at(2) && at(2) == right.at(1) && at(3) == right.at(0)) && at(0) == right.at(3)
                 || (at(1) == right.at(1) && at(2) == right.at(0) && at(3) == right.at(3)) && at(0) == right.at(2)

			     || (at(2) == right.at(0) && at(3) == right.at(3) && at(0) == right.at(2)) && at(1) == right.at(1)
			     || (at(2) == right.at(3) && at(3) == right.at(2) && at(0) == right.at(1)) && at(1) == right.at(0)
			     || (at(2) == right.at(2) && at(3) == right.at(1) && at(0) == right.at(0)) && at(1) == right.at(3)
                 || (at(2) == right.at(1) && at(3) == right.at(0) && at(0) == right.at(3)) && at(1) == right.at(2)

			     || (at(3) == right.at(0) && at(0) == right.at(3) && at(1) == right.at(2)) && at(2) == right.at(1)
			     || (at(3) == right.at(3) && at(0) == right.at(2) && at(1) == right.at(1)) && at(2) == right.at(0)
			     || (at(3) == right.at(2) && at(0) == right.at(1) && at(1) == right.at(0)) && at(2) == right.at(3)
                 || (at(3) == right.at(1) && at(0) == right.at(0) && at(1) == right.at(3)) && at(2) == right.at(2)*/
				 );
		else if (size() == 3)
			return ((at(0) == right.at(0) && at(1) == right.at(2) && at(2) == right.at(1))
			     || (at(0) == right.at(2) && at(1) == right.at(1) && at(2) == right.at(0))
			     || (at(0) == right.at(1) && at(1) == right.at(0) && at(2) == right.at(2))

			     || (at(1) == right.at(0) && at(2) == right.at(2) && at(0) == right.at(1))
				 || (at(1) == right.at(2) && at(2) == right.at(1) && at(0) == right.at(0))
				 || (at(1) == right.at(1) && at(2) == right.at(0) && at(0) == right.at(2))

				 || (at(2) == right.at(0) && at(0) == right.at(2) && at(1) == right.at(1))
				 || (at(2) == right.at(2) && at(0) == right.at(1) && at(1) == right.at(0))
				 || (at(2) == right.at(1) && at(0) == right.at(0) && at(1) == right.at(2))

				 || (at(0) == right.at(0) && at(1) == right.at(1) && at(2) == right.at(2)));
		else
			return *this == right;
	}
};

class Mesh
{
public:
	Mesh(const std::vector<Vertex>& v, const std::vector<Cell>& c, const ElementType cellType = HEXAHEDRA);
	Mesh(const Mesh& mesh);
	Mesh();
	virtual ~Mesh();
	std::vector<Vertex> V;
	std::vector<Edge> E;
	std::vector<Face> F;
	std::vector<Cell> C;

	std::multimap<unsigned long, unsigned long> VI_VI;   // Vid_Vid pair
	std::multimap<unsigned long, unsigned long> VI_EI;   // Vid_Eid pair
	std::multimap<unsigned long, unsigned long> VI_FI;   // Vid_Fid pair
	std::multimap<unsigned long, unsigned long> VI_CI;   // Vid_Fid pairAAAA
	std::multimap<unsigned long, unsigned long> EI_FI;   // Vid_Fid pair
	std::multimap<unsigned long, unsigned long> EI_CI;   // Eid_Cid pair
	std::multimap<unsigned long, unsigned long> FI_CI;   // Fid_Cid pair

    std::vector<glm::vec3> N_F;  // normal of each face in F
    std::vector<glm::vec3> N_V;  // normal of each vertex in V
	std::multimap<unsigned long, Edge> V_E;
	std::multimap<unsigned long, SFace> V_F;
	//std::multimap<unsigned long, Cell> V_C;
//	std::multimap<unsigned long, unsigned long> V_CI;

	std::multimap<unsigned long, unsigned long> V_V;
	std::map<unsigned long, std::vector<unsigned long> > V_Vs;
	std::map<unsigned long, std::vector<unsigned long> > V_Fs;

	std::multimap<unsigned long, unsigned long> C_C; // Cell and its neighbor Cells
	std::multimap<unsigned long, unsigned long> F_F; // Face and its neighbor Faces

	//std::vector<VertexInfo> vertexInfo;
	inline const size_t Size() const
	{
		return V.size();
	}

	std::vector<Face> surface;
	std::vector<Face> sortedSurface;

	int m_xLabel;
	int m_yLabel;
	int m_zLabel;

	std::vector<FACE_TYPE> faceType;            // faceType of face on the surface
	std::vector<int> faceLabel;                 // face Label of face on the surface
	std::vector<double> faceArea;
	std::vector<std::vector<int> > vertexLabel;  // size is number of Vertices
	std::vector<std::vector<int> > patchLabel;
	int m_patchNumber;
	std::vector<std::vector<unsigned long> > facePatches;  // the begin patch's lable is 1, next is 2, and so on. unsigned long is the face index in surface
	std::vector<std::vector<Edge> > edgePatches;  // the begin patch's lable is 1, next is 2, and so on
	std::vector<std::vector<unsigned long> > vertexPatches; // the begin patch's lable is 1, next is 2, and so on
	//std::vector<std::vector<int> > vertexLabel; // the begin patch's lable is 1, next is 2, and so on
	std::vector<unsigned long> surfaceVertexIndices;    // get surface Vertex Indices in V
	MESH_TYPE m_meshType;
	ElementType m_cellType;

	Vertex m_maxVertex;
	Vertex m_minVertex;

	std::vector<Vertex> m_vecPolycubeHexVertex;
	std::vector<Cell> m_vecPolycubeHexCell;

	std::vector<Vertex> m_vecHexVertex;
	std::vector<Cell> m_vecHexCell;

	void ExtractSurface();
	glm::vec3 GetFaceNormal(const Face& f) const;
	void GetNormalOfSurfaceFaces();             // must ExtractSurface(); first
	void GetNormalOfSurfaceVertices();          // must ExtractSurface(); first
	void NormalizeCoordinateValue();
	void GetCenterPointOfEdge(const Edge& e, glm::vec3& centerPoint);
	void GetCenterPointOfCell(const Cell& c, glm::vec3& centerPoint);
	void GetCurvatureCenterPointOfEdge(const Edge& e, glm::vec3& centerPoint);
	void GetCurvatureCenterPointOfCell(const Cell& c, glm::vec3& centerPoint);
	void GetCurvatureCenterPointOfCell_P(const Cell& c, glm::vec3& centerPoint);
	bool IsEdgeOnBondaryLine(const Edge& e);
	bool IsEdgeOnPatch(const Edge& e, const int patchLabel);
	bool IsFaceOnSurface(const Face& f);
	bool IsFaceOnPatch(const Face& f, const int patchLabel);
	const DIRECTION GetFaceAxis(const Face& f);
	void GetCurvature(const char* surfaceCurvatureFilename = NULL);
	double GetCenterCurvatureOfCell(const Cell& c);

	void CurvatureWeightingSmoothSurface(const unsigned int time = 1);
	virtual void GetMaxMinCoordinates();
	virtual void GetMaxMinCoordinates_Patch();

	void WriteHexahedralmesh(const char* outputFilename);
	void GenerateSmallCubes(const float& delta, const char* outputFilename);
	void GenerateSmallCubes_plus(const float& delta, const char* outputFilename, const float plus = 0.001f);
	void GenerateSmallCubes_Magifier(const float& delta);
	void removeInvalidHexahedronsInBox(const float delta = 0.008f, const DIRECTION dir = Z_AXIS);
	void removeInvalidHexahedronsInBox_ref(std::vector<Vertex>& Vs, const float delta = 0.008f, const DIRECTION dir = Z_AXIS);
	void removeInvalidHexahedronsInBox_new(const float delta = 0.008f, const DIRECTION dir = Z_AXIS);
	void removeInvalidHexahedronsInBox_fast(const float delta = 0.008f, const DIRECTION dir = Z_AXIS);
	void removeInvalidHexahedronsInBox_Magnifier(const float delta, const DIRECTION dir);
	void GenerateHexahedralMesh(const Mesh& orgTet,
			const char* outputPolycubeHexFilename = NULL, const char* outputHexFilename = NULL,
			const float delta = 0.008f, const DIRECTION = Z_AXIS);
	void GenerateHexahedralMesh_Magnifier(const Mesh& orgTet,
				const char* outputPolycubeHexFilename = NULL, const char* outputHexFilename = NULL,
				const float delta = 0.008f, const DIRECTION = Z_AXIS);
	void GenerateHexahedralMesh_plus(const Mesh& orgTet,
				const char* outputPolycubeHexFilename = NULL, const char* outputHexFilename = NULL,
				const float delta = 0.008f, const DIRECTION = Z_AXIS);
	void GenerateHexahedralMesh_align(const Mesh& orgTet,
				const char* outputPolycubeHexFilename = NULL, const char* outputHexFilename = NULL,
				const float delta = 0.008f, const DIRECTION = Z_AXIS);
	void GenerateHexahedralMesh_fast(const Mesh& orgTet,
				const char* outputPolycubeHexFilename = NULL, const char* outputHexFilename = NULL,
				const float delta = 0.008f, const DIRECTION = Z_AXIS);
	unsigned long GetCubeIndex(const Vertex& v, const float delta = 0.08f) const;
	unsigned int DivideHexahedronIntoSmallerOnes(const Vertex& originPoint, const Vector& axis, const float& xInterval,
			const float& yInterval, const float& zInterval, std::vector<Vertex>& vecVertex,	std::vector<Cell>& vecHexahedron);
	void CreateHexMesh(const Mesh& origTetMesh, std::vector<Vertex>& vecHexVertex, std::vector<Cell>& vecHexCell);
	void AddVertex(const std::vector<Vertex>& vecOrgTetVertex, Vertex& hexVertex, const Cell& tet,
			const std::vector<Vertex>& vecPolycubeVertex, std::vector<Vertex>& vecHexVertex);

	virtual void SmoothVolume(const unsigned int time = 1);
	virtual void SmoothSurface(const unsigned int time = 1);

	void GetParameterAndCubeIndicesOfVertexIndices();
	void MapbackFromDeformedCubes();
	void GetDensityField();
	void GetSurfaceDensityField();
	void OutputScaledJacobianDataFile(const std::vector<Cell>& C, const char* filename = "scaled_jacobian.csv");
	const float GetScaledJacobian(const Cell& c);
	bool JudgeDirection(const Cell& c);

	void GetSurfaceVertexIndices();
	void VTKExtractSurface();

	bool IsPointInsideMe(const Vertex& p);
	static bool IsPointInsideMesh(const Vertex& p, const Mesh& mesh);
	bool IsAllVerticesInsideMesh(const Mesh& mesh);

	void GetNeighboringInfo();
	void FixMesh();

	void InitGeodesicDistance(const char* polydataFilename);
	double ComputeGeodesicDistance(vtkPolyData *polyData, vtkIdType beginVertId = -1, vtkIdType endVertId = -1);

	//////////////////////////////////////////
	// cone optimization
	std::vector<DirectedEdge> DE;
//	std::multimap<DirectedEdge, DirectedEdge> m_directedEdge_surroundingDirectedEdge;
//	std::multimap<Edge, unsigned long/*SFace*/> E_F;
	std::multimap<DirectedEdge, DirectedEdge> E_E;
	std::vector<SFace> SF;
	std::vector<Edge> uniqueE;
	std::vector<SFace> uniqueF;
	unsigned int m_genus;
	void GetConeConfiguration();
	void GetLocalMinimalConeVertex(const std::multimap<DirectedEdge, DirectedEdge>::const_iterator it, glm::vec3& opt_v);
	double GetSingleEcone(const DirectedEdge& e_ij, const DirectedEdge& u_k_kplus);
	double GetSingleEcone(const DirectedEdge& e_ij, const DirectedEdge& u_k_kplus, const glm::vec3& nk_c);
	void OptimizeConeShape(int iterMaxNum = 1);
	void WriteConeShape(const char* filename, const std::multimap<DirectedEdge, DirectedEdge>::const_iterator it);

	bool VerifyStructure(const int genus = 0);
	bool VerifyElements() const;
	int GetGenus();
	//
	void SetFixedVertices(std::vector<unsigned long>& fixed_points);

public:
	void ExtractSurfaceStep(const Face& face);
	void GetF();
	void GetE();
    void GetVI_VI();
    void GetVI_EI();
    void GetVI_FI();
    void GetVI_CI();
    void GetEI_FI();
    void GetEI_CI();
    void GetFI_CI();
	unsigned long GetEI(const Edge& e);

//public:
	void Init();
private:
	void InitVertexInfo();
	void InitFace();
public:
	std::vector<Cell> tets;
	std::vector<glm::vec3> paras;
	std::vector<bool> paras_flag;
	std::vector<unsigned long> vec_cubeIndex;
public:
	std::vector<float> vec_densityFiled;
	std::vector<float> vec_averageLength_V;
	std::vector<unsigned long> invertedCellIndex;
	std::vector<unsigned long> lowQualityCellIndex;

	vtkUnstructuredGrid* pUnstructuredGrid;

	std::vector<unsigned long> outsidePoints;
	std::vector<unsigned long> highlightPatchLabel;
	std::vector<std::vector<std::vector<unsigned long> > > m_vertex_patches;
	std::vector<std::vector<unsigned long> > m_total_vertex_patches;
	int m_vertex_patches_num;
	int m_vertex_patches_num_x;
	int m_vertex_patches_num_y;
	int m_vertex_patches_num_z;

	std::vector<Patch> m_patches;

	std::vector<unsigned long> overlappingVertices;
	////////////////////////////////////////////////
	void GetCellAndNeighborCells();
	void GetFaceAndNeighborFaces();
	void GetNeighborCells(const unsigned long cellId, std::vector<unsigned long>& neighborCellIds);
	void GetNeighborFaces(const unsigned long faceId, std::vector<unsigned long>& neighborFaceIds);

	const glm::vec3 GetNormal(unsigned long faceIndex);
	void GetFaceType();

	void SmoothPatches();
	////////////////////////////////////////////////
	std::vector<std::vector<double> > scalarFields;
	std::vector<std::vector<double> > vectorFields;
	bool IsPointInPoints(unsigned long pointId, const std::vector<unsigned long>& Ids);
	bool IsCutPoint(unsigned long pointId, const std::vector<unsigned long>& Ids);
	int GroupCutPoints(const std::vector<unsigned long>& Ids, std::vector<std::vector<unsigned long> >& groupIds);

	////////////////////////////////////////////////
	// The output are edges indices for volumetric mesh or vertices indices for surface mesh
	std::vector<unsigned long> singularities;
	void GetSingularities();
};

extern float __plus;
#endif /* __MESH_H__ */
