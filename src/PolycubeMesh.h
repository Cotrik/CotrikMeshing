/*
 * PolycubeMesh.h
 *
 *  Created on: Jan 3, 2015
 *      Author: cotrik
 */

#ifndef __POLYCUBEMESH_H__
#define __POLYCUBEMESH_H__

#include "Mesh.h"
#include "MeshFileReader.h"
#include "Patch.h"

typedef std::multimap<Plane, Face*> pf;
typedef std::multimap<double, Face*> df;
typedef std::multimap<double, Line*> dl;
typedef std::multimap<Vector, Line*> vl;

extern unsigned int THREAD_NUM;
extern double g_upper_size;
extern double g_lower_size;
//void GetMaxMinCoordinateValueOfVertex(const Vertex& vertex, Vertex& vertexMax, Vertex& vertexMin);
class PolycubeMesh : public Mesh
{
public:
	PolycubeMesh(const std::vector<Vertex>& v, const std::vector<Cell>& c, const ElementType cellType = HEXAHEDRA);
	PolycubeMesh(const Mesh& mesh);
	PolycubeMesh();
	virtual ~PolycubeMesh();

//	void NormalizeCoordinateValue();

	void GetPlane_Face_MultiMap();
	void GetXYZPlane_Face_MultiMap();
	void GetPlane_Face_Multimap(df& m_plane_Face_Multimap, const DIRECTION dir = X_AXIS);
	void GetSinglePlaneFaceMultimap(const pf& plane_Face_Multimap, df& singlePlane_Face_Multimap, const float delta = 0.005);
	void GetBoundaryLinesOfMeshPlane(const std::vector<Face*>& vecSinglePlaneFace, std::vector<Line>& vecBoundaryLine);
	void GetPlaneBoundaryLine(const df& singlePlane_Face_Multimap, std::vector<Line>& vecBoundaryLine);

	void GetPlaneMultimap(const std::vector<Line>& vecBoundaryLine, vl& planeLineMultimap);

	void Getxy_xzPlaneLineMultimap(const vl& xPlaneLineMultimap, dl& xyPlaneLineMultimap, dl& xzPlaneLineMultimap);
	void Getyx_yzPlaneLineMultimap(const vl& yPlaneLineMultimap, dl& yxPlaneLineMultimap, dl& yzPlaneLineMultimap);
	void Getzx_zyPlaneLineMultimap(const vl& zPlaneLineMultimap, dl& zxPlaneLineMultimap, dl& zyPlaneLineMultimap);

	void GetXYPlaneCorner(const dl& xyPlaneLineMultimap, std::vector<unsigned long>& vec_xyPlaneCorner);
	void GetCorners();
	void GenerateSmallervolumetricCubeMesh(const float& delta, const char* outputFilename, const MESH_TYPE meshType);
	void AdjustPolycubeSurfaceToHoldAllSmallCubesVolume(const float& delta, const Vertex& minVertex, const char* outputFilename = NULL);
	void AddCompensation();
	void AddCompensation_round();

	void AdjustPolycubeSurfaceToHoldAllSmallCubesVolume_quad(const float& delta, const Vertex& minVertex, const char* outputFilename = NULL);
	void AdjustPolycubeSurfaceToHoldAllSmallCubesVolume1(const float& delta, const Vertex& minVertex, const char* outputFilename = NULL);
	void AdjustPolycubeSurfaceToHoldAllSmallCubesVolume1_Reverse(const float& delta, const Vertex& minVertex, const char* outputFilename = NULL);
	void AdjustPolycubeSurfaceToHoldAllSmallCubesVolume1_ROUND(const float& delta, const Vertex& minVertex, const char* outputFilename = NULL);
	void AdjustPolycubeSurfaceToHoldAllSmallCubesVolume2_ROUND(const float& delta, const Vertex& minVertex, const char* outputFilename = NULL);
	void AdjustPolycubeSurfaceToHoldAllSmallCubesVolume_round(const float& delta, const Vertex& minVertex, const char* outputFilename = NULL);
	void AddCompensation_quad();
	void AddCompensation_quad_round();

	void RoundSurface(const double eps = 1e-3); // eps = 1e-3 so that there are 10^3 cubes along x, y, z axis;
	void RoundSurfaceAxis(const double eps = 1e-3);
	void Magnify(const double eps = 1e-3);

	void GenerateHexVertices(std::vector<Vertex>& hexV); //Hex vertices are generated first for the PolyCube corners, then edges, then facets and finally volume
	void GenerateHexVerticesForCorners(std::vector<Vertex>& hexV, std::vector<int>& cornerIds);
	void GenerateHexVerticesForEdges(std::vector<Vertex>& hexV, std::vector<Edge>& edges);
	void GenerateHexVerticesForFacets(std::vector<Vertex>& hexV, std::vector<int>& facetIds);

	void HexMeshing(const char* outputPolycubeHexFilename, const float delta, const DIRECTION dir);


	void GetVertexInfo();
	void Smooth(const unsigned int time = 1);
	void SmoothBoundaryLine(const unsigned int time = 1);
	void SmoothBoundaryLine(const int patch, const unsigned int time = 1);
	virtual void SmoothSurface(const unsigned int time = 1);
	void SmoothVolume(const unsigned int time = 1);
	void SmoothVolumeP(const unsigned int time = 1);
	void SmoothBoundaryCell(const unsigned int time = 1);
	void AjustSurfaceToIntegerPlane(const float delta = 0.005);

	void ParticleSmooth(const unsigned int time = 1);
	void GetSurfaceNeignbor(const unsigned long vertexIndex, std::vector<unsigned long>& vecNeighbor);
	void GetBoundaryLineNeignbor(const unsigned long vertexIndex, std::vector<unsigned long>& vecNeighbor);
	void GetE_ij_and_F_ij(const unsigned long i, const unsigned long j,	double& E_ij, glm::vec3& F_ij);

	void CurvatureWeightingSmooth(const unsigned int time = 1);
	void CurvatureWeightingSmoothBoundaryLine(const unsigned int time = 1);
	void CurvatureWeightingSmoothSurface(const unsigned int time = 1);
	void CurvatureWeightingSmoothVolume(const unsigned int time = 1);

	void LaplacianWeightingSmoothBoundaryLine(const unsigned int time = 1);

	void SmoothCurvature(const unsigned int time = 1);

	///*virtual */void GetMaxMinCoordinates();
	bool IsVerticesOnTheSamePlane(const std::vector<unsigned long>& surfaceVertices, DIRECTION& d);
	////////////////////////////////////////////
	// deal with patches
	void LabelSurfaceFace();
	void LabelSinglePlaneFace(const df& singlePlane_Face_Multimap, int& label);
	void GetAreaPatches();
	void GetFacePatches();
	void GetEdgePatches();
	void GetVertexPatches();
	void CheckPatchesConnection();
	void SortPatches(const PolycubeMesh& polycubeTetMesh);
	void SortPatches();
	int DividePatch(const std::vector<unsigned long>& vertexPatch, std::vector<std::vector<unsigned long> >& vertex_patches, const DIRECTION dir = X_AXIS);

	void GetBoundaryEdge(const std::vector<unsigned long>& faceIndex, std::vector<Edge>& edgePatch);
	void GetPatchesLabel();
	void Cleanup();
	void RelabelSurfaceFacePerFace();
	void RelabelIsolatedPatch();
	void RelabelWedgePatch();
	void ComputeWedgeTwoEdgePathlengths(const int patchLabelIndex, double* pathLength, int* path_label);
	////////////////////////////////////////////
	pf  m_plane_Face_Multimap;
	pf m_xPlane_Face_Multimap;
	pf m_yPlane_Face_Multimap;
	pf m_zPlane_Face_Multimap;

	df m_xSinglePlane_Face_Multimap;
	df m_ySinglePlane_Face_Multimap;
	df m_zSinglePlane_Face_Multimap;

	df m_x_Face_Multimap;
	df m_y_Face_Multimap;
	df m_z_Face_Multimap;

	std::vector<Line> m_vecXBoundaryLine;
	std::vector<Line> m_vecYBoundaryLine;
	std::vector<Line> m_vecZBoundaryLine;

	vl m_xPlaneLineMultimap;
	vl m_yPlaneLineMultimap;
	vl m_zPlaneLineMultimap;

	std::multimap<double, Line*> m_xyPlaneLineMultimap;
	std::multimap<double, Line*> m_xzPlaneLineMultimap;
	std::multimap<double, Line*> m_yzPlaneLineMultimap;

	std::vector<unsigned long> m_vec_xyPlaneCorner;
	std::vector<unsigned long> m_vec_xzPlaneCorner;
	std::vector<unsigned long> m_vec_yzPlaneCorner;

	std::vector<Vertex*> m_vec_xyPlaneCornerVertex;
	std::vector<Vertex*> m_vec_xzPlaneCornerVertex;
	std::vector<Vertex*> m_vec_yzPlaneCornerVertex;

	std::vector<LineVertices> m_vec_xyPlaneCornerLine;
	std::vector<LineVertices> m_vec_xzPlaneCornerLine;
	std::vector<LineVertices> m_vec_yzPlaneCornerLine;

	std::vector<Vertex*> vecCubeCorner;

	std::vector<glm::vec3> m_vecVertexCompensation;

	std::map<unsigned long, glm::vec3> m_map_point_meanCoordinates;

	std::vector<double> m_patchPos;
	std::vector<double> m_xpatchInterval;
	std::vector<double> m_ypatchInterval;
	std::vector<double> m_zpatchInterval;
	double m_min_distance;
	double m_min_distance_x;
	double m_min_distance_y;
	double m_min_distance_z;

	void GetPatchMinDistance();
	void GetPatchMinDistance_round();
	unsigned int DividePolycube(const Vertex& originPoint, std::vector<Vertex>& vecVertex,	std::vector<Cell>& vecHexahedron);
	unsigned int DividePolycube_aniso(const Vertex& originPoint, std::vector<Vertex>& vecVertex,  std::vector<Cell>& vecHexahedron,
	        const glm::vec3& cellScale, const float delta = 0.002f);
	unsigned int DividePolycube_round(const Vertex& originPoint, std::vector<Vertex>& vecVertex,  std::vector<Cell>& vecHexahedron);
	//////////////////////////////////////////
	void GetFacePatches_N(); // using neighboring method
	void GetFacePatches_N(const unsigned long faceId, int& patchIndex); // using neighboring method
	void GetCellPatches_N(); // using neighboring method
	void GetCellPatches_N(const unsigned long cellId, int& patchIndex); // using neighboring method

	void GetPatches();
	void ModifyPatchesPosition();
	void ModifyPatchesPosition_ROUND();
	void ModifyPatchesPosition_ROUND(float delta);
	void ModifyPatchesPosition_ROUND2(float delta);

	void AssignPatches();
};

#endif /* __POLYCUBEMESH_H__ */
