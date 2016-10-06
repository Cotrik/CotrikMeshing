/*
 * GeoUtil.h
 *
 *  Created on: Feb 6, 2015
 *      Author: cotrik
 */

#ifndef __GEOUTIL_H__
#define __GEOUTIL_H__

#include "GeometricStruct.h"

class GeoUtil
{
public:
	GeoUtil();
	virtual ~GeoUtil();

	static void GetMaxCoordinateValueOfVertex(const Vertex& vertex, Vertex& vertexMax);
	static void GetMinCoordinateValueOfVertex(const Vertex& vertex, Vertex& vertexMin);
	static void GetMaxMinCoordinateValueOfVertex(const Vertex& vertex, Vertex& vertexMax, Vertex& vertexMin);
	static void GetMaxMinCoordinateValueOfTriangle(const std::vector<Vertex*>& m_vecVertex, const Triangle& triangle, Vertex& vertexMax, Vertex& vertexMin);
	static bool IsPointInTriangle(const glm::vec3& P, const glm::vec3& A, const glm::vec3& B, const glm::vec3& C, const int i = 0, const int j = 0);
	static bool IsPointCanProjectToTriangle(const glm::vec3& orig, const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2, glm::vec3& p);
	static bool IsPointInTetrahedron(const Vertex& p, const Vertex& p0, const Vertex& p1, const Vertex& p2, const Vertex& p3);
	static bool IsPointInTetrahedron(const Vertex& p, const Vertex& p0, const Vertex& p1, const Vertex& p2, const Vertex& p3, glm::vec3& uvw);
	static bool IsPointInTetrahedron_Robust(const Vertex& p, const Vertex& p0, const Vertex& p1, const Vertex& p2, const Vertex& p3, glm::vec4& uvwx);
	static bool IsPointInTetrahedron_Robust_test(const Vertex& p, const Vertex& p0, const Vertex& p1, const Vertex& p2, const Vertex& p3, glm::vec4& uvwx);
	static bool IsVertexInsideTetrahedron(const Vertex& vertex, const Cell& tet, const std::vector<Vertex>& vecVertex, glm::vec3& lambda);
	static bool IsVertexInsideTetrahedron_Robust(const Vertex& p, const Vertex& p0, const Vertex& p1, const Vertex& p2, const Vertex& p3, glm::vec3& lambda);
	static bool IsTetrahedronDegenerated(const Cell& tet, const std::vector<Vertex>& vecVertex, const double eps = 1e-9);
	static FACE_TYPE GetFaceType(const glm::vec3& normal);
};

#endif /* __GEOUTIL_H__ */
