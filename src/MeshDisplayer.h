/*
 * MeshDisplayer.h
 *
 *  Created on: Feb 13, 2015
 *      Author: cotrik
 */

#ifndef __MESHDISPLAYER_H__
#define __MESHDISPLAYER_H__

#include "Mesh.h"
#include "PolycubeMesh.h"

class MeshDisplayer
{
public:
	MeshDisplayer(const Mesh& mesh);
	virtual ~MeshDisplayer();

	void Draw() const;
private:
	/*static void init()*/;
	/*static  void display()*/;
	//void DrawCell(const Cell& cell);
	const Mesh& m_mesh;

};

void DrawLine(const Vertex& point0, const Vertex& point1);
void DrawLine(const Vertex& point0, const Vertex& point1, const glm::vec3& color);
void DrawPolygon(const Vertex& point0, const Vertex& point1, const Vertex& point2, const Vertex& point3);
void DrawHexahedron(const Vertex& point0, const Vertex& point1, const Vertex& point2, const Vertex& point3,
					const Vertex& point4, const Vertex& point5, const Vertex& point6, const Vertex& point7);
void DrawHexahedron_S(const Vertex& point0, const Vertex& point1, const Vertex& point2, const Vertex& point3,
					const Vertex& point4, const Vertex& point5, const Vertex& point6, const Vertex& point7);
void DrawTetrahedron(const Vertex& point0, const Vertex& point1, const Vertex& point2, const Vertex& point3);
void DrawTriangle(const Vertex& point0, const Vertex& point1, const Vertex& point2);
void DrawTriangle(const Vertex& point0, const Vertex& point1, const Vertex& point2, const glm::vec3& color);
void DrawQuad(const Vertex& point0, const Vertex& point1, const Vertex& point2, const Vertex& point3, const glm::vec3& color);
void DrawCell(const Cell& cell);
void DrawWireframe(const Cell& cell);
void DrawWireFrameLine_S(const Vertex& point0, const Vertex& point1);
void DrawWireFrameLine_Yellow(const Vertex& point0, const Vertex& point1);

extern Mesh* pMesh;
extern Mesh* pTetMesh;
extern Mesh* pTetSurfaceMesh;
extern Mesh* pTriMesh;
extern Mesh* pTetMagnifiedMesh;
extern Mesh* pHexMagnifiedMesh;
extern PolycubeMesh* pTetPolycubeMesh;
extern PolycubeMesh* pTriPolycubeMesh;
extern Mesh* pTetPolycubeSurfaceMesh;
extern PolycubeMesh* pHexPolycubeMesh;
extern Mesh* pHexMesh;
extern Mesh* pResultMesh;
extern int SELECTED_CUBE_INDEX;
#endif /* __MESHDISPLAYER_H__ */
