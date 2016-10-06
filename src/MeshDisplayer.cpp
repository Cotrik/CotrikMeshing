/*
 * MeshDisplayer.cpp
 *
 *  Created on: Feb 13, 2015
 *      Author: cotrik
 */

#include "MeshDisplayer.h"
// Include standard headers
#include <stdio.h>
#include <stdlib.h>
#include <vector>

//// Include GLEW
#include <GL/glew.h>
//
//// Include GLFW
#include <glfw3.h>

void DrawLine(const Vertex& point0, const Vertex& point1)
{
//	glLineWidth(2.0f);
//	glColor3f(0.0f,0.0f,0.5f);
	glBegin(GL_LINES);
	glVertex3f((GLfloat)point0.x, (GLfloat)point0.y, (GLfloat)point0.z);
	glVertex3f((GLfloat)point1.x, (GLfloat)point1.y, (GLfloat)point1.z);
	glEnd();
}

void DrawLine(const Vertex& point0, const Vertex& point1, const glm::vec3& color)
{
	glLineWidth(2.0f);
	glColor3f(color.r, color.g, color.b);
	glBegin(GL_LINES);
	glVertex3f((GLfloat)point0.x, (GLfloat)point0.y, (GLfloat)point0.z);
	glVertex3f((GLfloat)point1.x, (GLfloat)point1.y, (GLfloat)point1.z);
	glEnd();
}

void DrawWireFrameLine(const Vertex& point0, const Vertex& point1)
{
	glLineWidth(1.0f);
	if (pMesh->GetCubeIndex(point0, 0.08) == SELECTED_CUBE_INDEX
		|| pMesh->GetCubeIndex(point1, 0.08) == SELECTED_CUBE_INDEX)
		glColor3f(1.0f,0.0f,0.0f);
	else
		glColor3f(1.0f,1.0f,1.0f);
	glBegin(GL_LINES);
	glVertex3f((GLfloat)point0.x, (GLfloat)point0.y, (GLfloat)point0.z);
	glVertex3f((GLfloat)point1.x, (GLfloat)point1.y, (GLfloat)point1.z);
	glEnd();
}

void DrawWireFrameLine_S(const Vertex& point0, const Vertex& point1)
{
	glLineWidth(1.0f);
	glColor3f(0.0f,1.0f,0.0f);
	glBegin(GL_LINES);
	glVertex3f((GLfloat)point0.x, (GLfloat)point0.y, (GLfloat)point0.z);
	glVertex3f((GLfloat)point1.x, (GLfloat)point1.y, (GLfloat)point1.z);
	glEnd();
}

void DrawWireFrameLine_Yellow(const Vertex& point0, const Vertex& point1)
{
	glLineWidth(1.0f);
	glColor3f(1.0f,1.0f,0.0f);
	glBegin(GL_LINES);
	glVertex3f((GLfloat)point0.x, (GLfloat)point0.y, (GLfloat)point0.z);
	glVertex3f((GLfloat)point1.x, (GLfloat)point1.y, (GLfloat)point1.z);
	glEnd();
}

void DrawDashLine(const Vertex& point0, const Vertex& point1)
{
	glEnable(GL_LINE_STIPPLE);
	glLineStipple(2, 0x0F0F);
	glLineWidth(1.0f);
	glColor3f(0.1f, 0.1f, 0.1f);
	glBegin(GL_LINES);
	glVertex3f((GLfloat) point0.x, (GLfloat) point0.y, (GLfloat) point0.z);
	glVertex3f((GLfloat) point1.x, (GLfloat) point1.y, (GLfloat) point1.z);
	glEnd();
	glDisable(GL_LINE_STIPPLE);
}

void DrawPolygon(const Vertex& point0, const Vertex& point1, const Vertex& point2, const Vertex& point3)
{
	glBegin(GL_POLYGON); // Draw A Quad
	glVertex3f((GLfloat)point0.x, (GLfloat)point0.y, (GLfloat)point0.z); // Top Left
	glVertex3f((GLfloat)point1.x, (GLfloat)point1.y, (GLfloat)point1.z); // Top Right
	glVertex3f((GLfloat)point2.x, (GLfloat)point2.y, (GLfloat)point2.z); // Bottom Right
	glVertex3f((GLfloat)point3.x, (GLfloat)point3.y, (GLfloat)point3.z); // Bottom Left
	glVertex3f((GLfloat)point0.x, (GLfloat)point0.y, (GLfloat)point0.z); // Top Left
	glEnd();
};

void DrawHexahedron(const Vertex& point0, const Vertex& point1, const Vertex& point2, const Vertex& point3,
					const Vertex& point4, const Vertex& point5, const Vertex& point6, const Vertex& point7)
{
	glColor3f(0.5f, 0.5f, 0.5f);
//	glColor3f(1.0f, 0.0f, 0.0f); // Set Right Point Of Triangle To Red
	DrawPolygon(point0, point1, point2, point3);
//	glColor3f(1.0f, 1.0f, 0.0f); // Set Left Point Of Triangle To Yellow
	DrawPolygon(point0, point4, point5, point1);
//	glColor3f(0.0f, 1.0f, 0.0f); // Set Top Point Of Triangle To Green
	DrawPolygon(point0, point4, point7, point3);
//	glColor3f(1.0f, 1.0f, 0.0f); // Set Left Point Of Triangle To Yellow
	DrawPolygon(point3, point7, point6, point2);
//	glColor3f(0.0f, 1.0f, 0.0f); // Set Top Point Of Triangle To Green
	DrawPolygon(point1, point5, point6, point2);
//	glColor3f(1.0f, 0.0f, 0.0f); // Set Right Point Of Triangle To Red
	DrawPolygon(point4, point5, point6, point7);

	glLineWidth(2.0f);
	glColor3f(0.0f,0.0f,0.5f);
	DrawLine(point0, point1);
	DrawLine(point1, point2);
	DrawLine(point2, point3);
	DrawLine(point3, point0);
	DrawLine(point4, point5);
	DrawLine(point5, point6);
	DrawLine(point6, point7);
	DrawLine(point7, point4);
	DrawLine(point0, point4);
	DrawLine(point1, point5);
	DrawLine(point2, point6);
	DrawLine(point3, point7);
};

void DrawHexahedron_S(const Vertex& point0, const Vertex& point1, const Vertex& point2, const Vertex& point3,
					const Vertex& point4, const Vertex& point5, const Vertex& point6, const Vertex& point7)
{
	glColor3f(0.5f, 0.5f, 0.5f);
	glColor3f(1.0f, 0.0f, 0.0f); // Set Right Point Of Triangle To Red
	DrawPolygon(point0, point1, point2, point3);
	glColor3f(1.0f, 1.0f, 0.0f); // Set Left Point Of Triangle To Yellow
	DrawPolygon(point0, point4, point5, point1);
	glColor3f(0.0f, 1.0f, 0.0f); // Set Top Point Of Triangle To Green
	DrawPolygon(point0, point4, point7, point3);
	glColor3f(1.0f, 1.0f, 0.0f); // Set Left Point Of Triangle To Yellow
	DrawPolygon(point3, point7, point6, point2);
	glColor3f(0.0f, 1.0f, 0.0f); // Set Top Point Of Triangle To Green
	DrawPolygon(point1, point5, point6, point2);
	glColor3f(1.0f, 0.0f, 0.0f); // Set Right Point Of Triangle To Red
	DrawPolygon(point4, point5, point6, point7);

	glColor3f(0.5f,0.1f,0.7f);
	DrawLine(point0, point1);
	DrawLine(point1, point2);
	DrawLine(point2, point3);
	DrawLine(point3, point0);
	DrawLine(point4, point5);
	DrawLine(point5, point6);
	DrawLine(point6, point7);
	DrawLine(point7, point4);
	DrawLine(point0, point4);
	DrawLine(point1, point5);
	DrawLine(point2, point6);
	DrawLine(point3, point7);
};
void DrawTetrahedron(const Vertex& point0, const Vertex& point1, const Vertex& point2, const Vertex& point3)
{
//	glColor3f(1.0f,0.0f,0.0f); // Set Top Point Of Triangle To Red
	if (pMesh->GetCubeIndex(point0, 0.08) == SELECTED_CUBE_INDEX
		|| pMesh->GetCubeIndex(point1, 0.08) == SELECTED_CUBE_INDEX
		|| pMesh->GetCubeIndex(point2, 0.08) == SELECTED_CUBE_INDEX
		|| pMesh->GetCubeIndex(point3, 0.08) == SELECTED_CUBE_INDEX)
		glColor3f(1.0f,0.0f,0.0f);
	else
	glColor3f(0.5f,0.5f,0.5f);
	glBegin(GL_TRIANGLES); // Draw A Triangle
	glVertex3f((GLfloat)point0.x, (GLfloat)point0.y, (GLfloat)point0.z);
	glVertex3f((GLfloat)point1.x, (GLfloat)point1.y, (GLfloat)point1.z);
	glVertex3f((GLfloat)point3.x, (GLfloat)point3.y, (GLfloat)point3.z);
	glEnd();

//	glColor3f(0.0f,1.0f,0.0f); // Set Left Point Of Triangle To Green
	glBegin(GL_TRIANGLES); // Draw A Triangle
	glVertex3f((GLfloat)point0.x, (GLfloat)point0.y, (GLfloat)point0.z);
	glVertex3f((GLfloat)point2.x, (GLfloat)point2.y, (GLfloat)point2.z);
	glVertex3f((GLfloat)point3.x, (GLfloat)point3.y, (GLfloat)point3.z);
	glEnd();

//	glColor3f(0.0f,0.0f,1.0f); // Set Right Point Of Triangle To Blue
	glBegin(GL_TRIANGLES); // Draw A Triangle
	glVertex3f((GLfloat)point0.x, (GLfloat)point0.y, (GLfloat)point0.z);
	glVertex3f((GLfloat)point1.x, (GLfloat)point1.y, (GLfloat)point1.z);
	glVertex3f((GLfloat)point2.x, (GLfloat)point2.y, (GLfloat)point2.z);
	glEnd();

	glBegin(GL_TRIANGLES); // Draw A Triangle
	glVertex3f((GLfloat)point1.x, (GLfloat)point1.y, (GLfloat)point1.z);
	glVertex3f((GLfloat)point2.x, (GLfloat)point2.y, (GLfloat)point2.z);
	glVertex3f((GLfloat)point3.x, (GLfloat)point3.y, (GLfloat)point3.z);
	glEnd();

	glLineWidth(2.0f);
	glColor3f(0.0f,0.0f,0.5f);
	DrawLine(point0, point1);
	DrawLine(point0, point2);
	DrawLine(point0, point3);
	DrawLine(point1, point2);
	DrawLine(point1, point3);
	DrawLine(point2, point3);
};

void DrawTriangle(const Vertex& point0, const Vertex& point1, const Vertex& point2)
{
	glColor3f(0.0f,0.0f,1.0f); // Set Right Point Of Triangle To Blue
	glBegin(GL_TRIANGLES); // Draw A Triangle
	glVertex3f((GLfloat)point0.x, (GLfloat)point0.y, (GLfloat)point0.z);
	glVertex3f((GLfloat)point1.x, (GLfloat)point1.y, (GLfloat)point1.z);
	glVertex3f((GLfloat)point2.x, (GLfloat)point2.y, (GLfloat)point2.z);
	glEnd();
};

void DrawTriangle(const Vertex& point0, const Vertex& point1, const Vertex& point2, const glm::vec3& color)
{
	glColor3f(color.r,color.g,color.b); // Set Right Point Of Triangle To Blue
	glBegin(GL_TRIANGLES); // Draw A Triangle
	glVertex3f((GLfloat)point0.x, (GLfloat)point0.y, (GLfloat)point0.z);
	glVertex3f((GLfloat)point1.x, (GLfloat)point1.y, (GLfloat)point1.z);
	glVertex3f((GLfloat)point2.x, (GLfloat)point2.y, (GLfloat)point2.z);
	glEnd();
};

void DrawQuad(const Vertex& point0, const Vertex& point1, const Vertex& point2, const Vertex& point3, const glm::vec3& color)
{
	glColor3f(color.r,color.g,color.b); // Set Right Point Of Triangle To Blue
	glBegin(GL_QUADS); // Draw A Triangle
	glVertex3f((GLfloat)point0.x, (GLfloat)point0.y, (GLfloat)point0.z);
	glVertex3f((GLfloat)point1.x, (GLfloat)point1.y, (GLfloat)point1.z);
	glVertex3f((GLfloat)point2.x, (GLfloat)point2.y, (GLfloat)point2.z);
	glVertex3f((GLfloat)point3.x, (GLfloat)point3.y, (GLfloat)point3.z);
	glEnd();
};

Mesh* pMesh = NULL;
Mesh* pTetMesh = NULL;
Mesh* pTriMesh = NULL;
Mesh* pTetSurfaceMesh = NULL;
PolycubeMesh* pTetPolycubeMesh = NULL;
PolycubeMesh* pTriPolycubeMesh = NULL;
Mesh* pTetPolycubeSurfaceMesh = NULL;
PolycubeMesh* pHexPolycubeMesh = NULL;
Mesh* pHexMesh = NULL;
Mesh* pResultMesh = NULL;
Mesh* pTetMagnifiedMesh;
Mesh* pHexMagnifiedMesh;
int SELECTED_CUBE_INDEX = 0;
void DrawCell(const Cell& cell)
{
	const size_t cellSize = cell.size();
	if (cellSize > 9)
		return;
	const std::vector<Vertex>& v = pMesh->V;
	if (cellSize == 3)
		DrawTriangle(v[cell[0]], v[cell[1]], v[cell[2]]);
	else if (cellSize == 4)
	{
		int count = 0;
		if (pMesh->vertexInfo.at(cell[0]).bSurface) count++;
		if (pMesh->vertexInfo.at(cell[1]).bSurface) count++;
		if (pMesh->vertexInfo.at(cell[2]).bSurface) count++;
		if (pMesh->vertexInfo.at(cell[3]).bSurface) count++;
		if (count >= 3)
			DrawTetrahedron(v[cell[0]], v[cell[1]], v[cell[2]], v[cell[3]]);
	}
	else if (cellSize == 8)
	{
		if (pMesh->vertexInfo.at(cell[0]).bSurface
			|| pMesh->vertexInfo.at(cell[1]).bSurface
			|| pMesh->vertexInfo.at(cell[2]).bSurface
			|| pMesh->vertexInfo.at(cell[3]).bSurface
			|| pMesh->vertexInfo.at(cell[4]).bSurface
			|| pMesh->vertexInfo.at(cell[5]).bSurface
			|| pMesh->vertexInfo.at(cell[6]).bSurface
			|| pMesh->vertexInfo.at(cell[7]).bSurface)
		DrawHexahedron(v[cell[0]], v[cell[1]], v[cell[2]], v[cell[3]], v[cell[4]], v[cell[5]], v[cell[6]], v[cell[7]]);
	}
	else if (cellSize == 9)
	{
		for (unsigned int i = 0; i < 12; i++)
		{
			DrawDashLine(v[cell[HexEdge[i][0]]], v[cell[HexEdge[i][1]]]);
		}
	}
}

void DrawWireframe(const Cell& cell)
{
	const size_t cellSize = cell.size();
	if (cellSize > 9)
		return;
	const std::vector<Vertex>& v = pMesh->V;
	if (cellSize == 3)
	{
		DrawWireFrameLine(v[cell[0]], v[cell[1]]);
		DrawWireFrameLine(v[cell[1]], v[cell[2]]);
		DrawWireFrameLine(v[cell[2]], v[cell[0]]);
	}
	else if (cellSize == 4)
	{
		const std::vector<VertexInfo>& vInfo = pMesh->vertexInfo;
		for (unsigned int i = 0; i < 6; i++)
		{
			if (vInfo[cell[TetEdge[i][0]]].bSurface && vInfo[cell[TetEdge[i][1]]].bSurface)
			DrawWireFrameLine(v[cell[TetEdge[i][0]]], v[cell[TetEdge[i][1]]]);
		}
	}
	else if (cellSize == 8)
	{
		const std::vector<VertexInfo>& vInfo = pMesh->vertexInfo;
		for (unsigned int i = 0; i < 12; i++)
		{
			if (vInfo[cell[HexEdge[i][0]]].bSurface && vInfo[cell[HexEdge[i][1]]].bSurface)
			DrawWireFrameLine(v[cell[HexEdge[i][0]]], v[cell[HexEdge[i][1]]]);
		}
//		const Cell& c = pHexMesh->C.at(SELECTED_CUBE_INDEX);
//		for (unsigned int i = 0; i < 12; i++)
//		{
//			DrawWireFrameLine_S(v[c[HexEdge[i][0]]], v[c[HexEdge[i][1]]]);
//		}
	}
	else if (cellSize == 9)
	{
		for (unsigned int i = 0; i < 12; i++)
		{
			DrawDashLine(v[cell[HexEdge[i][0]]], v[cell[HexEdge[i][1]]]);
		}
	}
}

MeshDisplayer::MeshDisplayer(const Mesh& mesh)
: m_mesh(mesh)
{
	// TODO Auto-generated constructor stub
	pMesh = (Mesh *)&m_mesh;
}

MeshDisplayer::~MeshDisplayer()
{
	// TODO Auto-generated destructor stub
}
