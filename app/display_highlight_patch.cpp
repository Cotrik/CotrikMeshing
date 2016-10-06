/*
 * display_highlight_patch.cpp
 *
 *  Created on: Jun 15, 2015
 *      Author: cotrik
 */

#include "include.h"

int main(int argc, char* argv[])
{
	if (argc != 2)
	{
		std::cout << "Usage: display_highlight_patch mesh_filename!" << std::endl;
		return -1;
	}

	std::string strInputFileName(argv[1]);
	MeshFileReader meshFileReader(strInputFileName.c_str());
	PolycubeMesh mesh(meshFileReader.GetMesh());
	mesh.ExtractSurface();

	////////////////////////////////////////////////////////////////////////////////
	mesh.NormalizeCoordinateValue();
	mesh.ExtractSurface();
	mesh.GetNormalOfSurfaceFaces();
	mesh.GetNormalOfSurfaceVertices();
	mesh.GetCorners();
	mesh.GetVertexInfo();

	mesh.LabelSurfaceFace();
	mesh.GetFacePatches();
	mesh.GetEdgePatches();
	mesh.GetVertexPatches();
	mesh.GetPatchesLabel();
	mesh.GetFaceType();

	for (int i = 0; i < mesh.surface.size(); i++)
	{
		const Face& f = mesh.surface.at(i);
		//PFT.push_back(GetFaceType(mesh.N_F.at(i)));
		const Vertex* pVertex0 = &mesh.V[f.at(0)];
		const Vertex* pVertex1 = &mesh.V[f.at(1)];
		const Vertex* pVertex2 = &mesh.V[f.at(2)];

		const glm::vec3 v0(pVertex0->x, pVertex0->y, pVertex0->z);
		const glm::vec3 v1(pVertex1->x, pVertex1->y, pVertex1->z);
		const glm::vec3 v2(pVertex2->x, pVertex2->y, pVertex2->z);

		const glm::vec3 v10 = v0 - v1;
		const glm::vec3 v12 = v2 - v1;
		const glm::vec3 normal = glm::cross(v12, v10);

		//mesh.faceType.push_back(GetFaceType(normal));

		double area = 0.5f * glm::length(normal);
		mesh.faceArea.push_back(area);
	}

	if (mesh.m_cellType == HEXAHEDRA)
	{
		mesh.CheckPatchesConnection();

		pHexMesh         = (Mesh*)&mesh;
		pHexPolycubeMesh = (PolycubeMesh*)&mesh;
		DisplayMesh(argc, argv, mesh);
	}
	else if (mesh.m_cellType == TETRAHEDRA)
	{
		mesh.CheckPatchesConnection();

		pTetMesh         = (Mesh*)&mesh;
		pTetPolycubeMesh = (PolycubeMesh*)&mesh;
		DisplayMesh(argc, argv, mesh);
	}
	else if (mesh.m_cellType == TRIANGLE)
	{
		mesh.CheckPatchesConnection();

		pTriMesh         = (Mesh*)&mesh;
		pTriPolycubeMesh = (PolycubeMesh*)&mesh;
		DisplayMesh(argc, argv, mesh);
	}
	return 0;
}


