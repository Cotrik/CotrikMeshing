/*
 * hexgen_s.cpp
 *
 *  Created on: Feb 3, 2016
 *      Author: cotrik
 */

#include "include.h"

int main(int argc, char* argv[])
{
	if (argc != 3)
	{
		std::cout << "Usage: hexgen_s <input.off> <output.vtk>" << std::endl;
		return -1;
	}
	MeshFileReader triMeshFileReader(argv[1]);
	PolycubeMesh triPolycubeMesh(triMeshFileReader.GetMesh());

	triPolycubeMesh.GetMaxMinCoordinates();
	triPolycubeMesh.ExtractSurface();
	triPolycubeMesh.GetFaceType();
	triPolycubeMesh.GetFaceAndNeighborFaces();
	triPolycubeMesh.GetFacePatches_N();
//	for (int i = 0; i < triPolycubeMesh.faceLabel.size(); i++)
//	{
//		std::cout << triPolycubeMesh.faceLabel.at(i) << std::endl;
//	}

	triPolycubeMesh.GetNormalOfSurfaceFaces();
	triPolycubeMesh.GetNormalOfSurfaceVertices();
	triPolycubeMesh.GetFacePatches();
	triPolycubeMesh.GetEdgePatches();
	triPolycubeMesh.GetVertexPatches();
	//////////////////////////////
	triPolycubeMesh.GetPatches();
	triPolycubeMesh.SortPatches();
	for (int i = 0; i < triPolycubeMesh.m_patches.size(); i++)
	{
		triPolycubeMesh.facePatches.at(i) = triPolycubeMesh.m_patches.at(i).m_f;
		for (int j = 0; j < triPolycubeMesh.facePatches.at(i).size(); j++)
		{
			triPolycubeMesh.faceLabel.at(triPolycubeMesh.facePatches.at(i).at(j)) = i + 1;
		}
		triPolycubeMesh.edgePatches.at(i) = triPolycubeMesh.m_patches.at(i).m_e;
		triPolycubeMesh.vertexPatches.at(i) = triPolycubeMesh.m_patches.at(i).m_v;
	}
	triPolycubeMesh.GetMaxMinCoordinates_Patch();
	//////////////////////////////

	triPolycubeMesh.GetPatchesLabel();
	triPolycubeMesh.GetFaceType();

	for (int i = 0; i < triPolycubeMesh.surface.size(); i++)
	{
		const Face& f = triPolycubeMesh.surface.at(i);
		//PFT.push_back(GetFaceType(mesh.N_F.at(i)));
		const Vertex* pVertex0 = &triPolycubeMesh.V[f.at(0)];
		const Vertex* pVertex1 = &triPolycubeMesh.V[f.at(1)];
		const Vertex* pVertex2 = &triPolycubeMesh.V[f.at(2)];

		const glm::vec3 v0(pVertex0->x, pVertex0->y, pVertex0->z);
		const glm::vec3 v1(pVertex1->x, pVertex1->y, pVertex1->z);
		const glm::vec3 v2(pVertex2->x, pVertex2->y, pVertex2->z);

		const glm::vec3 v10 = v0 - v1;
		const glm::vec3 v12 = v2 - v1;
		const glm::vec3 normal = glm::cross(v12, v10);

		//triPolycubeMesh.faceType.push_back(GetFaceType(normal));

		double area = 0.5f * glm::length(normal);
		triPolycubeMesh.faceArea.push_back(area);
	}


//	else if (mesh_type == MESH_TYPE_TRIANGLE_OFF)
//	{
//		//mesh.CheckPatchesConnection();
//		pTriMesh = (Mesh*)&triPolycubeMesh;
//		pTriPolycubeMesh = (PolycubeMesh*)&triPolycubeMesh;
//		DisplayMesh(argc, argv, triPolycubeMesh);
//	}
	return 0;
}
