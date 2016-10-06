/*
 * display_mesh.cpp
 *
 *  Created on: Jun 13, 2015
 *      Author: cotrik
 */

#include "include.h"

int main(int argc, char* argv[])
{
	if (argc != 2)
	{
		std::cout << "Usage: display_mesh mesh_filename!" << std::endl;
		return -1;
	}

	std::string strInputFileName(argv[1]);
	MESH_TYPE mesh_type = GetMeshType(strInputFileName.c_str());
	MeshFileReader meshFileReader(strInputFileName.c_str(), mesh_type);
	Mesh mesh(meshFileReader.GetMesh());
	mesh.ExtractSurface();

	if (mesh_type == MESH_TYPE_HEXAHEDRON_VTK || mesh_type == MESH_TYPE_HEXAHEDRON)
	{
		pHexMesh         = (Mesh*)&mesh;
		DisplayMesh(argc, argv, mesh);
	}
	else if (mesh_type == MESH_TYPE_TETRAHEDRON_VTK)
	{
		pTetMesh         = (Mesh*)&mesh;
		DisplayMesh(argc, argv, mesh);
	}
	//	pTetMesh         = (Mesh*)&orgTetMeshFileReader.GetMesh();
	//	pTetMesh         = (Mesh*)&origTetNormalizedMesh;
	//	pTetSurfaceMesh  = (Mesh*)&tetSurfaceMesh;
	//	pTetPolycubeMesh = (Mesh*)& ;
	//	pHexPolycubeMesh = (Mesh*)&;
	//	pHexMesh         = (Mesh*)&hexMesh;
	//	pHexMesh         = (Mesh*)&hexMesh;
	//	pResultMesh      = (Mesh*)&origTetMesh;
	//	DisplayMesh(argc, argv, origTetMesh);
	return 0;
}


