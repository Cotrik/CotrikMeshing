/*
 * verify_structure.cpp
 *
 *  Created on: Oct 1, 2015
 *      Author: cotrik
 */

#include "include.h"

int main(int argc, char* argv[])
{
	if (argc != 3)
	{
		std::cout << "Usage: verify_structure input_tet_meshfile verify_hex_mesh_file" << std::endl;
		return -1;
	}

	MeshFileReader tetMeshReader(argv[1]);
	Mesh tetMesh(tetMeshReader.GetMesh());
	int genus = tetMesh.GetGenus();
	std::cout << "tet mesh genus = " << genus << std::endl;
	MeshFileReader meshReader(argv[2]);
	Mesh mesh(meshReader.GetMesh());
	if (mesh.VerifyStructure(genus))
	{
		std::cout << "structure of " << argv[2] << " is correct!" << std::endl;
	}

	if (!mesh.VerifyElements())
	{
		std::cout << "Elements of " << argv[2] << " is not correct!" << std::endl;
	}
	return 0;
}
