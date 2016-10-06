/*
 * write_boundery.cpp
 *
 *  Created on: Dec 8, 2015
 *      Author: cotrik
 */

#include "include.h"

int main(int argc, char* argv[])
{
    if (argc != 3)
	{
		std::cout << "Usage: vtk2mesh <input_file> <output_file>" << std::endl;
		return -1;
	}

    MeshFileReader vtk(argv[1]);
    Mesh mesh(vtk.GetMesh());
    mesh.ExtractSurface();
    MeshFileWriter meshWriter(mesh, argv[2]);
    meshWriter.WriteMeshFile();
    meshWriter.WriteMeshBoundaryInfo(argv[2]);

	return 0;
}


