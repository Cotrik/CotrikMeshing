/*
 * cone_optimize.cpp
 *
 *  Created on: Sep 28, 2015
 *      Author: cotrik
 */

#include "include.h"

int main(int argc, char* argv[])
{
	if (argc < 3)
	{
		std::cout << "Usage: cone_opt <input.vtk> <output.vtk> <iter_max>" << std::endl;
		return -1;
	}

	int iter_max = 1;
	if (argc == 4)
	{
		std::stringstream ss(argv[3]);
		ss >> iter_max;
	}

	MeshFileReader vtk(argv[1]);
	Mesh mesh(vtk.GetMesh());

	mesh.ExtractSurface();
	mesh.GetConeConfiguration();
	mesh.OptimizeConeShape(iter_max);

	MeshFileWriter meshWriter(mesh, argv[2]);
	meshWriter.WriteMeshFile();

	return 0;
}
