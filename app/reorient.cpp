/*
 * reorient.cpp
 *
 *  Created on: Oct 21, 2015
 *      Author: cotrik
 */
#include "include.h"

int main(int argc, char* argv[])
{
	if (argc != 3)
	{
		std::cout << "Usage: reorient <input_hex_file> <output_hex_file>" << std::endl;
		return -1;
	}

	MeshFileReader vtk(argv[1]);
	Mesh mesh(vtk.GetMesh());

	for (size_t i = 0; i < mesh.C.size(); i++)
	{
		Cell& c = mesh.C.at(i);
		if (!mesh.JudgeDirection(c))
		{
			swap(c[0], c[3]);
			swap(c[1], c[2]);
			swap(c[4], c[7]);
			swap(c[5], c[6]);
		}
	}

	MeshFileWriter meshWriter(mesh, argv[2]);
	meshWriter.WriteFile();

	return 0;
}


