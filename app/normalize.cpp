/*
 * normalize.cpp
 *
 *  Created on: Jan 6, 2016
 *      Author: cotrik
 */

#include "include.h"

int main(int argc, char* argv[])
{
    if (argc != 3)
	{
		std::cout << "Usage: normalize <input> <output>" << std::endl;
		return -1;
	}
    MeshFileReader reader(argv[1]);
    Mesh mesh(reader.GetMesh());
    mesh.NormalizeCoordinateValue();
    MeshFileWriter meshWriter(mesh, argv[2]);
    meshWriter.WriteFile();
	return 0;
}
