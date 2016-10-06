/*
 * boundary.cpp
 *
 *  Created on: Jan 4, 2016
 *      Author: cotrik
 */

#include "include.h"

int main(int argc, char* argv[])
{
    if (argc != 3)
	{
		std::cout << "Usage: boudary <input_file> <output_file>" << std::endl;
		return -1;
	}

    MeshFileReader reader(argv[1]);
    Mesh mesh(reader.GetMesh());
    mesh.ExtractSurface();
    std::ofstream os(argv[2]);
    for (int i = 0; i < mesh.V.size(); i++)
        os << mesh.V.at(i).vinfo.bSurface << std::endl;

	return 0;
}
