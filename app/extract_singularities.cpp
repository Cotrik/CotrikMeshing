/*
 * extract_surface.cpp
 *
 *  Created on: Feb 3, 2016
 *      Author: cotrik
 */

#include "include.h"

int main(int argc, char* argv[])
{
    if (argc != 3)
    {
        std::cout << "Usage: singularities <input_file> <output_vtk_file>" << std::endl;
        return -1;
    }

    MeshFileReader reader(argv[1]);
    Mesh mesh(reader.GetMesh());
    mesh.GetSingularities();
    MeshFileWriter writer(mesh, argv[2]);
    writer.WriteSingularitiesVtkFile();

	return 0;
}
