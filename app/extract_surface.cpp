/*
 * extract_surface.cpp
 *
 *  Created on: Feb 3, 2016
 *      Author: cotrik
 */

#include "include.h"

int main(int argc, char* argv[])
{
    if (argc != 3 && argc != 4)
    {
        std::cout << "Usage: extract_surface <input> <output> [1]" << std::endl;
        return -1;
    }

    MeshFileReader reader(argv[1]);
    Mesh mesh(reader.GetMesh());
    mesh.ExtractSurface();
    ElementType surfaceCellType = TRIANGLE;
    if (mesh.m_cellType == HEXAHEDRA)    surfaceCellType = QUAD;
    MeshFileWriter writer(mesh.V, mesh.surface, argv[2], surfaceCellType);
    if (argc == 4)
        writer.SetFixFlag();
    writer.WriteFile();

	return 0;
}
