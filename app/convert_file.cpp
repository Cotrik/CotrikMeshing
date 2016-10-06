/*
 * convert_file.cpp
 *
 *  Created on: Sep 30, 2015
 *      Author: cotrik
 */


#include "include.h"

int main(int argc, char* argv[])
{
	if (argc != 3)
	{
		std::cout << "Usage: convert_file <vtk2mesh | vtk2off | mesh2vtk | mesh2off | off2mesh | off2vtk | off2obj | off2stl | off2ply | obj2vtk | obj2mesh | obj2off | obj2stl | obj2ply | stl2vtk | stl2mesh | stl2off | stl2obj | stl2ply> <input_file> <output_file>" << std::endl;
		return -1;
	}
	MeshFileReader reader(argv[1]);
	MeshFileWriter writer(reader.GetMesh(), argv[2]);
	writer.WriteFile();
	return 0;
}

