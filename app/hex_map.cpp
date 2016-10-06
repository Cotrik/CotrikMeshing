/*
 * hex_map.cpp
 *
 *  Created on: May 21, 2016
 *      Author: cotrik
 */

#include "include.h"
#include "Parametrizer.h"

int main(int argc, char* argv[])
{
	if (argc != 5)
	{
		std::cout << "Usage: hex_map <input_hex_file> <input_tri_file> <input_deformed_tri_file> <output_hex_file>" << std::endl;
		return -1;
	}

	MeshFileReader orgHexMeshFileReader(argv[1]);
	Mesh origHexMesh(orgHexMeshFileReader.GetMesh());

	MeshFileReader orgTriMeshFileReader(argv[2]);
	Mesh origTriMesh(orgTriMeshFileReader.GetMesh());

	MeshFileReader deformedTriMeshFileReader(argv[3]);
	Mesh deformedMesh(deformedTriMeshFileReader.GetMesh());

	Parametrizer parametrizer;
	std::vector<Vertex> deformedHexV(origHexMesh.V);
	for (int i = 0; i < origHexMesh.V.size(); i++)
	{
		std::vector<glm::vec3> w;
		parametrizer.Parametrize3DInnerPoint(origHexMesh.V.at(i), origTriMesh.V, origTriMesh.C, w, 1e-6);
		double total_w = 0;
		glm::vec3 total_p(0.0, 0.0, 0.0);
		for (int j = 0; j < origTriMesh.C.size(); j++)
		{
			const Cell& c = origTriMesh.C.at(j);
			const double wi[3] = {w[j].x, w[j].y, w[j].z};
			for (int k = 0; k < 3; k++)
			{
				total_w += wi[k];
				total_p.x += wi[k] * deformedMesh.V.at(c.at(k)).x;
				total_p.y += wi[k] * deformedMesh.V.at(c.at(k)).y;
				total_p.z += wi[k] * deformedMesh.V.at(c.at(k)).z;
			}
		}
		const Vertex v(total_p.x/total_w, total_p.y/total_w, total_p.z/total_w);
		deformedHexV.at(i).x = v.x;
		deformedHexV.at(i).y = v.y;
		deformedHexV.at(i).z = v.z;
		if (isnan(v.x))
			std::cout << "vertex " << i << " is nan" << std::endl;
	}

	MeshFileWriter deformedHexMeshFileWriter(deformedHexV, origHexMesh.C, argv[4]);
	deformedHexMeshFileWriter.WriteFile();

	return 0;
}



