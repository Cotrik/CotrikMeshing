/*
 * mean_value_corrdinates_parametrizaton.cpp
 *
 *  Created on: May 20, 2016
 *      Author: cotrik
 */

#include "include.h"
#include "Parametrizer.h"

int main(int argc, char* argv[])
{
	if (argc != 4)
	{
		std::cout << "Usage: mean_value_coordinates_parametrization <input_mesh_file> <input_tri_file> <out_put_scalar_field>" << std::endl;
		return -1;
	}

	MeshFileReader origMeshFileReader(argv[1]);
	PolycubeMesh origMesh(origMeshFileReader.GetMesh());
	origMesh.ExtractSurface();
    std::vector<double> scalars(origMesh.V.size());
    for (int i = 0; i < origMesh.V.size(); i++)
    {
        if (origMesh.V.at(i).vinfo.bSurface)
            scalars.at(i) = 1.0;
    }
	MeshFileReader orgTriMeshFileReader(argv[2]);
	Mesh origTriMesh(orgTriMeshFileReader.GetMesh());

	Parametrizer parametrizer;
	std::ofstream ofs(argv[3]);
	ofs << "POINT_DATA " << origMesh.V.size() << std::endl;
	ofs << "SCALARS para float 1\nLOOKUP_TABLE default\n";
	for (int i = 0; i < origMesh.V.size(); i++)
	{
		std::vector<glm::vec3> w;
		parametrizer.Parametrize3DInnerPoint(origMesh.V.at(i), origTriMesh.V, origTriMesh.C, w, 1e-6);
		double total_w = 0;
		glm::vec3 total_p(0.0, 0.0, 0.0);
		double total_s = 0;
		for (int j = 0; j < origTriMesh.C.size(); j++)
		{
			const Cell& c = origTriMesh.C.at(j);
			const double wi[3] = {w[j].x, w[j].y, w[j].z};
			for (int k = 0; k < 3; k++)
			{
				total_w += wi[k];
				total_s += wi[k] * scalars.at(c.at(k));
			}
		}
		float s = total_s/total_w;
        if (origMesh.V.at(i).vinfo.bSurface)
            s = 1;
		ofs << s << std::endl;
	}

	return 0;
}
