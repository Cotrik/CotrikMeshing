/*
 * hexgen_surface.cpp
 *
 *  Created on: Feb 18, 2016
 *      Author: cotrik
 */


#include "include.h"
#include "Parametrizer.h"
#include <math.h>

int main(int argc, char* argv[])
{
	std::string strInputOrigTriFileName;
	std::string strInputPolyTriFileName;
	std::string strInputPolyHexFileName;
	std::string strOutputHexFileName;
	float delta;
	try
	{
#ifdef CGAL_USE_BOOST_PROGRAM_OPTIONS
		po::options_description desc("Allowed options");
		desc.add_options()
			("help,h", "produce help message.")
			("thread,j", po::value<unsigned int>(&THREAD_NUM)->default_value(4), "threads number")
			// ########## input #########
			// -i
			("input-tri-file,i", po::value<std::string>(&strInputOrigTriFileName)->default_value("orig.tri.off"),
				"name of the orig input triangle off file")
			// -t
			("polycube-tri-file,t", po::value<std::string>(&strInputPolyTriFileName)->default_value("polycube.tri.vtk"),
				"name of the polycube input triangle vtk file")
			// -p
			("polycube-hex-file,p", po::value<std::string>(&strInputPolyHexFileName)->default_value("polycube.hex.vtk"),
				"name of the polycube input hexahedral vtk file")
			// ########## parameter #########
			// -s
			("hex-parameterization-size,s", po::value<float>(&delta)->default_value(0.02),
				"cube size, 0.001f <= size <= 0.1f")
			// ########## output #########
			// -o
			("output-hex-file,o", po::value<std::string>(&strOutputHexFileName)->default_value("hex.vtk"),
				"name of the output hexahedral vtk file")

			;

		po::variables_map vm;
		po::store(po::parse_command_line(argc, argv, desc), vm);
		po::notify(vm);

		if (vm.count("help"))
		{
			std::cout << desc << "\n";
			return 1;
		}

		if (vm.count("input-tri-file"))		std::cout << "strInputOrigTriFileName = " << strInputOrigTriFileName << std::endl;
		if (vm.count("polycube-tri-file"))		std::cout << "strInputPolyTriFileName = " << strInputPolyTriFileName << std::endl;
		if (vm.count("polycube-hex-file"))		std::cout << "strInputPolyHexFileName = " << strInputPolyHexFileName << std::endl;
		if (vm.count("output-hex-file"))		std::cout << "strOutputHexFileName = " << strOutputHexFileName << std::endl;
#else
		strInputOrigTriFileName = "orig.tri.off";
		strInputPolyTriFileName = "polycube.tri.off";
		strInputPolyHexFileName = "polycube.hex.vtk";
		strOutputHexFileName = "hex.vtk";
#endif
	} catch (std::exception& e)
	{
		std::cerr << "error: " << e.what() << "\n";
		return 1;
	} catch (...)
	{
		std::cerr << "Exception of unknown type!\n";
	}

	MeshFileReader orgTriMeshFileReader(strInputOrigTriFileName.c_str(), MESH_TYPE_TRIANGLE_OFF);
	Mesh origTriMesh(orgTriMeshFileReader.GetMesh());

	MeshFileReader polycubeTriMeshFileReader(strInputPolyTriFileName.c_str(), MESH_TYPE_TRIANGLE_VTK);
	PolycubeMesh polycubeTriMesh(polycubeTriMeshFileReader.GetMesh());
	////////
	polycubeTriMesh.GetMaxMinCoordinates();
	polycubeTriMesh.ExtractSurface();
	polycubeTriMesh.GetFaceType();
	polycubeTriMesh.GetFaceAndNeighborFaces();
	polycubeTriMesh.GetFacePatches_N();
	polycubeTriMesh.GetNormalOfSurfaceFaces();
	polycubeTriMesh.GetNormalOfSurfaceVertices();
	polycubeTriMesh.GetFacePatches();
	polycubeTriMesh.GetEdgePatches();
	polycubeTriMesh.GetVertexPatches();
	//////////////////////////////
	polycubeTriMesh.GetPatches();
	polycubeTriMesh.SortPatches();
	polycubeTriMesh.ModifyPatchesPosition();

	polycubeTriMesh.GetMaxMinCoordinates_Patch();
	polycubeTriMesh.AdjustPolycubeSurfaceToHoldAllSmallCubesVolume1(delta, polycubeTriMesh.m_minVertex);
	polycubeTriMesh.AddCompensation();

	MeshFileWriter polycubeTetMeshFileWriter(polycubeTriMesh.V, polycubeTriMesh.C, "polycube.tri.n.a.vtk", MESH_TYPE_TRIANGLE_VTK);
	polycubeTetMeshFileWriter.WriteMeshFile("polycube.tri.n.a.vtk", MESH_TYPE_TRIANGLE_VTK);
	////////

	MeshFileReader polycubeHexMeshFileReader(strInputPolyHexFileName.c_str(), MESH_TYPE_HEXAHEDRON_VTK);
	PolycubeMesh polycubeHexMesh(polycubeHexMeshFileReader.GetMesh());
	/////////////////////////
	polycubeHexMesh.GetMaxMinCoordinates();
	polycubeHexMesh.ExtractSurface();
	polycubeHexMesh.GetFaceType();
	polycubeHexMesh.GetFaceAndNeighborFaces();
	polycubeHexMesh.GetFacePatches_N();
	polycubeHexMesh.GetNormalOfSurfaceFaces();
	polycubeHexMesh.GetNormalOfSurfaceVertices();
	polycubeHexMesh.GetFacePatches();
	polycubeHexMesh.GetEdgePatches();
	polycubeHexMesh.GetVertexPatches();

	polycubeHexMesh.GetPatches();
	polycubeHexMesh.SortPatches();
	polycubeHexMesh.ModifyPatchesPosition();
	/////////////////////////
	std::vector<Vertex> hexV(polycubeHexMesh.V);

	polycubeTriMesh.GetPatches();
	polycubeTriMesh.SortPatches();
	polycubeTriMesh.ModifyPatchesPosition();

	polycubeHexMesh.GetPatches();
	polycubeHexMesh.SortPatches();
	polycubeHexMesh.ModifyPatchesPosition();

	Parametrizer parametrizer;
	for (int i = 0; i < polycubeHexMesh.V.size(); i++)
	{
		std::vector<glm::vec3> w;
		parametrizer.Parametrize3DInnerPoint(polycubeHexMesh.V.at(i), polycubeTriMesh.V, polycubeTriMesh.C, w, 1e-6);
		double total_w = 0;
		glm::vec3 total_p(0.0, 0.0, 0.0);
		for (int j = 0; j < polycubeTriMesh.C.size(); j++)
		{
			const Cell& c = polycubeTriMesh.C.at(j);
			const double wi[3] = {w[j].x, w[j].y, w[j].z};
			for (int k = 0; k < 3; k++)
			{
				total_w += wi[k];
				total_p.x += wi[k] * origTriMesh.V.at(c.at(k)).x;
				total_p.y += wi[k] * origTriMesh.V.at(c.at(k)).y;
				total_p.z += wi[k] * origTriMesh.V.at(c.at(k)).z;
			}
		}
		const Vertex v(total_p.x/total_w, total_p.y/total_w, total_p.z/total_w);
		hexV.at(i).x = v.x;
		hexV.at(i).y = v.y;
		hexV.at(i).z = v.z;
	}
	MeshFileWriter origHexMeshFileWriter(hexV, polycubeHexMesh.C, strOutputHexFileName.c_str(), MESH_TYPE_HEXAHEDRON_VTK);
	origHexMeshFileWriter.WriteMeshFile(strOutputHexFileName.c_str(), MESH_TYPE_HEXAHEDRON_VTK);

	return 0;
}
