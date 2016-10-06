/*
 * ajust_polycube.cpp
 *
 *  Created on: Aug 4, 2015
 *      Author: cotrik
 */
#include "include.h"

int AjustPolycubeSurface(int argc, char* argv[]);

int main(int argc, char* argv[])
{
	return AjustPolycubeSurface(argc, argv);
}

int AjustPolycubeSurface(int argc, char* argv[])
{
	std::string strOutputResultHexFileName;
	std::string strInputFileName;
	float discrete = 0.001f;
	try
	{
#ifdef CGAL_USE_BOOST_PROGRAM_OPTIONS
		po::options_description desc("Allowed options");
		desc.add_options()
			("help,h", "produce help message.")
			// ########## input #########
			// -i
			("input-polycube-file,i", po::value<std::string>(&strInputFileName)->default_value("polycube_aligned2.hex.vtk"),
				"name of the polycube input hexahedral vtk file")
			// -o
			("out-polycube-file,o", po::value<std::string>(&strOutputResultHexFileName)->default_value("polycube_ajusted.hex.vtk"),
				"name of the output polycube hexahedral vtk/mesh file")
			// -d
			("discrete,d", po::value<float>(&discrete)->default_value(0.001), "discrete length")
			;

		po::variables_map vm;
		po::store(po::parse_command_line(argc, argv, desc), vm);
		po::notify(vm);

		if (vm.count("help"))
		{
			std::cout << desc << "\n";
			return 1;
		}

		std::cout << "strInputFileName = " << strInputFileName << std::endl;
		std::cout << "strOutputResultHexFileName = " << strOutputResultHexFileName << std::endl;
#else
		strInputOrigTetFileName = "orig.tet.vtk";
		strInputPolyTetFileName = "polycube.tet.vtk";
		strOutputPolyHexFileName = "poly.hex.vtk";
		strOutputResultHexFileName = "result.hex.vtk";
#endif
	} catch (std::exception& e)
	{
		std::cerr << "error: " << e.what() << "\n";
		return 1;
	} catch (...)
	{
		std::cerr << "Exception of unknown type!\n";
	}

	MeshFileReader polycubeMeshFileReader(strInputFileName.c_str(), MESH_TYPE_TETRAHEDRON_VTK);
	PolycubeMesh polycubeMesh(polycubeMeshFileReader.GetMesh());

	polycubeMesh.ExtractSurface();
	polycubeMesh.NormalizeCoordinateValue();
	polycubeMesh.ExtractSurface();
	polycubeMesh.GetNormalOfSurfaceFaces();
	polycubeMesh.GetNormalOfSurfaceVertices();
	polycubeMesh.GetCorners();
	polycubeMesh.GetVertexInfo();

	polycubeMesh.LabelSurfaceFace();
	polycubeMesh.GetFacePatches();
	polycubeMesh.GetEdgePatches();
	polycubeMesh.GetVertexPatches();
	polycubeMesh.GetPatchesLabel();
	GetFaceType(polycubeMesh);
	polycubeMesh.GetPatchMinDistance();

	////////////////////////////////////////////////////////////////

	polycubeMesh.GetMaxMinCoordinates();
	polycubeMesh.AdjustPolycubeSurfaceToHoldAllSmallCubesVolume_quad(0.001, polycubeMesh.m_minVertex);
	polycubeMesh.AddCompensation_quad();

	MeshFileWriter polycubeFileWriter(polycubeMesh.V, polycubeMesh.C, strOutputResultHexFileName.c_str(), MESH_TYPE_HEXAHEDRON_VTK);
	polycubeFileWriter.WriteMeshFile(strOutputResultHexFileName.c_str(), MESH_TYPE_HEXAHEDRON_VTK);
	////////////////////////////////////////////////////////////////

	return 0;
}

