/*
 * detect_overlapping_area.cpp
 *
 *  Created on: Jul 2, 2015
 *      Author: cotrik
 */

#include "include.h"

int main(int argc, char* argv[])
{
	std::string strInputOrigTetFileName;
	std::string strInputPolyTetFileName;
	std::string strInputPolyHexFileName;
	std::string strOutputHexFileName;

	try
	{
#ifdef CGAL_USE_BOOST_PROGRAM_OPTIONS
		po::options_description desc("Allowed options");
		desc.add_options()
			("help,h", "produce help message.")
			("thread,j", po::value<unsigned int>(&THREAD_NUM)->default_value(4), "threads number")
			// ########## input #########
			// -i
			("input-tet-file,i", po::value<std::string>(&strInputOrigTetFileName)->default_value("orig.tet.vtk"),
				"name of the orig input tetrahedral vtk file")
			// -t
			("polycube-tet-file,t", po::value<std::string>(&strInputPolyTetFileName)->default_value("polycube.tet.vtk"),
				"name of the polycube input tetrahedral vtk file")
			// -p
			("polycube-hex-file,p", po::value<std::string>(&strInputPolyHexFileName)->default_value("polycube.hex.vtk"),
				"name of the polycube input hexahedral vtk file")
			// -o
			("output-hex-file,o", po::value<std::string>(&strOutputHexFileName)->default_value("hex.vtk"),
				"name of the output hexahedral vtk file")
			// ########## output #########
			;

		po::variables_map vm;
		po::store(po::parse_command_line(argc, argv, desc), vm);
		po::notify(vm);

		if (vm.count("help"))
		{
			std::cout << desc << "\n";
			return 1;
		}

		if (vm.count("input-tet-file"))		std::cout << "strInputOrigTetFileName = " << strInputOrigTetFileName << std::endl;
		if (vm.count("polycube-tet-file"))		std::cout << "strInputPolyTetFileName = " << strInputPolyTetFileName << std::endl;
		if (vm.count("polycube-hex-file"))		std::cout << "strInputPolyHexFileName = " << strInputPolyHexFileName << std::endl;
		if (vm.count("output-hex-file"))		std::cout << "strOutputHexFileName = " << strOutputHexFileName << std::endl;
#else
		strInputOrigTetFileName = "orig.tet.vtk";
		strInputPolyTetFileName = "polycube.tet.vtk";
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

	MeshFileReader orgTetMeshFileReader(strInputOrigTetFileName.c_str(), MESH_TYPE_TETRAHEDRON_VTK);
	Mesh origTetMesh(orgTetMeshFileReader.GetMesh());

	MeshFileReader polycubeTetMeshFileReader(strInputPolyTetFileName.c_str(), MESH_TYPE_TETRAHEDRON_VTK);
	PolycubeMesh polycubeTetMesh(polycubeTetMeshFileReader.GetMesh());

	MeshFileReader polycubeHexMeshFileReader(strInputPolyHexFileName.c_str(), MESH_TYPE_HEXAHEDRON_VTK);
	PolycubeMesh polycubeHexMesh(polycubeHexMeshFileReader.GetMesh());

	origTetMesh.ExtractSurface();
	polycubeTetMesh.ExtractSurface();
	polycubeHexMesh.ExtractSurface();

	polycubeTetMesh.GetSurfaceVertexIndices();
	MeshFileWriter origTetMeshSurfaceFileWriter(origTetMesh.V, origTetMesh.surface, "orig.tet.surface.vtk", MESH_TYPE_TRIANGLE_VTK);
	origTetMeshSurfaceFileWriter.WriteMeshFile("orig.tet.surface.vtk", MESH_TYPE_TRIANGLE_VTK);
	origTetMesh.InitGeodesicDistance("orig.tet.surface.vtk");

	TetHexParas p;
	std::vector<OverlappingParas> overlappingParas;
	std::vector<unsigned long> overlapping_hex_indices;
	GetParameterAndDetectIntersection(origTetMesh, polycubeTetMesh, polycubeHexMesh, p, overlappingParas, overlapping_hex_indices);

	std::vector<Cell> hexC;
	for (int i = 0; i < polycubeHexMesh.C.size(); i++)
	{
		const Cell& hexCell = polycubeHexMesh.C.at(i);
		bool overlappingCell = false;
		for (int j = 0; j < overlapping_hex_indices.size(); j++)
		{
			const unsigned long hexIndex = overlapping_hex_indices.at(j);
			if (i == hexIndex)
			{
				overlappingCell = true;
				break;
			}
		}
		if (!overlappingCell)
		{
			hexC.push_back(hexCell);
		}
	}

	std::vector<Vertex> hexV;
	MapbackToOrigTet(origTetMesh, polycubeHexMesh, p, hexV);

	MeshFileWriter origHexMeshFileWriter(hexV, hexC/*polycubeHexMesh.C*/, strOutputHexFileName.c_str(), MESH_TYPE_HEXAHEDRON_VTK);
	origHexMeshFileWriter.WriteMeshFile(strOutputHexFileName.c_str(), MESH_TYPE_HEXAHEDRON_VTK);

	Mesh hexMesh(hexV, polycubeHexMesh.C, MESH_TYPE_HEXAHEDRON_VTK);

	pHexMesh         = (Mesh*)&hexMesh;

	////////////////////////////////////////////////////////////////////////////////
	polycubeHexMesh.NormalizeCoordinateValue();
	polycubeHexMesh.ExtractSurface();
	polycubeHexMesh.GetNormalOfSurfaceFaces();
	polycubeHexMesh.GetNormalOfSurfaceVertices();
	polycubeHexMesh.GetCorners();
	polycubeHexMesh.GetVertexInfo();

	polycubeHexMesh.LabelSurfaceFace();
	polycubeHexMesh.GetFacePatches();
	polycubeHexMesh.GetEdgePatches();
	polycubeHexMesh.GetVertexPatches();
	polycubeHexMesh.GetPatchesLabel();

	for (int i = 0; i < polycubeHexMesh.surface.size(); i++)
	{
		const Face& f = polycubeHexMesh.surface.at(i);
		//PFT.push_back(GetFaceType(mesh.N_F.at(i)));
		const Vertex* pVertex0 = &polycubeHexMesh.V[f.at(0)];
		const Vertex* pVertex1 = &polycubeHexMesh.V[f.at(1)];
		const Vertex* pVertex2 = &polycubeHexMesh.V[f.at(2)];

		const glm::vec3 v0(pVertex0->x, pVertex0->y, pVertex0->z);
		const glm::vec3 v1(pVertex1->x, pVertex1->y, pVertex1->z);
		const glm::vec3 v2(pVertex2->x, pVertex2->y, pVertex2->z);

		const glm::vec3 v10 = v0 - v1;
		const glm::vec3 v12 = v2 - v1;
		const glm::vec3 normal = glm::cross(v12, v10);

		polycubeHexMesh.faceType.push_back(GetFaceType(normal));

		double area = 0.5f * glm::length(normal);
		polycubeHexMesh.faceArea.push_back(area);
	}
	polycubeHexMesh.CheckPatchesConnection();


	pTetPolycubeMesh = (PolycubeMesh*)&polycubeTetMesh;
	pHexPolycubeMesh = (PolycubeMesh*)&polycubeHexMesh;
	DisplayMesh(argc, argv, polycubeHexMesh);
	return 0;
}

