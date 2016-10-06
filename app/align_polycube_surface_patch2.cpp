/*
 * align_polycube_surface_patch2.cpp
 *
 *  Created on: Jun 1, 2015
 *      Author: cotrik
 */

#include "include.h"

int AlignPolycubeSurface(int argc, char* argv[]);

int main(int argc, char* argv[])
{
	return AlignPolycubeSurface(argc, argv);
}

int AlignPolycubeSurface(int argc, char* argv[])
{
	std::string strInputPolyTetFileName;
	std::string strInputPolyHexFileName;
	std::string strOutputPolyHexFileName;

	try
	{
#ifdef CGAL_USE_BOOST_PROGRAM_OPTIONS
		po::options_description desc("Allowed options");
		desc.add_options()
		("help,h", "produce help message.")
		// ########## input #########
		// -i
		("input-polycube-tet-file,i", po::value<std::string>(&strInputPolyTetFileName)->default_value("polycube_aligned.tet.vtk"),
			"name of the input polycube tetrahedral vtk file")
		// -p
		("input-polycube-hex-file,p", po::value<std::string>(&strInputPolyHexFileName)->default_value("polycube.hex.vtk"),
			"name of the input polycube hexahedral vtk file")
		// -o
		("output-polycube-hex-file,o", po::value<std::string>(&strOutputPolyHexFileName)->default_value("polycube_aligned.hex.vtk"),
			"name of the output polycube hexahedral vtk file")
		;

		po::variables_map vm;
		po::store(po::parse_command_line(argc, argv, desc), vm);
		po::notify(vm);

		if (vm.count("help"))		{			std::cout << desc << "\n";			return 1;		}
		if (vm.count("input-polycube-tet-file"))	std::cout << "strInputPolyTetFileName = " << strInputPolyTetFileName << std::endl;
		if (vm.count("input-polycube-hex-file"))	std::cout << "strInputPolyHexFileName = " << strInputPolyHexFileName << std::endl;
		if (vm.count("output-polycube-hex-file"))	std::cout << "strOutputPolyHexFileName = " << strOutputPolyHexFileName << std::endl;
#else
		strInputPolyTetFileName = "polycube.tet.vtk";
		strInputPolyHexFileName = "polycube.hex.vtk";
		strOutputPolyHexFileName = "aligned_polycube.hex.vtk";
#endif
	} catch (std::exception& e)
	{
		std::cerr << "error: " << e.what() << "\n";
		return 1;
	} catch (...)
	{
		std::cerr << "Exception of unknown type!\n";
	}

	//////////////////////////////////////////
	// Read Orig Tet Mesh
	//////////////////////////////////////////
	MeshFileReader polycubeTetMeshFileReader(strInputPolyTetFileName.c_str());
	PolycubeMesh polycubeTetMesh(polycubeTetMeshFileReader.GetMesh());
	polycubeTetMesh.NormalizeCoordinateValue();
	polycubeTetMesh.ExtractSurface();
	polycubeTetMesh.GetNormalOfSurfaceFaces();
	polycubeTetMesh.GetNormalOfSurfaceVertices();
	polycubeTetMesh.GetCorners();

	polycubeTetMesh.LabelSurfaceFace();
	polycubeTetMesh.GetFacePatches();
	polycubeTetMesh.GetEdgePatches();
	polycubeTetMesh.GetVertexPatches();
	polycubeTetMesh.GetPatchesLabel();

	polycubeTetMesh.GetFaceType();

	for (int i = 0; i < polycubeTetMesh.patchLabel.size(); i++)
	{
		std::cout << "patch " << i << " label: ";
		for (int j = 0; j < polycubeTetMesh.patchLabel.at(i).size(); j++)
		{
			std::cout << polycubeTetMesh.patchLabel.at(i).at(j) << " ";
		}
		std::cout << std::endl;
	}
	std::cout << "---------------------------------------------------------------- " << std::endl;
	MeshFileReader polycubeHexMeshFileReader(strInputPolyHexFileName.c_str());
	PolycubeMesh polycubeHexMesh(polycubeHexMeshFileReader.GetMesh());
	polycubeHexMesh.NormalizeCoordinateValue();
	polycubeHexMesh.ExtractSurface();
	polycubeHexMesh.GetNormalOfSurfaceFaces();
	polycubeHexMesh.GetNormalOfSurfaceVertices();
	polycubeHexMesh.GetCorners();

	polycubeHexMesh.LabelSurfaceFace();
	polycubeHexMesh.GetFacePatches();
	polycubeHexMesh.GetEdgePatches();
	polycubeHexMesh.GetVertexPatches();
	polycubeHexMesh.GetPatchesLabel();

	for (int i = 0; i < polycubeHexMesh.patchLabel.size(); i++)
	{
		std::cout << "patch " << i << " label: ";
		for (int j = 0; j < polycubeHexMesh.patchLabel.at(i).size(); j++)
		{
			std::cout << polycubeHexMesh.patchLabel.at(i).at(j) << " ";
		}
		std::cout << std::endl;
	}

	std::vector<double> patchPos;
	for (int i = 0; i < polycubeTetMesh.patchLabel.size(); i++)
	{
		std::vector<unsigned long>& facePatch = polycubeTetMesh.facePatches.at(i);
		std::vector<unsigned long>& vetexPatch = polycubeTetMesh.vertexPatches.at(i);
		double x = 0.0, y = 0.0, z = 0.0;
		const size_t vetexPatchSize = vetexPatch.size();
		for (int j = 0; j < vetexPatchSize; j++)
		{
			const Vertex& v = polycubeTetMesh.V.at(vetexPatch.at(j));
			x += v.x;
			y += v.y;
			z += v.z;
		}
		if (vetexPatchSize != 0)
		{
			x /= vetexPatchSize;
			y /= vetexPatchSize;
			z /= vetexPatchSize;
		    //----------------------------------------
			std::vector<unsigned long>& faceHexPatch = polycubeHexMesh.facePatches.at(i);
			std::vector<unsigned long>& vetexHexPatch = polycubeHexMesh.vertexPatches.at(i);
			const size_t vetexHexPatchSize = vetexHexPatch.size();
			const FACE_TYPE faceType = polycubeTetMesh.faceType.at(facePatch.at(0));

			if (faceType == FACE_X)
			{
				patchPos.push_back(x);
				double oldx = 0;
				for (int j = 0; j < vetexHexPatchSize; j++)
				{
					Vertex& v = polycubeHexMesh.V.at(vetexHexPatch.at(j));
					oldx = v.x;
					v.x = x;
				}
				for (int j = 0; j < polycubeHexMesh.V.size(); j++)
				{
					Vertex& v = polycubeHexMesh.V.at(j);
					if (fabs(v.x - oldx) < 5e-4)
						v.x = x;
				}
			}
			else if (faceType == FACE_Y)
			{
				patchPos.push_back(y);
				double oldy = 0;
				for (int j = 0; j < vetexHexPatchSize; j++)
				{
					Vertex& v = polycubeHexMesh.V.at(vetexHexPatch.at(j));
					oldy = v.y;
					v.y = y;
				}
				for (int j = 0; j < polycubeHexMesh.V.size(); j++)
				{
					Vertex& v = polycubeHexMesh.V.at(j);
					if (fabs(v.y - oldy) < 5e-4)
						v.y = y;
				}
			}
			else if (faceType == FACE_Z)
			{
				patchPos.push_back(z);
				double oldz = 0;
				for (int j = 0; j < vetexHexPatchSize; j++)
				{
					Vertex& v = polycubeHexMesh.V.at(vetexHexPatch.at(j));\
					oldz = v.z;
					v.z = z;
				}
				for (int j = 0; j < polycubeHexMesh.V.size(); j++)
				{
					Vertex& v = polycubeHexMesh.V.at(j);
					if (fabs(v.z - oldz) < 5e-4)
						v.z = z;
				}
			}
		}
	}
	double min_distance = 1;
	for (int i = 1; i < polycubeTetMesh.m_xLabel; i++)
	{
		double distance = patchPos.at(i) - patchPos.at(i - 1);
		std::cout << "patch distance = " << patchPos.at(i) << " - " << patchPos.at(i - 1) << " = " << distance << std::endl;
		if (distance < min_distance)
			min_distance = distance;
	}
	for (int i = polycubeTetMesh.m_xLabel + 1; i < polycubeTetMesh.m_yLabel; i++)
	{
		double distance = patchPos.at(i) - patchPos.at(i - 1);
		std::cout << "patch distance = " << patchPos.at(i) << " - " << patchPos.at(i - 1) << " = " << distance << std::endl;
		if (distance < min_distance)
			min_distance = distance;
	}
	for (int i = polycubeTetMesh.m_yLabel + 1; i < polycubeTetMesh.m_zLabel; i++)
	{
		double distance = patchPos.at(i) - patchPos.at(i - 1);
		std::cout << "patch distance = " << patchPos.at(i) << " - " << patchPos.at(i - 1) << " = " << distance << std::endl;
		if (distance < min_distance)
			min_distance = distance;
	}
	cout << "patch min_distance = " << min_distance << std::endl;

	MeshFileWriter polycubeFileWriter(polycubeHexMesh.V, polycubeHexMesh.C, strOutputPolyHexFileName.c_str());
	polycubeFileWriter.WriteMeshFile();
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//	pTetPolycubeMesh = (PolycubeMesh*)&polycubeTetMesh;
//	pHexPolycubeMesh = (PolycubeMesh*)&polycubeHexMesh;
//	DisplayMesh(argc, argv, *pHexPolycubeMesh);
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	return 0;
}
