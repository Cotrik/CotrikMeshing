/*
 * align_polycube_surface_patch.cpp
 *
 *  Created on: May 30, 2015
 *      Author: cotrik
 */

#include "include.h"

int AlignSurface_Patch(int argc, char* argv[]);

int main(int argc, char* argv[])
{
	return AlignSurface_Patch(argc, argv);
}

int AlignSurface_Patch(int argc, char* argv[])
{
	std::string strInputPolyTetFileName;
	std::string strOutputPolyTetFileName;
	float discrete = 0.001f;
	try
	{
#ifdef CGAL_USE_BOOST_PROGRAM_OPTIONS
		po::options_description desc("Allowed options");
		desc.add_options()
		("help,h", "produce help message.")
		// ########## input #########
		// -i
		("input-tet-polycube,i", po::value<std::string>(&strInputPolyTetFileName)->default_value("polycube.tet.vtk"),
			"name of the input polycube tetrahedral vtk file")
		// -o
		("output-tet-polycube,o", po::value<std::string>(&strOutputPolyTetFileName)->default_value("polycube_aligned.tet.vtk"),
			"name of the input output tetrahedral vtk file")
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

		if (vm.count("input-tet-polycube"))	std::cout << "strInputPolyTetFileName = " << strInputPolyTetFileName << std::endl;
		if (vm.count("output-tet-polycube"))std::cout << "strOutputPolyTetFileName = " << strOutputPolyTetFileName << std::endl;
#else
		strInputPolyTetFileName = "polycube.tet.vtk";
		strOutputPolyTetFileName = "polycube_aligned.tet.vtk";
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
	// Read Polycube Tet Mesh
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
	//////////////////////////////////////////////
	polycubeTetMesh.GetVertexInfo();
	// Polycube Mesh cleaning
	bool hasLabeledAllSurfaceFace = false;
	std::vector<Vertex>& PV = polycubeTetMesh.V;
	const std::vector<Face>& PF = polycubeTetMesh.surface;
	std::vector<FACE_TYPE>& PFT = polycubeTetMesh.faceType;
	std::vector<int>& PFL = polycubeTetMesh.faceLabel;
	const size_t PF_SIZE = PF.size();
	polycubeTetMesh.GetFaceType();

//	for (int i = 0; i < PF_SIZE; i++)
//	{
//		const Face& f = PF.at(i);
//		//PFT.push_back(GetFaceType(polycubeTetMesh.N_F.at(i)));
//		const Vertex* pVertex0 = &polycubeTetMesh.V[f.at(0)];
//		const Vertex* pVertex1 = &polycubeTetMesh.V[f.at(1)];
//		const Vertex* pVertex2 = &polycubeTetMesh.V[f.at(2)];
//
//		const glm::vec3 v0(pVertex0->x, pVertex0->y, pVertex0->z);
//		const glm::vec3 v1(pVertex1->x, pVertex1->y, pVertex1->z);
//		const glm::vec3 v2(pVertex2->x, pVertex2->y, pVertex2->z);
//
//		const glm::vec3 v10 = v0 - v1;
//		const glm::vec3 v12 = v2 - v1;
//		const glm::vec3 normal = glm::cross(v12, v10);
//
//		PFT.push_back(GetFaceType(normal));
//	}

	for (int i = 0; i < polycubeTetMesh.patchLabel.size(); i++)
	{
		std::cout << "patch " << i << " label: ";
		for (int j = 0; j < polycubeTetMesh.patchLabel.at(i).size(); j++)
		{
			std::cout << polycubeTetMesh.patchLabel.at(i).at(j) << " ";
		}
		std::cout << std::endl;
	}

	for (int i = 0; i < polycubeTetMesh.facePatches.size(); i++)
	{
		std::vector<unsigned long>& facePatch = polycubeTetMesh.facePatches.at(i);
		std::vector<unsigned long>& vetexPatch = polycubeTetMesh.vertexPatches.at(i);
		double x = 0.0, y = 0.0, z = 0.0;
		const size_t vetexPatchSize = vetexPatch.size();
		for (int j = 0; j < vetexPatchSize; j++)
		{
			const Vertex& v = PV.at(vetexPatch.at(j));
			x += v.x; y += v.y; z += v.z;
		}
		if (vetexPatchSize != 0)
		{
			x /= vetexPatchSize; y /= vetexPatchSize; z /= vetexPatchSize;
//			double rx = (round(x/discrete))*discrete;
//			double ry = (round(y/discrete))*discrete;
//			double rz = (round(z/discrete))*discrete;
//
//			const double minx = x < rx ? x : rx;
//			const double miny = y < ry ? y : ry;
//			const double minz = z < rz ? z : rz;
//
//			const double maxx = x > rx ? x : rx;
//			const double maxy = y > ry ? y : ry;
//			const double maxz = z > rz ? z : rz;

			const FACE_TYPE faceType = PFT.at(facePatch.at(0));
			if (faceType == FACE_X)
			{
				for (int j = 0; j < vetexPatchSize; j++)
				{
					Vertex& v = PV.at(vetexPatch.at(j));
					v.x = x;
					//v.x = rx;
				}
			}
			else if (faceType == FACE_Y)
			{
				for (int j = 0; j < vetexPatchSize; j++)
				{
					Vertex& v = PV.at(vetexPatch.at(j));
					v.y = y;
					//v.y = ry;
				}
			}
			else if (faceType == FACE_Z)
			{
				for (int j = 0; j < vetexPatchSize; j++)
				{
					Vertex& v = PV.at(vetexPatch.at(j));
					v.z = z;
					//v.z = rz;
				}
			}
		}
	}

	MeshFileWriter polycubeFileWriter(polycubeTetMesh.V, polycubeTetMesh.C, strOutputPolyTetFileName.c_str());
	polycubeFileWriter.WriteMeshFile();
	return 0;
}

