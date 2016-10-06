#include "include.h"

int MeshDilate(int argc, char* argv[]);

int main(int argc, char* argv[])
{
	return MeshDilate(argc, argv);
}

int MeshDilate(int argc, char* argv[])
{
	std::string strInputFileName;
	std::string strOutputFileName;
	int iterNum = 3;
	try
	{
#ifdef CGAL_USE_BOOST_PROGRAM_OPTIONS
		po::options_description desc("Allowed options");
		desc.add_options()
			("help,h", "produce help message.")
			// -i
			("input-file,i", po::value<std::string>(&strInputFileName)->default_value("orig.tet.vtk"),
				"name of the orig input tetrahedral vtk file")
			// -o
			("output-file,o", po::value<std::string>(&strOutputFileName)->default_value("orig.density.tet.vtk"),
				"name of the output density vtk file")
			// -l
			("loop,l", po::value<int>(&iterNum)->default_value(1),
					"dilating iteration number")
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
		std::cout << "strOutputFileName = " << strOutputFileName << std::endl;
		std::cout << "dilating loop = " << iterNum << std::endl;
#else
		strInputFileName = "orig.tet.vtk";
		strOutputFileName = "orig.density.tet.vtk";
#endif
	} catch (std::exception& e)
	{
		std::cerr << "error: " << e.what() << "\n";
		return 1;
	} catch (...)
	{
		std::cerr << "Exception of unknown type!\n";
	}

	MeshFileReader tetMeshFileReader(strInputFileName.c_str(), MESH_TYPE_TETRAHEDRON_VTK);
	Mesh tetMesh(tetMeshFileReader.GetMesh());

	tetMesh.ExtractSurface();
	int count = 0;
	while(iterNum--)
	{
		count++;
	tetMesh.GetSurfaceDensityField();
	//tetMesh.SmoothSurface();
	tetMesh.GetNormalOfSurfaceVertices();
	// update vertex positions according to density filed and normale direction of surface
	for (int i = 0; i < tetMesh.V.size(); i++)
	{
		if (tetMesh.vertexInfo[i].bSurface && tetMesh.vec_densityFiled[i] > (0.4 - count*0.02))
		{
			Vertex& v = tetMesh.V.at(i);
//			v.x -= 0.2*tetMesh.vec_averageLength_V[i]*tetMesh.N_V[i].x;
//			v.y -= 0.2*tetMesh.vec_averageLength_V[i]*tetMesh.N_V[i].y;
//			v.z -= 0.2*tetMesh.vec_averageLength_V[i]*tetMesh.N_V[i].z;
			v.x -= 0.025*tetMesh.N_V[i].x;
			v.y -= 0.025*tetMesh.N_V[i].y;
			v.z -= 0.025*tetMesh.N_V[i].z;
			tetMesh.vertexInfo[i].bNeedSmoothing = true;
		}
//		else
//		{
//			tetMesh.vertexInfo[i].bNeedSmoothing = false;
//		}
	}
		if ((iterNum)%5 == 0)
		{
			tetMesh.SmoothSurface();
		}
	}
//	tetMesh.SmoothSurface();
	MeshFileWriter tetFileWriter(tetMesh, strOutputFileName.c_str(), MESH_TYPE_TETRAHEDRON_VTK);
	tetFileWriter.WriteMeshDesityFieldFile(strOutputFileName.c_str(), MESH_TYPE_TETRAHEDRON_VTK);

	return 0;
}
