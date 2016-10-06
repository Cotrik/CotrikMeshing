#include "include.h"

int TetDensity(int argc, char* argv[]);

int main(int argc, char* argv[])
{
	return TetDensity(argc, argv);
}


int TetDensity(int argc, char* argv[])
{
	std::string strInputFileName;
	std::string strOutputFileName;
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

	tetMesh.GetDensityField();

	MeshFileWriter tetFileWriter(tetMesh, strOutputFileName.c_str(), MESH_TYPE_TETRAHEDRON_VTK);
	tetFileWriter.WriteMeshDesityFieldFile(strOutputFileName.c_str(), MESH_TYPE_TETRAHEDRON_VTK);

	return 0;
}

