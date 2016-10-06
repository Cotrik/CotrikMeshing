#include "include.h"

int HexDensity(int argc, char* argv[]);

int main(int argc, char* argv[])
{
	return HexDensity(argc, argv);
}

int HexDensity(int argc, char* argv[])
{
	std::string strInputOrigTetFileName;
	std::string strInputPolyTetFileName;
	std::string strOutputPolyHexFileName;
	std::string strOutputResultHexFileName;

	std::string strAjustedPolycubeTetMeshFilename;
	std::string strOrigNormalTetMeshFileName;
	std::string strPolycubeNormalTetMeshFileName;
	std::string strHexCpMeshFileName;
	std::string strMagnifiedHexCpMeshFileName;
	std::string strMagnifiedTetMeshFileName;

	std::string strInputFileName;
//	std::string strResultNormalTetMeshFileName;

	std::string strDir;
	DIRECTION dir = X_AXIS;
	float delta = 0.02f;
//	float cubesize = 0.08f;
	int iterNum = 3;
	int Func = 1;

	try
	{
#ifdef CGAL_USE_BOOST_PROGRAM_OPTIONS
		po::options_description desc("Allowed options");
		desc.add_options()
			("help,h", "produce help message.")
			("thread,j", po::value<unsigned int>(&THREAD_NUM)->default_value(4), "threads number")
			// ########## input #########
			// -i
			("input-hex-file,i", po::value<std::string>(&strInputFileName)->default_value("kitty.hex.vtk"),
				"name of the orig input hexahedral vtk file")
			// -o
			("out-hex-file,o", po::value<std::string>(&strOutputResultHexFileName)->default_value("kitty.smooth.hex.vtk"),
				"name of the output hexaheral vtk/mesh file")
			// -l
			("laplacian,l",	po::value<int>(&iterNum)->default_value(1),
				"laplacian smooth iteration number")
			;

		po::variables_map vm;
		po::store(po::parse_command_line(argc, argv, desc), vm);
		po::notify(vm);

		if (vm.count("help"))
		{
			std::cout << desc << "\n";
			return 1;
		}

		GetDir(strDir, dir);

		std::cout << "strInputFileName = " << strInputFileName << std::endl;
		std::cout << "strOutputResultHexFileName = " << strOutputResultHexFileName << std::endl;
		std::cout << "smooth iteration number = " << iterNum << std::endl;
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

	MeshFileReader hexMeshFileReader(strInputFileName.c_str(), MESH_TYPE_HEXAHEDRON_VTK);
	Mesh hexMesh(hexMeshFileReader.GetMesh());

	hexMesh.GetDensityField();

	MeshFileWriter hexFileWriter(hexMesh, strOutputResultHexFileName.c_str(), MESH_TYPE_HEXAHEDRON_VTK);
	hexFileWriter.WriteMeshDesityFieldFile(strOutputResultHexFileName.c_str(), MESH_TYPE_HEXAHEDRON_VTK);

	return 0;
}
