#include "include.h"

int HexScaledJacobian(int argc, char* argv[]);

int main(int argc, char* argv[])
{
	return HexScaledJacobian(argc, argv);
}

int HexScaledJacobian(int argc, char* argv[])
{
	std::string strOutputScaledJacobianFileName;
	std::string strInputFileName;
	bool bShowMesh = false;
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
			("out-csv-file,o", po::value<std::string>(&strOutputScaledJacobianFileName)->default_value("scale_jacobian.csv"),
				"output statistics file")
			// -s
			("showmesh,s", po::value<bool>(&bShowMesh)->default_value(false),
				"show mesh file")
			;
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
		std::cout << "strOutputScaledJacobianFileName = " << strOutputScaledJacobianFileName << std::endl;
#else
		strInputFileName = "hex.vtk";
		strOutputScaledJacobianFileName = "scale_jacobian.csv";
#endif
	} catch (std::exception& e)
	{
		std::cerr << "error: " << e.what() << "\n";
		return 1;
	} catch (...)
	{
		std::cerr << "Exception of unknown type!\n";
	}

	MeshFileReader hexMeshFileReader(strInputFileName.c_str());
	Mesh hexMesh(hexMeshFileReader.GetMesh());
	std::cout << "#Points: "<< hexMesh.V.size() << " #Cells: " << hexMesh.C.size() << std::endl;
//	MeshFileReader cubeMeshFileReader("polycube.hex.vtk", MESH_TYPE_HEXAHEDRON_VTK);
//	Mesh cubeMesh(cubeMeshFileReader.GetMesh());

	hexMesh.OutputScaledJacobianDataFile(hexMesh.C, strOutputScaledJacobianFileName.c_str());
	hexMesh.ExtractSurface();
//	cubeMesh.ExtractSurface();
//
//	//cubeMesh.vertexInfo = hexMesh.vertexInfo;
//	cubeMesh.invertedCellIndex = hexMesh.invertedCellIndex;
//	cubeMesh.lowQualityCellIndex = hexMesh.lowQualityCellIndex;

//	if (bShowMesh)
//	{
//		pHexMesh = (Mesh*) &hexMesh;
//		pHexPolycubeMesh = (PolycubeMesh*) &hexMesh;
//		DisplayMesh(argc, argv, hexMesh);
//	}
	return 0;
}
