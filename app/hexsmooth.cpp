#include "include.h"

void GetDir(const std::string& strDir, DIRECTION& dir)
{
    if (strDir.empty())     return;
    if (strDir[0] == 'X')       dir = X_AXIS;
    else if (strDir[0] == 'Y')  dir = Y_AXIS;
    else if (strDir[0] == 'Z')  dir = Z_AXIS;
}

int HexSmooth(int argc, char* argv[]);

int main(int argc, char* argv[])
{
	return HexSmooth(argc, argv);
}

int HexSmooth(int argc, char* argv[])
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
	bool isPolycube = false;
	bool bSmoothBoundaryLine = false;
	bool bSmoothSurface = true;
	bool bSmoothVolume = false;

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
			// -c
			("isPolycube,c",	po::value<bool>(&isPolycube)->default_value(false),
				"Is a polycube mesh")
			// -v
			("volume,v",	po::value<bool>(&bSmoothVolume)->default_value(false),
				"whether smooth volume")
			// -s
			("surface,s",	po::value<bool>(&bSmoothSurface)->default_value(true),
				"whether smooth surface")
			// -b
			("boundaryline,b",	po::value<bool>(&bSmoothBoundaryLine)->default_value(false),
				"whether smooth boundaryline")
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

	MeshFileReader hexMeshFileReader(strInputFileName.c_str());
	if (isPolycube)
	{
        PolycubeMesh hexMesh(hexMeshFileReader.GetMesh());

        hexMesh.ExtractSurface();
        hexMesh.NormalizeCoordinateValue();
        hexMesh.ExtractSurface();
        hexMesh.GetNormalOfSurfaceFaces();
        hexMesh.GetNormalOfSurfaceVertices();
        hexMesh.GetCorners();
        hexMesh.GetVertexInfo();

        hexMesh.LabelSurfaceFace();
        hexMesh.GetFacePatches();
        hexMesh.GetEdgePatches();
        hexMesh.GetVertexPatches();
        hexMesh.GetPatchesLabel();
        hexMesh.GetFaceType();
        if (isPolycube)
        hexMesh.GetCorners();
        if (bSmoothBoundaryLine)
        {
            hexMesh.LaplacianWeightingSmoothBoundaryLine(iterNum);
        }
        if (bSmoothSurface)
        {
            hexMesh.SmoothSurface(iterNum);
        }
        if (bSmoothVolume)
            hexMesh.SmoothVolume(iterNum);
        hexMesh.GetMaxMinCoordinates();
        hexMesh.AdjustPolycubeSurfaceToHoldAllSmallCubesVolume_quad(0.02, hexMesh.m_minVertex);

        std::string strHexSurfaceFileName = strOutputResultHexFileName + std::string(".hex.surface.off");
        MeshFileWriter hexSurfaceMeshFileWriter(hexMesh.V, hexMesh.surface, strHexSurfaceFileName.c_str());
        hexSurfaceMeshFileWriter.WriteMeshFile();

        MeshFileWriter hexFileWriter(hexMesh.V, hexMesh.C, strOutputResultHexFileName.c_str());
        hexFileWriter.WriteMeshFile();
	}
	else
	{
		Mesh hexMesh(hexMeshFileReader.GetMesh());

		hexMesh.ExtractSurface();
//		if (isPolycube)
//		hexMesh.GetCorners();
		if (bSmoothSurface)
		hexMesh.SmoothSurface(iterNum);
		if (bSmoothVolume)
		hexMesh.SmoothVolume(iterNum);
		std::string strHexSurfaceFileName = strOutputResultHexFileName + std::string(".hex.surface.off");
		MeshFileWriter hexSurfaceMeshFileWriter(hexMesh.V, hexMesh.surface, strHexSurfaceFileName.c_str());
		hexSurfaceMeshFileWriter.WriteMeshFile();

		MeshFileWriter hexFileWriter(hexMesh.V, hexMesh.C, strOutputResultHexFileName.c_str());
		hexFileWriter.WriteMeshFile();
	}

	return 0;
}
