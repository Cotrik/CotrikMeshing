#include "include.h"

int HexifyMesh_plus(int argc, char* argv[]);

int main(int argc, char* argv[])
{
	return HexifyMesh_plus(argc, argv);
}

int HexifyMesh_plus(int argc, char* argv[])
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
			("orig-input-tet-file,i", po::value<std::string>(&strInputOrigTetFileName)->default_value("orig.tet.vtk"),
				"name of the orig input tetrahedral vtk file")
			// -p
			("polycube-input-tet-file,p", po::value<std::string>(&strInputPolyTetFileName)->default_value("polycube.tet.vtk"),
				"name of the polycube input tetrahedral vtk file")
			// ########## output #########
			// -I
			("orig-normalize-tet-file,I", po::value<std::string>(&strOrigNormalTetMeshFileName)->default_value("orig.normal.tet.vtk"),
				"name of the orig-normalize-tet-file vtk/mesh file")
			// -P
			("polycube-normalize-tet-file,P", po::value<std::string>(&strPolycubeNormalTetMeshFileName)->default_value("orig.normal.tet.vtk"),
				"name of the orig-normal-tet-file vtk/mesh file")
			// -n
			("out-polycube-hex-file,n",	po::value<std::string>(&strOutputPolyHexFileName)->default_value("poly.hex.vtk"),
				"name of the output polycube hexaheral vtk/mesh file")
			// -o
			("out-hex-file,o", po::value<std::string>(&strOutputResultHexFileName)->default_value("result.hex.vtk"),
				"name of the output hexaheral vtk/mesh file")
            // -C
			("cube-cp-file,C", po::value<std::string>(&strHexCpMeshFileName)->default_value("Cp.hex.vtk"),
				"name of the cube Cp vtk/mesh file")
			// -M
			("magnified-hex-Cp-file,M", po::value<std::string>(&strMagnifiedHexCpMeshFileName)->default_value("Cp.Magnified.hex.vtk"),
				"name of the magnified hex Cp vtk/mesh file")
            // -m
			("magnified-tet-file,m", po::value<std::string>(&strMagnifiedTetMeshFileName)->default_value("magnified.tet.vtk"),
				"name of the magnified tet vtk/mesh file")
			// -A
			("Ajust-hex-file,A", po::value<std::string>(&strAjustedPolycubeTetMeshFilename)->default_value("ajusted.polycube.tet.vtk"),
				"name of the ajust tet polycube vtk/mesh file")
			// ########## parameters #########
			// -s
			("hex-parameterization-size,s", po::value<float>(&delta)->default_value(0.02),
				"cube size, 0.001f <= size <= 0.1f")
			// -S
			("Index-Scale-Parameters,S", po::value<std::string>(&strScaleParas)->default_value("scale_paras.dat"),
				"Index-Scale-Parameters file name")
			// -c
			("magnifier-cube-size,c", po::value<float>(&CUBESIZE)->default_value(0.08),
				"cube size, 0.05f <= size <= 0.1f")
			// -w
			("offset,w", po::value<float>(&__plus)->default_value(0.001),
				"0.01s <= offset <= s")
			// -F
			("Function,F", po::value<int>(&Func)->default_value(1),
				"0 Hexify, 1 Magnify")
			// -r
			("ray-direction,r",	po::value<std::string>(&strDir)->default_value("Z"),
				"ray-direction, X, Y, or Z")
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

		if (vm.count("orig-input-tet-file"))
		std::cout << "strInputOrigTetFileName = " << strInputOrigTetFileName << std::endl;
		if (vm.count("polycube-input-tet-file"))
		std::cout << "strInputPolyTetFileName = " << strInputPolyTetFileName << std::endl;
		if (vm.count("out-polycube-hex-file"))
		std::cout << "strOutputPolyHexFileName = " << strOutputPolyHexFileName << std::endl;
		if (vm.count("out-hex-file"))
		std::cout << "strOutputResultHexFileName = " << strOutputResultHexFileName << std::endl;
		if (vm.count("Ajust-hex-file"))
		std::cout << "strAjustedPolycubeTetMeshFilename = " << strAjustedPolycubeTetMeshFilename << std::endl;
		if (vm.count("orig-normalize-tet-file"))
		std::cout << "strOrigNormalTetMeshFileName = " << strOrigNormalTetMeshFileName << std::endl;
		if (vm.count("polycube-normalize-tet-file"))
		std::cout << "strPolycubeNormalTetMeshFileName = " << strPolycubeNormalTetMeshFileName << std::endl;
		if (vm.count("cube-cp-file"))
		std::cout << "strHexCpMeshFileName = " << strHexCpMeshFileName << std::endl;
		if (vm.count("magnified-hex-Cp-file"))
		std::cout << "strMagnifiedHexCpMeshFileName = " << strMagnifiedHexCpMeshFileName << std::endl;
		if (vm.count("magnified-tet-file"))
		std::cout << "strMagnifiedTetMeshFileName = " << strMagnifiedTetMeshFileName << std::endl;
		if (vm.count("Index-Scale-Parameter"))
		std::cout << "Index-Scale-Parameter File = " << strScaleParas << std::endl;
		if (vm.count("thread"))
		std::cout << "threadnum = " << THREAD_NUM << std::endl;
		if (vm.count("hex-parameterization-size"))
		std::cout << "hexsize = " << delta << std::endl;
		if (vm.count("magnifier-cube-size"))
		std::cout << "cubesize = " << CUBESIZE << std::endl;
		if (vm.count("ray-direction"))
		std::cout << "ray direction = " << strDir << std::endl;
		if (vm.count("laplacian"))
		std::cout << "smooth iteration number = " << iterNum << std::endl;
		if (vm.count("offset"))
		std::cout << "offset = " << __plus << std::endl;

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

	//////////////////////////////////////////
	// Read Orig Tet Mesh
	//////////////////////////////////////////
	MeshFileReader orgTetMeshFileReader(strInputOrigTetFileName.c_str(), MESH_TYPE_TETRAHEDRON_VTK);
	Mesh origTetMesh(orgTetMeshFileReader.GetMesh());

	//////////////////////////////////////////
	// Get Orig Tet Surface Curvature
	//////////////////////////////////////////
	origTetMesh.ExtractSurface();
	origTetMesh.NormalizeCoordinateValue();
	MeshFileWriter normalTetMeshFileWriter(origTetMesh);
	normalTetMeshFileWriter.WriteMeshFile(strOrigNormalTetMeshFileName.c_str(), MESH_TYPE_TETRAHEDRON_VTK);
	std::string strTetSurfaceFileName = strInputOrigTetFileName + std::string(".tet.surface.off");
	MeshFileWriter tetSurfaceMeshFileWriter(origTetMesh.V, origTetMesh.surface, strTetSurfaceFileName.c_str(), MESH_TYPE_TRIANGLE_OFF);
	tetSurfaceMeshFileWriter.WriteMeshFile(strTetSurfaceFileName.c_str(), MESH_TYPE_TRIANGLE_OFF);

	if (argc == 3)
	{
		return 0;
	}

	std::string strTetSurfaceCurvatureFileName = strInputOrigTetFileName + std::string(".tet.surface.off.cur.txt");
	MeshFileWriter orgTetMeshFileWriter(origTetMesh);
	orgTetMeshFileWriter.WriteTetrahedralSurfaceCurvature(strTetSurfaceCurvatureFileName.c_str(), strInputOrigTetFileName.c_str());
	//////////////////////////////////////////
	// Read Polycube Tet Mesh
	//////////////////////////////////////////
	MeshFileReader polycubeTetMeshFileReader(strInputPolyTetFileName.c_str(), MESH_TYPE_TETRAHEDRON_VTK);
	PolycubeMesh polycubeTetMesh(polycubeTetMeshFileReader.GetMesh());
	////////////////
	char filename1[64] = { 0 };
	sprintf(filename1, "%siter%d.vtk", strInputPolyTetFileName.c_str(), 1);
	MeshFileWriter meshFileWriter1(polycubeTetMesh.V, polycubeTetMesh.C, filename1, MESH_TYPE_TETRAHEDRON_VTK);
	meshFileWriter1.WriteMeshFile();
	////////////
	polycubeTetMesh.ExtractSurface();
	polycubeTetMesh.GetCorners();
	polycubeTetMesh.NormalizeCoordinateValue();
	MeshFileWriter polycubeNormalTetMeshFileWriter(polycubeTetMesh);
	polycubeNormalTetMeshFileWriter.WriteMeshFile(strPolycubeNormalTetMeshFileName.c_str(), MESH_TYPE_TETRAHEDRON_VTK);

	std::string strPolycubeSurfaceFileName = strInputPolyTetFileName + std::string(".tet.surface.off");
	MeshFileWriter polycubeSurfaceMeshFileWriter(polycubeTetMesh.V, polycubeTetMesh.surface, strPolycubeSurfaceFileName.c_str(), MESH_TYPE_TRIANGLE_OFF);
	polycubeSurfaceMeshFileWriter.WriteMeshFile(strPolycubeSurfaceFileName.c_str(), MESH_TYPE_TRIANGLE_OFF);
	std::string strPolyTetSurfaceCurvatureFileName = strInputOrigTetFileName + std::string(".tet.surface.off.cur.txt");
	MeshFileWriter polycubeTetMeshFileWriter(polycubeTetMesh);
	polycubeTetMeshFileWriter.WriteTetrahedralSurfaceCurvature(strPolyTetSurfaceCurvatureFileName.c_str(), strInputPolyTetFileName.c_str());

	// smoothing
	polycubeTetMesh.GetVertexInfo();
	polycubeTetMesh.GetCurvature(strPolyTetSurfaceCurvatureFileName.c_str());
	polycubeTetMesh.SmoothCurvature(10);
	origTetMesh.vertexInfo = polycubeTetMesh.vertexInfo;
	polycubeTetMesh.GetMaxMinCoordinates();
	polycubeTetMesh.AdjustPolycubeSurfaceToHoldAllSmallCubesVolume(delta, polycubeTetMesh.m_minVertex);

	for (int i = 0; i < iterNum; i++)
	{
		int count = i + 1;
		polycubeTetMesh.CurvatureWeightingSmooth(count);
		char filename[64] = { 0 };
		sprintf(filename, "%siter%d.tet.vtk", strInputPolyTetFileName.c_str(), count);
		MeshFileWriter meshFileWriter(polycubeTetMesh.V, polycubeTetMesh.C, filename, MESH_TYPE_TETRAHEDRON_VTK);
		meshFileWriter.SetVertexInfo(polycubeTetMesh.vertexInfo);
		meshFileWriter.WriteMeshFile();
	}

	polycubeTetMesh.AddCompensation();
	MeshFileWriter ajustedmeshFileWriter(polycubeTetMesh, strAjustedPolycubeTetMeshFilename.c_str(), MESH_TYPE_TETRAHEDRON_VTK);
	ajustedmeshFileWriter.WriteMeshFile();
	strAjustedPolycubeTetMeshFilename = strInputPolyTetFileName + std::string(".ajusted.vtk");
	polycubeTetMesh.GenerateHexahedralMesh_plus(origTetMesh, strOutputPolyHexFileName.c_str(), strOutputResultHexFileName.c_str(), delta, dir);

//	MeshFileReader hexMeshFileReader(strOutputResultHexFileName.c_str(), MESH_TYPE_HEXAHEDRON_VTK);
//	Mesh hexMesh(hexMeshFileReader.GetMesh());
//	hexMesh.ExtractSurface();
//	std::string strHexSurfaceFileName = strOutputResultHexFileName + std::string(".hex.surface.off");
//	MeshFileWriter hexSurfaceMeshFileWriter(hexMesh.V, hexMesh.surface, strHexSurfaceFileName.c_str(),	MESH_TYPE_QUAD_OFF);
//	hexSurfaceMeshFileWriter.WriteMeshFile(strHexSurfaceFileName.c_str(), MESH_TYPE_QUAD_OFF);

	return 0;
}

