#include "include.h"

int AlignSurface(int argc, char* argv[]);

int main(int argc, char* argv[])
{
	return AlignSurface(argc, argv);
}

int AlignSurface(int argc, char* argv[])
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
			("input-polycube-file,i", po::value<std::string>(&strInputFileName)->default_value("polycube.tet.vtk"),
				"name of the polycube input tetrahedral vtk file")
			// -o
			("out-polycube-file,o", po::value<std::string>(&strOutputResultHexFileName)->default_value("polycube_aligned.tet.vtk"),
				"name of the output polycube tetrahedral vtk/mesh file")
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
		strInputFileName = "polycube.tet.vtk";
		strOutputResultHexFileName = "polycube_aligned.tet.vtk";
#endif
	} catch (std::exception& e)
	{
		std::cerr << "error: " << e.what() << "\n";
		return 1;
	} catch (...)
	{
		std::cerr << "Exception of unknown type!\n";
	}

	MeshFileReader polycubeMeshFileReader(strInputFileName.c_str());
	Mesh polycubeMesh(polycubeMeshFileReader.GetMesh());

	polycubeMesh.ExtractSurface();
	for (int i = 0; i < polycubeMesh.surface.size(); i++)
	{
		const Cell& tri = polycubeMesh.surface.at(i);
		for (int j = 0; j < 3; j++)
		{
			Vertex& v = polycubeMesh.V.at(tri.at(j));
			if (v.vinfo.bUsed)
				continue;
			v.x = (round(v.x/discrete))*discrete;
			v.y = (round(v.y/discrete))*discrete;
			v.z = (round(v.z/discrete))*discrete;
			v.vinfo.bUsed = true;
		}
	}

	MeshFileWriter polycubeFileWriter(polycubeMesh.V, polycubeMesh.C, strOutputResultHexFileName.c_str());
	polycubeFileWriter.WriteMeshFile();

	return 0;
}
