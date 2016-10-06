/*
 * polycubehexgen.cpp
 *
 *  Created on: Nov 3, 2015
 *      Author: cotrik
 */

#include "include.h"

void GetDir(const std::string& strDir, DIRECTION& dir)
{
    if (strDir.empty())     return;
    if (strDir[0] == 'X')       dir = X_AXIS;
    else if (strDir[0] == 'Y')  dir = Y_AXIS;
    else if (strDir[0] == 'Z')  dir = Z_AXIS;
}
static int ParseArg(int argc, char* argv[],
        std::string& strInputPolyTriFileName,
        std::string& strOutputPolyHexFileName,
        std::string& strAjustedPolycubeTriMeshFilename,
        std::string& strPolycubeNormalTriMeshFileName,
        std::string& strDir,
        DIRECTION& dir, float& delta, int& iterNum,
        float& x, float& y, float& z
        );

int main(int argc, char* argv[])
{
    std::string strInputPolyTriFileName;
    std::string strOutputPolyHexFileName;

    std::string strAjustedPolycubeTriMeshFilename;
    std::string strOrigNormalTriMeshFileName;
    std::string strPolycubeNormalTriMeshFileName;
    std::string strDir;
    DIRECTION dir = Y_AXIS;
    float delta = 0.02f;
    int iterNum = 0;
//    int Func = 1;
    float x = 1.0f, y = 1.0f, z = 1.0f;
    ParseArg(argc, argv, strInputPolyTriFileName,
            strOutputPolyHexFileName,
            strAjustedPolycubeTriMeshFilename,
            strPolycubeNormalTriMeshFileName,
            strDir, dir, delta, iterNum, x, y, z);

    //////////////////////////////////////////
    // Read Polycube Tet Mesh
    //////////////////////////////////////////
    std::cout << "------------- Processing Polycube Triangle Mesh -------------- " << std::endl;
    MeshFileReader polycubeTriMeshFileReader(strInputPolyTriFileName.c_str());
    PolycubeMesh polycubeTriMesh(polycubeTriMeshFileReader.GetMesh());
    polycubeTriMesh.ExtractSurface();
    polycubeTriMesh.NormalizeCoordinateValue();
    polycubeTriMesh.GetCorners();
    MeshFileWriter polycubeNormalTetMeshFileWriter(polycubeTriMesh, strPolycubeNormalTriMeshFileName.c_str());
    polycubeNormalTetMeshFileWriter.WriteFile();
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    polycubeTriMesh.GetVertexInfo();
    polycubeTriMesh.GetMaxMinCoordinates();
    polycubeTriMesh.GetFaceType();
    polycubeTriMesh.GetNormalOfSurfaceFaces();
    polycubeTriMesh.GetNormalOfSurfaceVertices();
    polycubeTriMesh.GetCorners();
    polycubeTriMesh.GetVertexInfo();
    ////////////////////////////////////////////
    polycubeTriMesh.GetFaceAndNeighborFaces();
    polycubeTriMesh.GetFacePatches_N();
    ////////////////////////////////////////////
    polycubeTriMesh.GetFacePatches();
    polycubeTriMesh.GetEdgePatches();
    polycubeTriMesh.GetVertexPatches();
    polycubeTriMesh.GetPatchesLabel();
    ////////////////////////////////////////////
    polycubeTriMesh.GetPatches();
    polycubeTriMesh.SortPatches();
    polycubeTriMesh.ModifyPatchesPosition();

    polycubeTriMesh.AssignPatches();

    std::cout << "-------------Extracting Polycube Hex Mesh -------------- " << std::endl;
    const glm::vec3 cellScale(x, y, z);
    polycubeTriMesh.GetMaxMinCoordinates_Patch();
    polycubeTriMesh.DividePolycube_aniso(polycubeTriMesh.m_minVertex, polycubeTriMesh.m_vecPolycubeHexVertex,
            polycubeTriMesh.m_vecPolycubeHexCell, cellScale, delta);
    polycubeTriMesh.removeInvalidHexahedronsInBox(polycubeTriMesh.m_min_distance, dir);
    polycubeTriMesh.FixMesh();
    polycubeTriMesh.WriteHexahedralmesh(strOutputPolyHexFileName.c_str());

    MeshFileWriter ajustedmeshFileWriter(polycubeTriMesh, strAjustedPolycubeTriMeshFilename.c_str());
    ajustedmeshFileWriter.WriteFile();

    return 0;
}

int ParseArg(int argc, char* argv[],
        std::string& strInputPolyTriFileName,
        std::string& strOutputPolyHexFileName,
        std::string& strAjustedPolycubeTriMeshFilename,
        std::string& strPolycubeNormalTriMeshFileName,
        std::string& strDir,
        DIRECTION& dir,
        float& delta,
        int& iterNum,
        float& x,
        float& y,
        float& z
        )
{
    try
    {
#ifdef CGAL_USE_BOOST_PROGRAM_OPTIONS
        po::options_description desc("Allowed options");
        desc.add_options()
            ("help,h", "produce help message.")
            ("thread,j", po::value<unsigned int>(&THREAD_NUM)->default_value(4), "threads number")
            // ########## input #########
            // -p
            ("polycube-input-tri-file,p", po::value<std::string>(&strInputPolyTriFileName)->default_value("polycube.tri.off"),
                "name of the polycube input triangle off file")
            // ########## output #########
            // -P
            ("polycube-normalize-tri-file,P", po::value<std::string>(&strPolycubeNormalTriMeshFileName)->default_value("polycube.normal.tri.vtk"),
                "name of the orig-normal-tri-file vtk file")
            // -n
            ("out-polycube-hex-file,n", po::value<std::string>(&strOutputPolyHexFileName)->default_value("polycube.aniso.hex.vtk"),
                "name of the output polycube hexaheral vtk file")
            // -A
            ("Ajust-tri-file,A", po::value<std::string>(&strAjustedPolycubeTriMeshFilename)->default_value("adjusted.polycube.tri.vtk"),
                "name of the ajust tri polycube vtk file")
            // ########## parameters #########
            // -s
            ("hex-parameterization-size,s", po::value<float>(&delta)->default_value(0.02), "cube size, 0.001f <= size <= 0.1f")
            // -w
            ("offset,w", po::value<float>(&__plus)->default_value(0.001), "0.01s <= offset <= s")
            // -r
            ("ray-direction,r", po::value<std::string>(&strDir)->default_value("Z"), "ray-direction, X, Y, or Z")
            // -l
            ("laplacian,l", po::value<int>(&iterNum)->default_value(1), "laplacian smooth iteration number")
            // -U
            ("Upper-size,U", po::value<double>(&g_upper_size)->default_value(0.008), "cube upper size")
            // -L
            ("Lower-size,L",po::value<double>(&g_lower_size)->default_value(0.004), "cube lower size")
            ("x-aniso,x", po::value<float>(&x)->default_value(1.0), "x-aniso ratio respect to cube size")
            ("y-aniso,y", po::value<float>(&y)->default_value(1.0), "y-aniso ratio respect to cube size")
            ("z-aniso,z", po::value<float>(&z)->default_value(1.0), "z-aniso ratio respect to cube size")
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

        if (vm.count("polycube-input-tri-file"))
            std::cout << "InputPolyTriFileName = " << strInputPolyTriFileName << std::endl;
        if (vm.count("out-polycube-hex-file"))
            std::cout << "OutputPolyHexFileName = " << strOutputPolyHexFileName << std::endl;
        if (vm.count("Ajust-hex-file"))
            std::cout << "AjustedPolycubeTriMeshFilename = " << strAjustedPolycubeTriMeshFilename << std::endl;
        if (vm.count("polycube-normalize-tri-file"))
            std::cout << "PolycubeNormalTriMeshFileName = " << strPolycubeNormalTriMeshFileName << std::endl;

        if (vm.count("thread"))
            std::cout << "threadnum = " << THREAD_NUM << std::endl;
        if (vm.count("hex-parameterization-size"))
            std::cout << "hexsize = " << delta << std::endl;
        if (vm.count("ray-direction"))
            std::cout << "ray direction = " << strDir << std::endl;
        if (vm.count("laplacian"))
            std::cout << "smooth iteration number = " << iterNum << std::endl;
        if (vm.count("offset"))
            std::cout << "offset = " << __plus << std::endl;
        if (vm.count("x-aniso"))
            std::cout << "x-aniso = " << x << std::endl;
        if (vm.count("y-aniso"))
            std::cout << "y-aniso = " << y << std::endl;
        if (vm.count("z-aniso"))
            std::cout << "z-aniso = " << z << std::endl;

#else
        strInputPolyTriFileName = "polycube.tri.vtk";
        strOutputPolyHexFileName = "polycube.hex.vtk";
#endif
    } catch (std::exception& e)
    {
        std::cerr << "error: " << e.what() << "\n";
        return 1;
    } catch (...)
    {
        std::cerr << "Exception of unknown type!\n";
    }
}

