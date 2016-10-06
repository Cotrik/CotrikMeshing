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
        DIRECTION& dir,
        float& delta,
        int& iterNum,
        int& Func)
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
            ("out-polycube-hex-file,n", po::value<std::string>(&strOutputPolyHexFileName)->default_value("polycube.hex.vtk"),
                "name of the output polycube hexaheral vtk file")
            // -A
            ("Ajust-tri-file,A", po::value<std::string>(&strAjustedPolycubeTriMeshFilename)->default_value("ajusted.polycube.tri.vtk"),
                "name of the ajust tri polycube vtk file")
            // ########## parameters #########
            // -s
            ("hex-parameterization-size,s", po::value<float>(&delta)->default_value(0.02),
                "cube size, 0.001f <= size <= 0.1f")
            // -w
            ("offset,w", po::value<float>(&__plus)->default_value(0.001),
                "0.01s <= offset <= s")
            // -F
            ("Function,F", po::value<int>(&Func)->default_value(1),
                "0 Hexify, 1 Magnify")
            // -r
            ("ray-direction,r", po::value<std::string>(&strDir)->default_value("Z"),
                "ray-direction, X, Y, or Z")
            // -l
            ("laplacian,l", po::value<int>(&iterNum)->default_value(1),
                "laplacian smooth iteration number")
            // -U
            ("Upper-size,U", po::value<double>(&g_upper_size)->default_value(0.008), "cube upper size")
            // -L
            ("Lower-size,L",po::value<double>(&g_lower_size)->default_value(0.004), "cube lower size")
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

#else
        strInputPolyTriFileName = "polycube.tri.off";
        strOutputPolyHexFileName = "poly.hex.vtk";
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
    int Func = 1;
    ParseArg(argc, argv, strInputPolyTriFileName,
            strOutputPolyHexFileName,
            strAjustedPolycubeTriMeshFilename,
            strPolycubeNormalTriMeshFileName,
            strDir, dir, delta, iterNum, Func);

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

    for (int i = 0; i < polycubeTriMesh.m_patches.size(); i++)
    {
        polycubeTriMesh.facePatches.at(i) = polycubeTriMesh.m_patches.at(i).m_f;
        for (int j = 0; j < polycubeTriMesh.facePatches.at(i).size(); j++)
        {
            polycubeTriMesh.faceLabel.at(polycubeTriMesh.facePatches.at(i).at(j)) = i + 1;
        }
        polycubeTriMesh.edgePatches.at(i) = polycubeTriMesh.m_patches.at(i).m_e;
        polycubeTriMesh.vertexPatches.at(i) = polycubeTriMesh.m_patches.at(i).m_v;
    }

    polycubeTriMesh.GetMaxMinCoordinates_Patch();
    polycubeTriMesh.AdjustPolycubeSurfaceToHoldAllSmallCubesVolume1_Reverse(delta, polycubeTriMesh.m_minVertex);
    polycubeTriMesh.AddCompensation();

    polycubeTriMesh.GetPatches();
    polycubeTriMesh.ModifyPatchesPosition();
    polycubeTriMesh.GetPatches();

    std::cout << "-------------Extracting Polycube Hex Mesh -------------- " << std::endl;
    polycubeTriMesh.GetMaxMinCoordinates_Patch();
    polycubeTriMesh.DividePolycube(polycubeTriMesh.m_minVertex, polycubeTriMesh.m_vecPolycubeHexVertex, polycubeTriMesh.m_vecPolycubeHexCell);
    polycubeTriMesh.removeInvalidHexahedronsInBox_ref((std::vector<Vertex>&)polycubeTriMeshFileReader.GetMesh().V, polycubeTriMesh.m_min_distance, dir);
    polycubeTriMesh.FixMesh();
    polycubeTriMesh.WriteHexahedralmesh(strOutputPolyHexFileName.c_str());

    MeshFileWriter ajustedmeshFileWriter(polycubeTriMesh, strAjustedPolycubeTriMeshFilename.c_str());
    ajustedmeshFileWriter.WriteFile();

    return 0;
}
