/*
 * hexgen_auto.cpp
 *
 *  Created on: Jun 8, 2015
 *      Author: cotrik
 */

#include "include.h"

static void GetDir(const std::string& strDir, DIRECTION& dir)
{
    if (strDir.empty())     return;
    if (strDir[0] == 'X')       dir = X_AXIS;
    else if (strDir[0] == 'Y')  dir = Y_AXIS;
    else if (strDir[0] == 'Z')  dir = Z_AXIS;
}

static int ParseArg(int argc, char* argv[],
        std::string& strInputOrigTetFileName,
        std::string& strInputPolyTetFileName,
        std::string& strOutputPolyHexFileName,
        std::string& strOutputResultHexFileName,

        std::string& strAjustedPolycubeTetMeshFilename,
        std::string& strOrigNormalTetMeshFileName,
        std::string& strPolycubeNormalTetMeshFileName,

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
            ("out-polycube-hex-file,n", po::value<std::string>(&strOutputPolyHexFileName)->default_value("poly.hex.vtk"),
                "name of the output polycube hexaheral vtk/mesh file")
            // -o
            ("out-hex-file,o", po::value<std::string>(&strOutputResultHexFileName)->default_value("result.hex.vtk"),
                "name of the output hexaheral vtk/mesh file")
            // -A
            ("Ajust-hex-file,A", po::value<std::string>(&strAjustedPolycubeTetMeshFilename)->default_value("ajusted.polycube.tet.vtk"),
                "name of the ajust tet polycube vtk/mesh file")
            // ########## parameters #########
            // -s
            ("hex-parameterization-size,s", po::value<float>(&delta)->default_value(0.02),
                "cube size, 0.001f <= size <= 0.1f")
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
}

static void AlignPatch(PolycubeMesh& polycubeTetMesh, PolycubeMesh& polycubeHexMesh)
{
    for (int i = 0; i < polycubeTetMesh.m_vertex_patches_num_x; i++)
    {
        std::vector<unsigned long>& vetexPatch = polycubeTetMesh.m_patches.at(i).m_v;
        double pos = 0.0;
        const size_t vetexPatchSize = vetexPatch.size();
        for (int j = 0; j < vetexPatchSize; j++)
        {
            const Vertex& v = polycubeTetMesh.V.at(vetexPatch.at(j));
            pos += v.x;
        }
        if (vetexPatchSize != 0)
        {
            pos /= vetexPatchSize;
            std::vector<unsigned long>& vetexHexPatch = polycubeHexMesh.m_patches.at(i).m_v;
            const size_t vetexHexPatchSize = vetexHexPatch.size();
            double oldx = 0;
            for (int j = 0; j < vetexHexPatchSize; j++)
            {
                Vertex& v = polycubeHexMesh.V.at(vetexHexPatch.at(j));
                oldx = v.x;
                v.x = pos;
            }
            for (int j = 0; j < polycubeHexMesh.V.size(); j++)
            {
                Vertex& v = polycubeHexMesh.V.at(j);
                if (fabs(v.x - oldx) < 5e-5)
                    v.x = pos;
            }
        }
    }
    for (int i = polycubeTetMesh.m_vertex_patches_num_x; i < polycubeTetMesh.m_vertex_patches_num_y; i++)
    {
        std::vector<unsigned long>& vetexPatch = polycubeTetMesh.m_patches.at(i).m_v;
        double pos = 0.0;
        const size_t vetexPatchSize = vetexPatch.size();
        for (int j = 0; j < vetexPatchSize; j++)
        {
            const Vertex& v = polycubeTetMesh.V.at(vetexPatch.at(j));
            pos += v.y;
        }
        if (vetexPatchSize != 0)
        {
            pos /= vetexPatchSize;
            std::vector<unsigned long>& vetexHexPatch = polycubeHexMesh.m_patches.at(i).m_v;
            const size_t vetexHexPatchSize = vetexHexPatch.size();
            double oldy = 0;
            for (int j = 0; j < vetexHexPatchSize; j++)
            {
                Vertex& v = polycubeHexMesh.V.at(vetexHexPatch.at(j));
                oldy = v.y;
                v.y = pos;
            }
            for (int j = 0; j < polycubeHexMesh.V.size(); j++)
            {
                Vertex& v = polycubeHexMesh.V.at(j);
                if (fabs(v.y - oldy) < 5e-5)
                    v.y = pos;
            }
        }
    }
    for (int i = polycubeTetMesh.m_vertex_patches_num_y; i < polycubeTetMesh.m_vertex_patches_num_z; i++)
    {
        std::vector<unsigned long>& vetexPatch = polycubeTetMesh.m_patches.at(i).m_v;
        double pos = 0.0;
        const size_t vetexPatchSize = vetexPatch.size();
        for (int j = 0; j < vetexPatchSize; j++)
        {
            const Vertex& v = polycubeTetMesh.V.at(vetexPatch.at(j));
            pos += v.z;
        }
        if (vetexPatchSize != 0)
        {
            pos /= vetexPatchSize;
            std::vector<unsigned long>& vetexHexPatch = polycubeHexMesh.m_patches.at(i).m_v;
            const size_t vetexHexPatchSize = vetexHexPatch.size();
            double oldz = 0;
            for (int j = 0; j < vetexHexPatchSize; j++)
            {
                Vertex& v = polycubeHexMesh.V.at(vetexHexPatch.at(j));
                oldz = v.z;
                v.z = pos;
            }
            for (int j = 0; j < polycubeHexMesh.V.size(); j++)
            {
                Vertex& v = polycubeHexMesh.V.at(j);
                if (fabs(v.z - oldz) < 5e-5)
                    v.z = pos;
            }
        }
    }
}

int main(int argc, char* argv[])
{
    std::string strInputOrigTetFileName;
    std::string strInputPolyTetFileName;
    std::string strOutputPolyHexFileName;
    std::string strOutputResultHexFileName;

    std::string strAjustedPolycubeTetMeshFilename;
    std::string strOrigNormalTetMeshFileName;
    std::string strPolycubeNormalTetMeshFileName;

    std::string strDir;
    DIRECTION dir = X_AXIS;
    float delta = 0.02f;
    int iterNum = 0;
    int Func = 1;
    ParseArg(argc, argv, strInputOrigTetFileName, strInputPolyTetFileName,
            strOutputPolyHexFileName, strOutputResultHexFileName,
            strAjustedPolycubeTetMeshFilename, strOrigNormalTetMeshFileName,
            strPolycubeNormalTetMeshFileName,
            strDir, dir, delta, iterNum, Func);
    //////////////////////////////////////////
    // Read Orig Tet Mesh
    //////////////////////////////////////////
    std::cout << "-------------Processing Orig Tet Mesh -------------- " << std::endl;
    MeshFileReader orgTetMeshFileReader(strInputOrigTetFileName.c_str());
    Mesh origTetMesh(orgTetMeshFileReader.GetMesh());
    origTetMesh.ExtractSurface();
    origTetMesh.NormalizeCoordinateValue();
    MeshFileWriter normalTetMeshFileWriter(origTetMesh, strOrigNormalTetMeshFileName.c_str());
    normalTetMeshFileWriter.WriteFile();

    std::cout << "Writing surface vtk file" << std::endl;
    MeshFileWriter surfaceTetMeshFileWriter(origTetMesh.V, origTetMesh.surface, "orig.tet.surface.vtk", TETRAHEDRA);
    surfaceTetMeshFileWriter.WriteFile();
    std::cout << "Finish writing surface vtk file" << std::endl;
    //////////////////////////////////////////
    // Read Polycube Tet Mesh
    //////////////////////////////////////////
    std::cout << "-------------Processing Polycube Tet Mesh -------------- " << std::endl;
    MeshFileReader polycubeTetMeshFileReader(strInputPolyTetFileName.c_str());
    PolycubeMesh polycubeTetMesh(polycubeTetMeshFileReader.GetMesh());
    polycubeTetMesh.ExtractSurface();
    polycubeTetMesh.NormalizeCoordinateValue();
    polycubeTetMesh.GetCorners();
    MeshFileWriter polycubeNormalTetMeshFileWriter(polycubeTetMesh, strPolycubeNormalTetMeshFileName.c_str());
    polycubeNormalTetMeshFileWriter.WriteFile();

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    polycubeTetMesh.GetVertexInfo();
    for (size_t i = 0; i < origTetMesh.V.size(); i++)
        origTetMesh.V.at(i).vinfo = polycubeTetMesh.V.at(i).vinfo;

    MeshFileWriter ajustedmeshFileWriter(polycubeTetMesh, strAjustedPolycubeTetMeshFilename.c_str());
    ajustedmeshFileWriter.WriteFile();
    strAjustedPolycubeTetMeshFilename = strInputPolyTetFileName + std::string(".ajusted.vtk");
    ////////////////////////////////////////////////////////////////////////////////
    polycubeTetMesh.GetMaxMinCoordinates();
    polycubeTetMesh.ExtractSurface();
    polycubeTetMesh.GetFaceType();
    polycubeTetMesh.GetFaceAndNeighborFaces();
    polycubeTetMesh.GetFacePatches_N();
    polycubeTetMesh.GetNormalOfSurfaceFaces();
    polycubeTetMesh.GetNormalOfSurfaceVertices();
    polycubeTetMesh.GetFacePatches();
    polycubeTetMesh.GetEdgePatches();
    polycubeTetMesh.GetVertexPatches();
    //////////////////////////////
    polycubeTetMesh.GetPatches();
    polycubeTetMesh.SortPatches();
    polycubeTetMesh.ModifyPatchesPosition();

    std::cout << "-------------Extracting Polycube Hex Mesh -------------- " << std::endl;
    polycubeTetMesh.GetMaxMinCoordinates_Patch();
    polycubeTetMesh.DividePolycube(polycubeTetMesh.m_minVertex,
            polycubeTetMesh.m_vecPolycubeHexVertex, polycubeTetMesh.m_vecPolycubeHexCell);
    polycubeTetMesh.removeInvalidHexahedronsInBox(polycubeTetMesh.m_min_distance, dir);
    polycubeTetMesh.FixMesh();
    polycubeTetMesh.WriteHexahedralmesh(strOutputPolyHexFileName.c_str());

    std::cout << "-------------Processing Polycube Hex Mesh -------------- " << std::endl;
    MeshFileReader polycubeHexMeshFileReader(strOutputPolyHexFileName.c_str());
    PolycubeMesh polycubeHexMesh(polycubeHexMeshFileReader.GetMesh());
    //polycubeHexMesh.NormalizeCoordinateValue();
    polycubeHexMesh.ExtractSurface();
    polycubeTetMesh.GetFaceType();
    polycubeHexMesh.GetMaxMinCoordinates();
    polycubeHexMesh.GetNormalOfSurfaceFaces();
    polycubeHexMesh.GetNormalOfSurfaceVertices();
    polycubeHexMesh.GetCorners();
    polycubeHexMesh.GetVertexInfo();

    polycubeHexMesh.LabelSurfaceFace();
    polycubeHexMesh.GetFacePatches();
    polycubeHexMesh.GetEdgePatches();
    polycubeHexMesh.GetVertexPatches();
    polycubeHexMesh.GetPatchesLabel();
    //GetFaceType(polycubeHexMesh);
    polycubeHexMesh.GetPatchMinDistance();

    polycubeTetMesh.CheckPatchesConnection();
    polycubeHexMesh.CheckPatchesConnection();
    polycubeHexMesh.SortPatches(polycubeTetMesh);
    //Align(polycubeTetMesh, polycubeHexMesh);
    AlignPatch(polycubeTetMesh, polycubeHexMesh);

    MeshFileWriter polycubeFileWriter(polycubeHexMesh.V, polycubeHexMesh.C, strOutputPolyHexFileName.c_str(), HEXAHEDRA);
    polycubeFileWriter.WriteFile();

//    MeshFileWriter polycubeMeshFileWriter(polycubeHexMesh.V, polycubeHexMesh.C, "polycube.hex.mesh", HEXAHEDRA);
//    polycubeMeshFileWriter.WriteFile();

    MeshFileWriter meshFileWriter1(polycubeTetMesh.V, polycubeTetMesh.C, strAjustedPolycubeTetMeshFilename.c_str(), TETRAHEDRA);
    meshFileWriter1.WriteFile();
    {
        polycubeHexMesh.GetMaxMinCoordinates();
        polycubeHexMesh.AdjustPolycubeSurfaceToHoldAllSmallCubesVolume_quad(delta, polycubeHexMesh.m_minVertex);
        polycubeHexMesh.AddCompensation_quad();

        MeshFileWriter polycubeFileWriter2(polycubeHexMesh.V, polycubeHexMesh.C, strOutputPolyHexFileName.c_str(), HEXAHEDRA);
        polycubeFileWriter2.WriteFile();

        polycubeHexMesh.IsAllVerticesInsideMesh(polycubeTetMesh);
    }
    std::cout << "-------------Extracting Hex Mesh -------------- " << std::endl;
    std::vector<Vertex> vecHexVertex;
    std::vector<Cell> vecHexCell;
    polycubeTetMesh.m_vecPolycubeHexVertex = polycubeHexMesh.V;

    polycubeTetMesh.CreateHexMesh(origTetMesh, vecHexVertex, vecHexCell);

    MeshFileWriter meshFileWriter(vecHexVertex, vecHexCell, strOutputResultHexFileName.c_str(), HEXAHEDRA);
    meshFileWriter.WriteFile();

    return 0;
}

