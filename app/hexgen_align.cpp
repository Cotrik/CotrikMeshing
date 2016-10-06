/*
 * hexgen_align.cpp
 *
 *  Created on: Jun 7, 2015
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

int HexGen_align(int argc, char* argv[])
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

    //////////////////////////////////////////
    // Read Orig Tet Mesh
    //////////////////////////////////////////
    MeshFileReader orgTetMeshFileReader(strInputOrigTetFileName.c_str());
    Mesh origTetMesh(orgTetMeshFileReader.GetMesh());
    origTetMesh.ExtractSurface();
    origTetMesh.NormalizeCoordinateValue();
    MeshFileWriter normalTetMeshFileWriter(origTetMesh, strOrigNormalTetMeshFileName.c_str());
    normalTetMeshFileWriter.WriteMeshFile();
    //////////////////////////////////////////
    // Read Polycube Tet Mesh
    //////////////////////////////////////////
    MeshFileReader polycubeTetMeshFileReader(strInputPolyTetFileName.c_str());
    PolycubeMesh polycubeTetMesh(polycubeTetMeshFileReader.GetMesh());
    polycubeTetMesh.ExtractSurface();
    polycubeTetMesh.NormalizeCoordinateValue();
    polycubeTetMesh.GetCorners();
    MeshFileWriter polycubeNormalTetMeshFileWriter(polycubeTetMesh, strPolycubeNormalTetMeshFileName.c_str());
    polycubeNormalTetMeshFileWriter.WriteMeshFile();

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    polycubeTetMesh.GetVertexInfo();
    for (int i = 0; i < origTetMesh.V.size(); i++)
        origTetMesh.V.at(i).vinfo = polycubeTetMesh.V.at(i).vinfo;

    MeshFileWriter ajustedmeshFileWriter(polycubeTetMesh, strAjustedPolycubeTetMeshFilename.c_str());
    ajustedmeshFileWriter.WriteMeshFile();
    strAjustedPolycubeTetMeshFilename = strInputPolyTetFileName + std::string(".ajusted.vtk");
//  polycubeTetMesh.GenerateHexahedralMesh_align(origTetMesh, strOutputPolyHexFileName.c_str(), strOutputResultHexFileName.c_str(), delta, dir);
    polycubeTetMesh.GenerateSmallCubes_plus(delta, NULL, __plus);
    polycubeTetMesh.removeInvalidHexahedronsInBox(delta, dir);
    polycubeTetMesh.WriteHexahedralmesh(strOutputPolyHexFileName.c_str());

    ////////////////////////////////////////////////////////////////////////////////
    polycubeTetMesh.NormalizeCoordinateValue();
    polycubeTetMesh.ExtractSurface();
    polycubeTetMesh.GetNormalOfSurfaceFaces();
    polycubeTetMesh.GetNormalOfSurfaceVertices();
    polycubeTetMesh.GetCorners();

    polycubeTetMesh.LabelSurfaceFace();
    polycubeTetMesh.GetFacePatches();
    polycubeTetMesh.GetEdgePatches();
    polycubeTetMesh.GetVertexPatches();
    polycubeTetMesh.GetPatchesLabel();
    polycubeTetMesh.GetFaceType();

//    for (int i = 0; i < polycubeTetMesh.surface.size(); i++)
//    {
//        const Face& f = polycubeTetMesh.surface.at(i);
//        //PFT.push_back(GetFaceType(polycubeTetMesh.N_F.at(i)));
//        const Vertex* pVertex0 = &polycubeTetMesh.V[f.at(0)];
//        const Vertex* pVertex1 = &polycubeTetMesh.V[f.at(1)];
//        const Vertex* pVertex2 = &polycubeTetMesh.V[f.at(2)];
//
//        const glm::vec3 v0(pVertex0->x, pVertex0->y, pVertex0->z);
//        const glm::vec3 v1(pVertex1->x, pVertex1->y, pVertex1->z);
//        const glm::vec3 v2(pVertex2->x, pVertex2->y, pVertex2->z);
//
//        const glm::vec3 v10 = v0 - v1;
//        const glm::vec3 v12 = v2 - v1;
//        const glm::vec3 normal = glm::cross(v12, v10);
//
//        polycubeTetMesh.faceType.push_back(GetFaceType(normal));
//    }

    for (int i = 0; i < polycubeTetMesh.patchLabel.size(); i++)
    {
        std::cout << "patch " << i << " label: ";
        for (int j = 0; j < polycubeTetMesh.patchLabel.at(i).size(); j++)
        {
            std::cout << polycubeTetMesh.patchLabel.at(i).at(j) << " ";
        }
        std::cout << std::endl;
    }
    std::cout << "---------------------------------------------------------------- " << std::endl;
    MeshFileReader polycubeHexMeshFileReader(strOutputPolyHexFileName.c_str());
    PolycubeMesh polycubeHexMesh(polycubeHexMeshFileReader.GetMesh());
    polycubeHexMesh.NormalizeCoordinateValue();
    polycubeHexMesh.ExtractSurface();
    polycubeHexMesh.GetNormalOfSurfaceFaces();
    polycubeHexMesh.GetNormalOfSurfaceVertices();
    polycubeHexMesh.GetCorners();

    polycubeHexMesh.LabelSurfaceFace();
    polycubeHexMesh.GetFacePatches();
    polycubeHexMesh.GetEdgePatches();
    polycubeHexMesh.GetVertexPatches();
    polycubeHexMesh.GetPatchesLabel();

    for (int i = 0; i < polycubeHexMesh.patchLabel.size(); i++)
    {
        std::cout << "patch " << i << " label: ";
        for (int j = 0; j < polycubeHexMesh.patchLabel.at(i).size(); j++)
        {
            std::cout << polycubeHexMesh.patchLabel.at(i).at(j) << " ";
        }
        std::cout << std::endl;
    }

    std::vector<double> patchPos;
    for (int i = 0; i < polycubeTetMesh.patchLabel.size(); i++)
    {
        std::vector<unsigned long>& facePatch = polycubeTetMesh.facePatches.at(i);
        std::vector<unsigned long>& vetexPatch = polycubeTetMesh.vertexPatches.at(i);
        double x = 0.0, y = 0.0, z = 0.0;
        const size_t vetexPatchSize = vetexPatch.size();
        for (int j = 0; j < vetexPatchSize; j++)
        {
            const Vertex& v = polycubeTetMesh.V.at(vetexPatch.at(j));
            x += v.x;
            y += v.y;
            z += v.z;
        }
        if (vetexPatchSize != 0)
        {
            x /= vetexPatchSize;
            y /= vetexPatchSize;
            z /= vetexPatchSize;
            //----------------------------------------
            std::vector<unsigned long>& faceHexPatch = polycubeHexMesh.facePatches.at(i);
            std::vector<unsigned long>& vetexHexPatch = polycubeHexMesh.vertexPatches.at(i);
            const size_t vetexHexPatchSize = vetexHexPatch.size();
            const FACE_TYPE faceType = polycubeTetMesh.faceType.at(facePatch.at(0));

            if (faceType == FACE_X)
            {
                patchPos.push_back(x);
                double oldx = 0;
                for (int j = 0; j < vetexHexPatchSize; j++)
                {
                    Vertex& v = polycubeHexMesh.V.at(vetexHexPatch.at(j));
                    oldx = v.x;
                    v.x = x;
                }
                for (int j = 0; j < polycubeHexMesh.V.size(); j++)
                {
                    Vertex& v = polycubeHexMesh.V.at(j);
                    if (fabs(v.x - oldx) < 5e-4)
                        v.x = x;
                }
            }
            else if (faceType == FACE_Y)
            {
                patchPos.push_back(y);
                double oldy = 0;
                for (int j = 0; j < vetexHexPatchSize; j++)
                {
                    Vertex& v = polycubeHexMesh.V.at(vetexHexPatch.at(j));
                    oldy = v.y;
                    v.y = y;
                }
                for (int j = 0; j < polycubeHexMesh.V.size(); j++)
                {
                    Vertex& v = polycubeHexMesh.V.at(j);
                    if (fabs(v.y - oldy) < 5e-4)
                        v.y = y;
                }
            }
            else if (faceType == FACE_Z)
            {
                patchPos.push_back(z);
                double oldz = 0;
                for (int j = 0; j < vetexHexPatchSize; j++)
                {
                    Vertex& v = polycubeHexMesh.V.at(vetexHexPatch.at(j));\
                    oldz = v.z;
                    v.z = z;
                }
                for (int j = 0; j < polycubeHexMesh.V.size(); j++)
                {
                    Vertex& v = polycubeHexMesh.V.at(j);
                    if (fabs(v.z - oldz) < 5e-4)
                        v.z = z;
                }
            }
        }
    }
    double min_distance = 1;
    for (int i = 1; i < polycubeTetMesh.m_xLabel; i++)
    {
        double distance = patchPos.at(i) - patchPos.at(i - 1);
        std::cout << "patch distance = " << patchPos.at(i) << " - " << patchPos.at(i - 1) << " = " << distance << std::endl;
        if (distance < min_distance)
            min_distance = distance;
    }
    for (int i = polycubeTetMesh.m_xLabel + 1; i < polycubeTetMesh.m_yLabel; i++)
    {
        double distance = patchPos.at(i) - patchPos.at(i - 1);
        std::cout << "patch distance = " << patchPos.at(i) << " - " << patchPos.at(i - 1) << " = " << distance << std::endl;
        if (distance < min_distance)
            min_distance = distance;
    }
    for (int i = polycubeTetMesh.m_yLabel + 1; i < polycubeTetMesh.m_zLabel; i++)
    {
        double distance = patchPos.at(i) - patchPos.at(i - 1);
        std::cout << "patch distance = " << patchPos.at(i) << " - " << patchPos.at(i - 1) << " = " << distance << std::endl;
        if (distance < min_distance)
            min_distance = distance;
    }
    cout << "patch min_distance = " << min_distance << std::endl;

    MeshFileWriter polycubeFileWriter(polycubeHexMesh.V, polycubeHexMesh.C, strOutputPolyHexFileName.c_str());
    polycubeFileWriter.WriteMeshFile();
    ///////////////////////////////////////////////////////////////////////////////
    polycubeTetMesh.GetMaxMinCoordinates();
    polycubeTetMesh.AdjustPolycubeSurfaceToHoldAllSmallCubesVolume1(delta, polycubeTetMesh.m_minVertex);
    polycubeTetMesh.AddCompensation();
    MeshFileWriter meshFileWriter1(polycubeTetMesh.V, polycubeTetMesh.C, strAjustedPolycubeTetMeshFilename.c_str());
    meshFileWriter1.WriteMeshFile();

//  polycubeHexMesh.GetMaxMinCoordinates();
//  polycubeHexMesh.AdjustPolycubeSurfaceToHoldAllSmallCubesVolume_quad(delta, polycubeTetMesh.m_minVertex);
//  polycubeHexMesh.AddCompensation_quad();

    ///////////////////////////////////////////////////////////////////////////////
    std::vector<Vertex> vecHexVertex;
    std::vector<Cell> vecHexCell;
    polycubeTetMesh.m_vecPolycubeHexVertex = polycubeHexMesh.V;
    polycubeTetMesh.CreateHexMesh(origTetMesh, vecHexVertex, vecHexCell);

    MeshFileWriter meshFileWriter(vecHexVertex, vecHexCell, strOutputResultHexFileName.c_str());
    meshFileWriter.WriteMeshFile();

    return 0;
}


int main(int argc, char* argv[])
{
	return HexGen_align(argc, argv);
}


