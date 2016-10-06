/*
 * polycubehexgen_round_aniso.cpp
 *
 *  Created on: Nov 3, 2015
 *      Author: cotrik
 */

#include "include.h"
#include "igl/colon.h"
#include "igl/directed_edge_orientations.h"
#include "igl/directed_edge_parents.h"
#include "igl/forward_kinematics.h"
#include "igl/PI.h"
#include "igl/lbs_matrix.h"
#include "igl/deform_skeleton.h"
#include "igl/dqs.h"
#include "igl/readDMAT.h"
#include "igl/readOFF.h"
#include "igl/arap.h"

#include <Eigen/Geometry>
#include <Eigen/StdVector>
#include <Eigen/Core>
#include <Eigen/Eigen>
#include <Eigen/Dense>
#include <Eigen/Sparse>
//#include <Eigen/MPRealSupport>
#include <Eigen/Cholesky>
#include <Eigen/LU>
#include <Eigen/SVD>
#include <Eigen/SparseCore>
static int ParseArg(int argc, char* argv[],
        std::string& strInputPolyTetFileName,
        std::string& strOutputPolyHexFileName,
        std::string& strOutputAdjustedPolyTetMeshFilename,
        float& delta, float& x, float& y, float& z
        );
void update_vertices(Mesh& mesh);
void set_handle_points(const Mesh& mesh, const std::vector<unsigned long>& selectedPointIds);
int main(int argc, char* argv[])
{
    std::string strInputPolyTetFileName;
    std::string strOutputPolyHexFileName;
    std::string strOutputAdjustedPolyTetMeshFilename;
    DIRECTION dir = Z_AXIS;
    float delta = 0.02f;
    float x = 1.0f, y = 1.0f, z = 1.0f;
    ParseArg(argc, argv, strInputPolyTetFileName,
            strOutputPolyHexFileName,
            strOutputAdjustedPolyTetMeshFilename,
            delta, x, y, z);

    //////////////////////////////////////////
    // Read Polycube Tet Mesh
    //////////////////////////////////////////
    std::cout << "------------- Processing Polycube Triangle Mesh -------------- " << std::endl;
    MeshFileReader polycubeTetMeshFileReader(strInputPolyTetFileName.c_str());
    PolycubeMesh polycubeTetMesh(polycubeTetMeshFileReader.GetMesh());
    polycubeTetMesh.ExtractSurface();
    polycubeTetMesh.NormalizeCoordinateValue();
    polycubeTetMesh.GetCorners();
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    polycubeTetMesh.GetVertexInfo();
    polycubeTetMesh.GetMaxMinCoordinates();
    polycubeTetMesh.GetFaceType();
    polycubeTetMesh.GetNormalOfSurfaceFaces();
    polycubeTetMesh.GetNormalOfSurfaceVertices();
    polycubeTetMesh.GetCorners();
    polycubeTetMesh.GetVertexInfo();
    ////////////////////////////////////////////
    polycubeTetMesh.GetFaceAndNeighborFaces();
    polycubeTetMesh.GetFacePatches_N();
    ////////////////////////////////////////////
    polycubeTetMesh.GetFacePatches();
    polycubeTetMesh.GetEdgePatches();
    polycubeTetMesh.GetVertexPatches();
    polycubeTetMesh.GetPatchesLabel();
    ////////////////////////////////////////////
    polycubeTetMesh.GetPatches();
    polycubeTetMesh.SortPatches();
    polycubeTetMesh.ModifyPatchesPosition();

    polycubeTetMesh.AssignPatches();
    MeshFileWriter adjustedmeshFileWriter1(polycubeTetMesh.V, polycubeTetMesh.surface, "n.vtk", TRIANGLE);
    adjustedmeshFileWriter1.WriteFile();
    /////////////////////////////////////////////////////////////
    polycubeTetMesh.GetMaxMinCoordinates_Patch();
    // Get selectedVerticesIds
    std::vector<unsigned long > selectedVerticesIds;
    for (int i = 0; i < polycubeTetMesh.surface.size(); i++)
        for (int j = 0; j < 3; j++)
            selectedVerticesIds.push_back(polycubeTetMesh.surface.at(i).at(j));
    std::sort(selectedVerticesIds.begin(), selectedVerticesIds.end());
    std::vector<unsigned long>::iterator iter = std::unique(selectedVerticesIds.begin(), selectedVerticesIds.end());
    selectedVerticesIds.resize(std::distance(selectedVerticesIds.begin(), iter));
    set_handle_points(polycubeTetMesh, selectedVerticesIds);
    polycubeTetMesh.ModifyPatchesPosition_ROUND(delta);
    update_vertices(polycubeTetMesh);
    polycubeTetMesh.GetPatches();
    polycubeTetMesh.ModifyPatchesPosition();
    polycubeTetMesh.GetPatches();
    MeshFileWriter adjustedmeshFileWriter(polycubeTetMesh, strOutputAdjustedPolyTetMeshFilename.c_str());
    adjustedmeshFileWriter.WriteFile();
    /////////////////////////////////////////////////////////////
    std::cout << "-------------Extracting Polycube Hex Mesh -------------- " << std::endl;
    const glm::vec3 cellScale(x, y, z);
    polycubeTetMesh.GetMaxMinCoordinates_Patch();
    polycubeTetMesh.DividePolycube_aniso(polycubeTetMesh.m_minVertex, polycubeTetMesh.m_vecPolycubeHexVertex,
            polycubeTetMesh.m_vecPolycubeHexCell, cellScale, delta);
    polycubeTetMesh.removeInvalidHexahedronsInBox(polycubeTetMesh.m_min_distance, dir);
    polycubeTetMesh.FixMesh();
    polycubeTetMesh.WriteHexahedralmesh(strOutputPolyHexFileName.c_str());

    return 0;
}

int ParseArg(int argc, char* argv[],
        std::string& strInputPolyTetFileName,
        std::string& strOutputPolyHexFileName,
        std::string& strOutputAdjustedPolyTetMeshFilename,
        float& delta, float& x, float& y, float& z)
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
            ("input-polycube-tet-file,i", po::value<std::string>(&strInputPolyTetFileName)->default_value("polycube.tet.vtk"),
                "input-polycube-tet-file")
            // ########## output #########
            // -o
            ("output-polycube-hex-file,o", po::value<std::string>(&strOutputPolyHexFileName)->default_value("polycube.hex.vtk"),
                "output-polycube-hex-file")
            // -a
            ("adjusted-tet-file,a", po::value<std::string>(&strOutputAdjustedPolyTetMeshFilename)->default_value("adjusted.polycube.tet.vtk"),
                "adjusted-tet-file")
            // ########## parameters #########
            // -s
            ("cube-size,s", po::value<float>(&delta)->default_value(0.02), "cube size, 0.001f <= size <= 0.1f")
            // -w
            ("offset,w", po::value<float>(&__plus)->default_value(0.000), "0.01s <= offset <= s")
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

        if (vm.count("input-polycube-tet-file"))
            std::cout << "InputPolyTriFileName = " << strInputPolyTetFileName << std::endl;
        if (vm.count("output-polycube-hex-file"))
            std::cout << "OutputPolyHexFileName = " << strOutputPolyHexFileName << std::endl;
        if (vm.count("adjusted-tet-file"))
            std::cout << "strOutputAdjustedPolyTetMeshFilename = " << strOutputAdjustedPolyTetMeshFilename << std::endl;
        if (vm.count("thread"))
            std::cout << "threadnum = " << THREAD_NUM << std::endl;
        if (vm.count("cube-size"))
            std::cout << "hexsize = " << delta << std::endl;
        if (vm.count("offset"))
            std::cout << "offset = " << __plus << std::endl;
        if (vm.count("x-aniso"))
            std::cout << "x-aniso = " << x << std::endl;
        if (vm.count("y-aniso"))
            std::cout << "y-aniso = " << y << std::endl;
        if (vm.count("z-aniso"))
            std::cout << "z-aniso = " << z << std::endl;

#else
        strInputPolyTetFileName = "polycube.tet.vtk";
        strOutputPolyHexFileName = "polycube.hex.vtk";
        strOutputAdjustedPolyTetMeshFilename = "adjusted.polycube.tet.vtk";
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

Eigen::MatrixXd V, U;
Eigen::MatrixXi F;
Eigen::VectorXi S, b;
Eigen::RowVector3d mid;
igl::ARAPData arap_data;

void update_vertices(Mesh& mesh)
{
    using namespace Eigen;
    using namespace std;
    MatrixXd bc(b.size(), V.cols());
    for (int i = 0; i < b.size(); i++)
    {
        const Vertex& v = mesh.V.at((unsigned long)b(i));
        bc(i, 0) = v.x;
        bc(i, 1) = v.y;
        bc(i, 2) = v.z;
    }
    igl::arap_solve(bc, arap_data, U);
    for (unsigned long i = 0; i < U.rows(); i++)
    {
        Vertex& v = mesh.V.at(i);
        v.x = U(i,0);
        v.y = U(i,1);
        v.z = U(i,2);
    }
}

void set_handle_points(const Mesh& mesh, const std::vector<unsigned long>& selectedPointIds)
{
    using namespace Eigen;
    using namespace std;
    //igl::readOFF("tri.off", V, F);
    V = Eigen::MatrixXd::Zero(mesh.V.size(), 3);
    F = Eigen::MatrixXi::Zero(mesh.C.size(), mesh.C.at(0).size());
    double p[3];
    for (int i = 0; i < V.rows(); i++)
    {
        const Vertex& v = mesh.V.at(i);
        V(i,0) = v.x;
        V(i,1) = v.y;
        V(i,2) = v.z;
    }
    for (int i = 0; i < F.rows(); i++)
        for (int j = 0; j < F.cols(); j++)
            F(i, j) = mesh.C.at(i).at(j);
    U = V;
    S = Eigen::VectorXi::Zero(mesh.V.size());
    for(int i = 0; i < V.rows(); i++) S(i) = -1;
    for(int i = 0; i < selectedPointIds.size(); i++) S(selectedPointIds[i]) = 1;
    // vertices in selection
    igl::colon<int>(0, V.rows() - 1, b);  // b is 0, 1, 2, ..., V.rows() - 1

    b.conservativeResize(stable_partition(b.data(), b.data() + b.size(), [](int i)->bool
    {   return S(i)>=0;}) - b.data());  // get the indices of selected points

    mid = 0.5 * (V.colwise().maxCoeff() + V.colwise().minCoeff());
    // Precomputation
    arap_data.max_iter = 1000;
    igl::arap_precomputation(V, F, V.cols(), b, arap_data);
}
