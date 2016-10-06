/*
 * arap_tet.cpp
 *
 *  Created on: Aug 11, 2016
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
#include <Eigen/Cholesky>
#include <Eigen/LU>
#include <Eigen/SVD>
#include <Eigen/SparseCore>

void ReadSelectedTetIndices(const char* filename, std::vector<unsigned long>& selectedTetIndices)
{
    std::ifstream ifs(filename);
    unsigned long index = 0;
    while (ifs >> index)
        selectedTetIndices.push_back(index);
}

void ReadSelectedPointIndices(const char* filename, std::vector<unsigned long>& selectedPointIndices)
{
    std::ifstream ifs(filename);
    unsigned long index = 0;
    while (ifs >> index)
        selectedPointIndices.push_back(index);

    std::sort(selectedPointIndices.begin(), selectedPointIndices.end());
    std::vector<unsigned long>::iterator iter = std::unique(selectedPointIndices.begin(), selectedPointIndices.end());
    selectedPointIndices.resize(std::distance(selectedPointIndices.begin(), iter));
}

void GetSelectedPointIds(const Mesh& mesh, const std::vector<unsigned long>& selectedTetIndices, std::vector<unsigned long>& selectedPointIds)
{
    for (size_t i = 0; i < selectedTetIndices.size(); i++)
    {
        const unsigned long tetIndex = selectedTetIndices.at(i);
        const Cell& tet = mesh.C.at(tetIndex);
        for (size_t j = 0; j < 4; j++)
            selectedPointIds.push_back(tet.at(j));
    }

    std::sort(selectedPointIds.begin(), selectedPointIds.end());
    std::vector<unsigned long>::iterator iter = std::unique(selectedPointIds.begin(), selectedPointIds.end());
    selectedPointIds.resize(std::distance(selectedPointIds.begin(), iter));
}

void ModifySelectedPointIdsLocation(Mesh& mesh, const std::vector<unsigned long>& selectedPointIds)
{
    for (size_t i = 0; i < selectedPointIds.size(); i++)
    {
        const unsigned long pointIndex = selectedPointIds.at(i);
        Vertex& v = mesh.V.at(pointIndex);
        v.x += -0.0348762;
        v.y += -0.0303857;
        v.z += 0.071173;
    }
}

Eigen::MatrixXd V, U;
Eigen::MatrixXi F;
Eigen::VectorXi S, b;
Eigen::RowVector3d mid;
igl::ARAPData arap_data;

//bool pre_draw(igl::viewer::Viewer & viewer)
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
    V = Eigen::MatrixXd::Zero(mesh.V.size(), 3);
    F = Eigen::MatrixXi::Zero(mesh.C.size(), mesh.C.at(0).size());

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

    // Centroid
    mid = 0.5 * (V.colwise().maxCoeff() + V.colwise().minCoeff());

    // Precomputation
    arap_data.max_iter = 100;
    igl::arap_precomputation(V, F, V.cols(), b, arap_data);
}


int main(int argc, char* argv[])
{
    if (argc != 3)
    {
        std::cout << "Usage: arap_tet <input.tet.vtk> <input.tet.indices>" << std::endl;
        return -1;
    }

    MeshFileReader reader(argv[1]);
    Mesh mesh(reader.GetMesh());
    //std::vector<unsigned long> selectedTetIndices;
    //ReadSelectedTetIndices(argv[2], selectedTetIndices);
    std::vector<unsigned long> selectedPointIds;
    //GetSelectedPointIds(mesh, selectedTetIndices, selectedPointIds);
    ReadSelectedPointIndices(argv[2], selectedPointIds);
    set_handle_points(mesh, selectedPointIds);
    ModifySelectedPointIdsLocation(mesh, selectedPointIds);
    update_vertices(mesh);
    ///////////////////////////////////////////////
    MeshFileWriter meshWriter(mesh, "arap.tet.vtk");
    meshWriter.WriteFile();

    return 0;
}
