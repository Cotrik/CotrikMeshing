/*
 * magnify_multi_segments_scales.cpp
 *
 *  Created on: May 18, 2016
 *      Author: cotrik
 */

#include "include.h"
#include <CGAL/Simple_cartesian.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/bounding_box.h>
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

#include <eigen3/Eigen/Geometry>
#include <eigen3/Eigen/StdVector>
#include "eigen3/Eigen/Core"
#include "eigen3/Eigen/Eigen"
#include "eigen3/Eigen/Dense"
#include "eigen3/Eigen/Sparse"
//#include "eigen3/Eigen/MPRealSupport"
#include "eigen3/Eigen/Cholesky"
#include "eigen3/Eigen/LU"
#include "eigen3/Eigen/SVD"
#include "eigen3/Eigen/SparseCore"

    ///////////////////////////////////////////////
    // fit a plance
    typedef double FT;
    typedef CGAL::Simple_cartesian<FT> K;
    typedef K::Line_3 LineC;
    typedef K::Plane_3 PlaneC;
    typedef K::Point_3 PointC;
    typedef K::Point_2 Point_2;
    typedef K::Triangle_3 TriangleC;

void Parametrize(const Mesh& mesh, const std::vector<unsigned long>& pointIds, const glm::vec3 cube[8], std::vector<int>& tetIds, std::vector<glm::vec3>& w);
void MapToCube(Mesh& mesh, const std::vector<unsigned long>& pointIds, const glm::vec3 cube[8], const std::vector<int>& tetId, const std::vector<glm::vec3>& w);
void update_vertices(Mesh& mesh);
void set_handle_points(const Mesh& mesh, const std::vector<unsigned long>& selectedPointIds);
void read_segmentids(const char* filename, std::vector<unsigned long>& segmentIds);
void read_scales(const char* filename, std::vector<double>& scales);
void get_cell_segment_of_ROI(const Mesh& triMesh, const unsigned long segmentId, std::vector<unsigned long>& cell_segment0, double& avg_seg_medial_size)
{
    double sum_seg_medial_size = 0.0;
    for (size_t i = 0; i < triMesh.scalarFields.at(0).size(); i++)
    {
        if (triMesh.scalarFields.at(0).at(i) == segmentId)
        {
            cell_segment0.push_back(i);
            sum_seg_medial_size += triMesh.scalarFields.at(1).at(i);
        }
    }
    avg_seg_medial_size = sum_seg_medial_size / cell_segment0.size();
}

void get_point_segment_of_ROI(const Mesh& triMesh, const std::vector<unsigned long>& cell_segment, std::vector<unsigned long>& point_segment)
{
    for (size_t i = 0; i < cell_segment.size(); i++)
    {
        const Cell& cell = triMesh.C.at(cell_segment.at(i));
        for (size_t j = 0; j < cell.size(); j++)
            point_segment.push_back(cell.at(j));
    }

    std::sort(point_segment.begin(), point_segment.end());
    std::vector<unsigned long>::iterator iter = std::unique(point_segment.begin(), point_segment.end());
    point_segment.resize(std::distance(point_segment.begin(), iter));
}

void get_cutplane_point_of_ROI(Mesh& triMesh, const std::vector<unsigned long>& point_segment, std::vector<unsigned long>& cut_segment)
{
    for (size_t i = 0; i < point_segment.size(); i++)
    {
        unsigned long pointId = point_segment.at(i);
        if (triMesh.IsCutPoint(pointId, point_segment))
        {
            std::cout << pointId << " ";
            cut_segment.push_back(pointId);
        }
    }
    std::cout << std::endl;

    //////////////////////////////
    std::vector<std::vector<unsigned long> > groupIds;
    int groups = triMesh.GroupCutPoints(cut_segment, groupIds);
    std::cout << "groups = " << groups << std::endl;
    cut_segment = groupIds.at(0);
    if (groupIds.size() > 1)
        cut_segment = groupIds.at(1);
    for (size_t i = 0; i < cut_segment.size(); i++)
    {
        std::cout << cut_segment.at(i) << " ";
    }
    std::cout << std::endl;
    if (groupIds.size() > 1)
        for (size_t i = 0; i < groupIds.at(1).size(); i++)
        {
            std::cout << groupIds.at(1).at(i) << " ";
        }
    std::cout << std::endl;
    //////////////////////////////
}

void get_cut_plane(const Mesh& triMesh, const std::vector<unsigned long>& cut_segment, PlaneC& plane)
{
    std::vector<PointC> points;
    for (size_t i = 0; i < cut_segment.size(); i++)
    {
        const Vertex& v = triMesh.V.at(cut_segment.at(i));
        PointC p(v.x, v.y, v.z);
        points.push_back(p);
    }
    //Plane plane;
    linear_least_squares_fitting_3(points.begin(), points.end(), plane, CGAL::Dimension_tag<0>());
    std::cout << plane << std::endl;
}

double get_max_distance_from_cut_plane_for_points_in_ROI(const Mesh& triMesh, const std::vector<unsigned long>& point_segment,
        double& a, double& b, double& c, double& d // cutplane ax + by + cz + d = 0;
        )
{
    double max_dis = 0.0;
    double max_dis_id = point_segment.at(0);
    double sum = 0.0;
    for (size_t i = 0; i < point_segment.size(); i++)
    {
        const unsigned long pointId = point_segment.at(i);
        const Vertex& v = triMesh.V.at(pointId);
        double distance = a * v.x + b * v.y + c * v.z + d;
        sum += distance;
        if (fabs(distance) > max_dis)
        {
            max_dis = fabs(distance);
            max_dis_id = pointId;
        }
//      const double t = distance; //distance * abc;
//      projectPoints.at(i) = Point(v.x - a*t, v.y - b*t, v.z - c*t);
    }

    // make sure plane's normal point to ROI.
    if (sum < 0)
    {
        a = -a;
        b = -b;
        c = -c;
        d = -d;
    }

    return max_dis;
}

void erase_invalid_points_beyond_ROI(const Mesh& triMesh, std::vector<unsigned long>& point_segment,
        const double a, const double b, const double c, const double d // cutplane ax + by + cz + d = 0;
        )
{
    for (std::vector<unsigned long>::iterator iter = point_segment.begin(); iter != point_segment.end();)
    {
        const Vertex& v = triMesh.V.at(*iter);
        double distance = a * v.x + b * v.y + c * v.z + d;
        if (distance < 0)
        {
            iter = point_segment.erase(iter);
        }
        else
        {
            ++iter;
        }
    }
}

void get_project_points_location_for_points_of_ROI(const Mesh& triMesh,
        const std::vector<unsigned long>& point_segment,
        std::vector<PointC> projectPoints,
        const double a, const double b, const double c, const double d // cutplane ax + by + cz + d = 0;
        )
{
    projectPoints.resize(point_segment.size());
    const double abc = 1.0 / (a * a + b * b + c * c);
    for (size_t i = 0; i < point_segment.size(); i++)
    {
        const unsigned long pointId = point_segment.at(i);
        const Vertex& v = triMesh.V.at(pointId);
        double distance = a * v.x + b * v.y + c * v.z + d;
        const double t = distance; //distance * abc;
        projectPoints.at(i) = PointC(v.x - a * t, v.y - b * t, v.z - c * t);
    }
}

double get_max_distance_for_projectPoints(const std::vector<unsigned long>& point_segment, const std::vector<PointC>& projectPoints,
        unsigned long& id_1, unsigned long& id_2) // id_1 and id_2 are the two indices in projectPoints that will have the long distance
{
    // get the maxdistance in the projectPoints
    unsigned long id1 = point_segment.at(0);
    unsigned long id2 = point_segment.at(1);
    id_1 = point_segment.at(0);
    id_2 = point_segment.at(1);
    glm::vec3 pt1(projectPoints[0].x(), projectPoints[0].y(), projectPoints[0].z());
    glm::vec3 pt2(projectPoints[1].x(), projectPoints[1].y(), projectPoints[1].z());
    double max_dis2 = glm::length(pt1 - pt2);
    for (size_t i = 0; i < projectPoints.size(); i++)
    {
        const PointC& p1 = projectPoints.at(i);
        for (size_t j = 0; j < projectPoints.size(); j++)
        {
            const PointC& p2 = projectPoints.at(j);
            if (i != j)
            {
                glm::vec3 p1(projectPoints[i].x(), projectPoints[i].y(), projectPoints[i].z());
                glm::vec3 p2(projectPoints[j].x(), projectPoints[j].y(), projectPoints[j].z());
                double distance = glm::length(p1 - p2);
                if (distance > max_dis2)
                {
                    max_dis2 = distance;
                    id1 = point_segment.at(i);
                    id2 = point_segment.at(j);
                    id_1 = i;
                    id_2 = j;
                }
            }
        }
    }

    return max_dis2;
}

double get_cubes(const std::vector<PointC>& projectPoints, const unsigned long id_1, const unsigned long id_2, const double max_dis, const double max_dis2,
        const double a, const double b, const double c, const double d, // cutplane ax + by + cz + d = 0;
        glm::vec3 cube[8], glm::vec3 cube2[8]
        )
{
    const glm::vec3 v1(projectPoints[id_1].x(), projectPoints[id_1].y(), projectPoints[id_1].z());
    const glm::vec3 v2(projectPoints[id_2].x(), projectPoints[id_2].y(), projectPoints[id_2].z());
    glm::vec3 dir1(a, b, c);
    glm::vec3 dir2(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z);
    double dir2_len = glm::length(dir2);
    dir2.x /= dir2_len;
    dir2.y /= dir2_len;
    dir2.z /= dir2_len;
    glm::vec3 dir3 = glm::cross(dir1, dir2);
    ///////////////////
    const Vertex v_mid = Vertex(0.5 * (v1.x + v2.x), 0.5 * (v1.y + v2.y), 0.5 * (v1.z + v2.z));
    const Vertex v_up(v_mid.x + a * max_dis, v_mid.y + b * max_dis, v_mid.z + c * max_dis);
    const double half_dis = max_dis2 * 0.5;
    //std::cout << "half_dis = " << half_dis << std::endl;
    //glm::vec3 cube[8] =
    {
        cube[0] = glm::vec3(v1.x + dir3.x * half_dis, v1.y + dir3.y * half_dis, v1.z + dir3.z * half_dis);
        cube[1] = glm::vec3(v2.x + dir3.x * half_dis, v2.y + dir3.y * half_dis, v2.z + dir3.z * half_dis);
        cube[2] = glm::vec3(v2.x - dir3.x * half_dis, v2.y - dir3.y * half_dis, v2.z - dir3.z * half_dis);
        cube[3] = glm::vec3(v1.x - dir3.x * half_dis, v1.y - dir3.y * half_dis, v1.z - dir3.z * half_dis);
    }
    for (int i = 4; i < 8; i++)
    {
        //std::cout << "dis = " << glm::length(cube[i-4] - v1) << std::endl;
        cube[i].x = cube[i - 4].x + dir1.x * max_dis;
        cube[i].y = cube[i - 4].y + dir1.y * max_dis;
        cube[i].z = cube[i - 4].z + dir1.z * max_dis;
    }

    ///////////////////////////////////////////////
    // magnify cube *2
    //double scale = 0.5* avg_medial_size / avg_seg_medial_size;
    //double scale = scales.at(seg_id);
    double scale1 = 3, scale2 = 2, scale3 = 1;
    std::cout << "scale1 = " << scale1 << std::endl;
    std::cout << "scale2 = " << scale2 << std::endl;
    std::cout << "scale3 = " << scale3 << std::endl;
//        if (scale < 1.0)
//            continue;
    double size1 = scale1 * max_dis;
    double size2 = scale2 * half_dis;
    double size3 = scale3 * half_dis;
//  if (groupIds.size() > 1)
//      scale = 0.5;
    glm::vec3 v1m(v_mid.x + dir2.x * size2, v_mid.y + dir2.y * size2, v_mid.z + dir2.z * size2);
    glm::vec3 v2m(v_mid.x - dir2.x * size2, v_mid.y - dir2.y * size2, v_mid.z - dir2.z * size2);
    //glm::vec3 cube2[8] =
    {
        cube2[0] = glm::vec3(v1m.x + dir3.x * size3, v1m.y + dir3.y * size3, v1m.z + dir3.z * size3);
        cube2[1] = glm::vec3(v2m.x + dir3.x * size3, v2m.y + dir3.y * size3, v2m.z + dir3.z * size3);
        cube2[2] = glm::vec3(v2m.x - dir3.x * size3, v2m.y - dir3.y * size3, v2m.z - dir3.z * size3);
        cube2[3] = glm::vec3(v1m.x - dir3.x * size3, v1m.y - dir3.y * size3, v1m.z - dir3.z * size3);
    };
    for (int i = 4; i < 8; i++)
    {
        cube2[i].x = cube2[i - 4].x + dir1.x * max_dis * scale1;
        cube2[i].y = cube2[i - 4].y + dir1.y * max_dis * scale1;
        cube2[i].z = cube2[i - 4].z + dir1.z * max_dis * scale1;
    }
    ///////////////////////////////////////////////
}

int main(int argc, char* argv[])
{
    if (argc != 6)
    {
        std::cout << "Usage: magnify_multi_segments_scales_volume_aniso input_tet_file input_tri_seg_file output_tet_file segmentid.txt scales.txt" << std::endl;
        return -1;
    }
    MeshFileReader tetReader(argv[1]);
    Mesh tetMesh(tetReader.GetMesh());
    tetMesh.ExtractSurface();

    MeshFileReader vtkReader(argv[2]);
    vtkReader.GetScalarFields();
    Mesh triMesh(vtkReader.GetMesh());
    triMesh.GetVI_VI();

    /////////////////////////////////////////
    std::vector<unsigned long> segmentIds;
    read_segmentids(argv[4], segmentIds);
    std::vector<double> scales;
    read_scales(argv[5], scales);
    //std::vector<unsigned long> selectedPoints;
    double sum_medial_size = 0.0;
    for (size_t i = 0; i < triMesh.scalarFields.at(0).size(); i++)
    {
        sum_medial_size += triMesh.scalarFields.at(1).at(i);
    }
    double avg_medial_size = sum_medial_size / triMesh.C.size();
    /////////////////////////////////////////
    for (int seg_id = 0; seg_id < segmentIds.size(); seg_id++)
    {
        const unsigned long segmentId = segmentIds.at(seg_id);
        std::cout << "Cell Ids of segment " << segmentId << " : " << std::endl;

        std::vector<unsigned long> cell_segment0;
        double avg_seg_medial_size = 0.0;
        get_cell_segment_of_ROI(triMesh, segmentId, cell_segment0, avg_seg_medial_size);

        std::vector<unsigned long> point_segment0;
        get_point_segment_of_ROI(triMesh, cell_segment0, cell_segment0);
        // std::cout << point_segment0.size() << " Vertex Ids of point segment : " << std::endl;

        std::vector<unsigned long> cut_segment0;
        get_cutplane_point_of_ROI(triMesh, cell_segment0, cut_segment0);

        PlaneC plane;
        get_cut_plane(triMesh, cut_segment0, plane);
        double a = plane.a();
        double b = plane.b();
        double c = plane.c();
        double d = plane.d();
        double max_dis = get_max_distance_from_cut_plane_for_points_in_ROI(triMesh, point_segment0, a, b, c, d);

        erase_invalid_points_beyond_ROI(triMesh, point_segment0, a, b, c, d);

        std::vector<PointC> projectPoints(point_segment0.size());
        get_project_points_location_for_points_of_ROI(triMesh, point_segment0, projectPoints, a, b, c, d);

        unsigned long id_1 = point_segment0.at(0);
        unsigned long id_2 = point_segment0.at(1);
        double max_dis2 = get_max_distance_for_projectPoints(point_segment0, projectPoints, id_1, id_2);

        max_dis *= 1.1;
        max_dis2 *= 1.1;
        glm::vec3 cube[8];
        glm::vec3 cube2[8];
        get_cubes(projectPoints, id_1, id_2, max_dis, max_dis2, a, b, c, d, cube, cube2);
        MeshFileReader cubeReader("cube.vtk");
        Mesh cubeMesh(cubeReader.GetMesh());
        for (int i = 0; i < 8; i++)
        {
            cubeMesh.V.at(i).x = cube[i].x;
            cubeMesh.V.at(i).y = cube[i].y;
            cubeMesh.V.at(i).z = cube[i].z;
        }
//        MeshFileWriter cubeWriter(cubeMesh, "/home/cotrik/Downloads/software/CGAL-4.7/build/examples/Surface_mesh_skeletonization/data/cube2.vtk");
//        cubeWriter.WriteFile();

        // parametrize vertex
        std::vector<int> tetIds;
        std::vector<glm::vec3> w;
        set_handle_points(triMesh, point_segment0);
        Parametrize(triMesh, point_segment0, cube, tetIds, w);
        MapToCube(triMesh, point_segment0, cube2, tetIds, w);
        update_vertices(triMesh);
    }
    ///////////////////////////////////////////////
    MeshFileWriter meshWriter(triMesh, argv[3]);
    meshWriter.WriteFile();

    std::cout << "Job Finished!" << std::endl;
    return 0;
}
//std::vector<glm::vec3> paras;
void Parametrize(const Mesh& mesh, const std::vector<unsigned long>& pointIds, const glm::vec3 cube[8], std::vector<int>& tetIds, std::vector<glm::vec3>& w)
{
    for (int i = 0; i < pointIds.size(); i++)
    {
        const Vertex& v = mesh.V.at(pointIds.at(i));
        // 1. judge v is in which tet
        bool flag = false;
        for (int k = 0; k < 5; k++)
        {
            const glm::vec3& p0 = cube[hexTet[k][0]];
            const glm::vec3& p1 = cube[hexTet[k][1]];
            const glm::vec3& p2 = cube[hexTet[k][2]];
            const glm::vec3& p3 = cube[hexTet[k][3]];

            glm::vec3 v03((p0.x - p3.x), (p0.y - p3.y), (p0.z - p3.z));
            glm::vec3 v13((p1.x - p3.x), (p1.y - p3.y), (p1.z - p3.z));
            glm::vec3 v23((p2.x - p3.x), (p2.y - p3.y), (p2.z - p3.z));

            glm::mat3x3 T(v03, v13, v23);
            glm::vec3 r(v.x, v.y, v.z);
            glm::vec3 r4(p3.x, p3.y, p3.z);
            glm::mat3x3 T_inverse = glm::inverse(T);
            glm::vec3 rambda = T_inverse * (r - r4);
            //paras.push_back(rambda);

            if (rambda.x > -1e-3 && rambda.y > -1e-3 && rambda.z > -1e-3 && (rambda.x + rambda.y + rambda.z) < 1.001)
            {
                w.push_back(rambda);
                tetIds.push_back(k);
                flag = true;
                break;
            }
        }
        if (!flag)
        {
            std::cout << pointIds.at(i) << " fail parametrization! \n";
        }
    }
}

void MapToCube(Mesh& mesh, const std::vector<unsigned long>& pointIds, const glm::vec3 cube[8], const std::vector<int>& tetId, const std::vector<glm::vec3>& w)
{
    std::cout << "MapToCube" << std::endl;
    for (int i = 0; i < pointIds.size(); i++)
    {
        Vertex& v = mesh.V.at(pointIds.at(i));
        int k = tetId.at(i);
        const glm::vec3& p0 = cube[hexTet[k][0]];
        const glm::vec3& p1 = cube[hexTet[k][1]];
        const glm::vec3& p2 = cube[hexTet[k][2]];
        const glm::vec3& p3 = cube[hexTet[k][3]];
        v.x = w[i].x * p0.x + w[i].y * p1.x + w[i].z * p2.x + (1.0 - w[i].x - w[i].y - w[i].z) * p3.x;
        v.y = w[i].x * p0.y + w[i].y * p1.y + w[i].z * p2.y + (1.0 - w[i].x - w[i].y - w[i].z) * p3.y;
        v.z = w[i].x * p0.z + w[i].y * p1.z + w[i].z * p2.z + (1.0 - w[i].x - w[i].y - w[i].z) * p3.z;
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
        const Vertex& v = mesh.V.at((unsigned long) b(i));
        bc(i, 0) = v.x;
        bc(i, 1) = v.y;
        bc(i, 2) = v.z;
    }
    igl::arap_solve(bc, arap_data, U);
    for (unsigned long i = 0; i < U.rows(); i++)
    {
        Vertex& v = mesh.V.at(i);
        v.x = U(i, 0);
        v.y = U(i, 1);
        v.z = U(i, 2);
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
        V(i, 0) = v.x;
        V(i, 1) = v.y;
        V(i, 2) = v.z;
    }
    for (int i = 0; i < F.rows(); i++)
        for (int j = 0; j < F.cols(); j++)
            F(i, j) = mesh.C.at(i).at(j);
    U = V;
    //igl::readDMAT("decimated-knight-selection.dmat", S);
    S = Eigen::VectorXi::Zero(mesh.V.size());
    for (int i = 0; i < V.rows(); i++)
        S(i) = -1;
    for (int i = 0; i < selectedPointIds.size(); i++)
        S(selectedPointIds[i]) = 1;
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

void read_segmentids(const char* filename, std::vector<unsigned long>& segmentIds)
{
    std::ifstream ifs(filename);
    unsigned long id;
    while (ifs >> id)
        segmentIds.push_back(id);
}

void read_scales(const char* filename, std::vector<double>& scales)
{
    std::ifstream ifs(filename);
    float s;
    while (ifs >> s)
        scales.push_back(s);
}
