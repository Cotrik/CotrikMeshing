/*
 * magnify_segments.cpp
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

void Parametrize(const Mesh& mesh, const std::vector<unsigned long>& pointIds, const glm::vec3 cube[8],
		std::vector<int>& tetIds, std::vector<glm::vec3>& w);
void MapToCube(Mesh& mesh, const std::vector<unsigned long>& pointIds, const glm::vec3 cube[8],
		const std::vector<int>& tetId, const std::vector<glm::vec3>& w);
void update_vertices(Mesh& mesh);
void set_handle_points(const Mesh& mesh, const std::vector<unsigned long>& selectedPointIds);
int main(int argc, char* argv[])
{
	if (argc != 3)
	{
		std::cout << "Usage: magnify_segments input.vtk output.vtk" << std::endl;
		return -1;
	}
	MeshFileReader vtkReader(argv[1]);
	vtkReader.GetScalarFields();
	Mesh triMesh(vtkReader.GetMesh());
	/////////////////////////////////////////
	std::cout << "Cell Ids of segment0 : " << std::endl;
	std::vector<unsigned long> cell_segment0;
	for (size_t i = 0; i < triMesh.scalarFields.at(0).size(); i++)
	{
		if (triMesh.scalarFields.at(0).at(i) == 0)
		{
			cell_segment0.push_back(i);
			std::cout << i << " ";
		}
	}

	std::vector<unsigned long> point_segment0;
	for (size_t i = 0; i < cell_segment0.size(); i++)
	{
		const Cell& cell = triMesh.C.at(cell_segment0.at(i));
		for (size_t j = 0; j < cell.size(); j++)
			point_segment0.push_back(cell.at(j));
	}
	std::cout << std::endl;

    std::sort(point_segment0.begin(), point_segment0.end());
    std::vector<unsigned long>::iterator iter = std::unique(point_segment0.begin(), point_segment0.end());
    point_segment0.resize(std::distance(point_segment0.begin(), iter));

	std::cout << point_segment0.size() << " Vertex Ids of cut of segment0 : " << std::endl;
    triMesh.GetVI_VI();
	std::vector<unsigned long> cut_segment0;
	for (size_t i = 0; i < point_segment0.size(); i++)
	{
		unsigned long pointId = point_segment0.at(i);
		if (triMesh.IsCutPoint(pointId, point_segment0))
		{
			std::cout << pointId << " ";
			cut_segment0.push_back(pointId);
		}
	}
	std::cout << std::endl;
	///////////////////////////////////////////////
	// fit a plance
	typedef double                      FT;
	typedef CGAL::Simple_cartesian<FT>  K;
	typedef K::Line_3                   Line;
	typedef K::Plane_3                  Plane;
	typedef K::Point_3                  Point;
	typedef K::Point_2                  Point_2;
	typedef K::Triangle_3               Triangle;

	std::vector<Point> points;
	for (size_t i = 0; i < cut_segment0.size(); i++)
	{
		const Vertex& v = triMesh.V.at(cut_segment0.at(i));
		Point p(v.x, v.y, v.z);
		points.push_back(p);
	}
	Plane plane;
	linear_least_squares_fitting_3(points.begin(),points.end(),plane,CGAL::Dimension_tag<0>());
	std::cout << plane << std::endl;
	double a = plane.a();
	double b = plane.b();
	double c = plane.c();
	double d = plane.d();

	double sum = 0.0;
	double max_dis = 0;
	double max_dis_id = point_segment0.at(0);

	for (size_t i = 0; i < point_segment0.size(); i++)
	{
		const unsigned long pointId = point_segment0.at(i);
		const Vertex& v = triMesh.V.at(pointId);
		double distance = a*v.x + b*v.y + c*v.z + d;
		sum += distance;
		if (fabs(distance) > max_dis)
		{
			max_dis = fabs(distance);
			max_dis_id = pointId;
		}

//		const double t = distance; //distance * abc;
//		projectPoints.at(i) = Point(v.x - a*t, v.y - b*t, v.z - c*t);
	}
	std::cout << "max_dis up = " << max_dis << std::endl;
	std::cout << "max_dis_id = " << max_dis_id << std::endl;
	if (sum < 0)
	{
		a = -a; b = -b; c = -c; d = -d;
	}
	for (iter = point_segment0.begin(); iter != point_segment0.end();)
	{
		const Vertex& v = triMesh.V.at(*iter);
		double distance = a * v.x + b * v.y + c * v.z + d;
		if (distance < 0)
		{
			iter = point_segment0.erase(iter);
		}
		else
		{
			++iter;
		}
	}
	std::vector<Point> projectPoints(point_segment0.size());
	const double abc = 1.0/(a*a + b*b + c*c);
	for (size_t i = 0; i < point_segment0.size(); i++)
	{
		const unsigned long pointId = point_segment0.at(i);
		const Vertex& v = triMesh.V.at(pointId);
		double distance = a * v.x + b * v.y + c * v.z + d;
		const double t = distance; //distance * abc;
		projectPoints.at(i) = Point(v.x - a*t, v.y - b*t, v.z - c*t);
	}
	///////////////////////////////////////////////
	// Bounding Box axis aligned box
//	K::Iso_cuboid_3 boundingBox = CGAL::bounding_box(projectPoints.begin(), projectPoints.end());
//	std::cout << "boundingBox: " << boundingBox << std::endl;
	// get the maxdistance in the projectPoints
	unsigned long id1 = point_segment0.at(0);
	unsigned long id2 = point_segment0.at(1);
	unsigned long id_1 = point_segment0.at(0);
	unsigned long id_2 = point_segment0.at(1);
	glm::vec3 pt1(projectPoints[0].x(),projectPoints[0].y(),projectPoints[0].z());
	glm::vec3 pt2(projectPoints[1].x(),projectPoints[1].y(),projectPoints[1].z());
	double max_dis2 = glm::length(pt1 - pt2);
	for (size_t i = 0; i < projectPoints.size(); i++)
	{
		const Point& p1 = projectPoints.at(i);
		for (size_t j = 0; j < projectPoints.size(); j++)
		{
			const Point& p2 = projectPoints.at(j);
			if (i != j)
			{
				glm::vec3 p1(projectPoints[i].x(),projectPoints[i].y(),projectPoints[i].z());
				glm::vec3 p2(projectPoints[j].x(),projectPoints[j].y(),projectPoints[j].z());
				double distance = glm::length(p1 - p2);
				if (distance > max_dis2)
				{
					max_dis2 = distance;
					id1 = point_segment0.at(i);
					id2 = point_segment0.at(j);
					id_1 = i;
					id_2 = j;
				}
			}
		}
	}
	max_dis *= 1.1;
	max_dis2 *= 1.1;
	std::cout << "max_dis2 = " << max_dis2 << std::endl;
	std::cout << "id1 = " << id1 << std::endl;
	std::cout << "id2 = " << id2 << std::endl;
//	const Vertex& v1 = triMesh.V.at(id1);
//	const Vertex& v2 = triMesh.V.at(id2);
	const glm::vec3 v1(projectPoints[id_1].x(), projectPoints[id_1].y(), projectPoints[id_1].z());
	const glm::vec3 v2(projectPoints[id_2].x(), projectPoints[id_2].y(), projectPoints[id_2].z());
	glm::vec3 dir1(a, b, c);
	glm::vec3 dir2(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z);
	double dir2_len = glm::length(dir2);
	dir2.x /= dir2_len; dir2.y /= dir2_len; dir2.z /= dir2_len;
	glm::vec3 dir3 = glm::cross(dir1, dir2);
	std::cout << "max_dis3 = " << glm::length(dir3) << std::endl;
	///////////////////
	std::cout << "dir1 = (" << dir1.x << ", " << dir1.y << ", " << dir1.z << ")"  << std::endl;
	std::cout << "dir2 = (" << dir2.x << ", " << dir2.y << ", " << dir2.z << ")"  << std::endl;
	std::cout << "dir3 = (" << dir3.x << ", " << dir3.y << ", " << dir3.z << ")"  << std::endl;
	glm::vec3 d1 = glm::cross(dir2, dir3);
	glm::vec3 d2 = glm::cross(dir3, dir1);
	glm::vec3 d3 = glm::cross(dir1, dir2);
	std::cout << "d1 = (" << d1.x << ", " << d1.y << ", " << d1.z << ")"  << std::endl;
	std::cout << "d2 = (" << d2.x << ", " << d2.y << ", " << d2.z << ")"  << std::endl;
	std::cout << "d3 = (" << d3.x << ", " << d3.y << ", " << d3.z << ")"  << std::endl;
	///////////////////
	const Vertex v_mid = Vertex(0.5*(v1.x + v2.x), 0.5*(v1.y + v2.y), 0.5*(v1.z + v2.z));
	const Vertex v_up(v_mid.x + a*max_dis, v_mid.y + b*max_dis, v_mid.z + c*max_dis);
	const double half_dis = max_dis2 * 0.5;
	std::cout << "half_dis = " << half_dis << std::endl;
	glm::vec3 cube[8] =
	{
			glm::vec3(v1.x + dir3.x*half_dis, v1.y + dir3.y*half_dis, v1.z + dir3.z*half_dis),
			glm::vec3(v2.x + dir3.x*half_dis, v2.y + dir3.y*half_dis, v2.z + dir3.z*half_dis),
			glm::vec3(v2.x - dir3.x*half_dis, v2.y - dir3.y*half_dis, v2.z - dir3.z*half_dis),
			glm::vec3(v1.x - dir3.x*half_dis, v1.y - dir3.y*half_dis, v1.z - dir3.z*half_dis),
	};
	for (int i = 4; i < 8; i++)
	{
		std::cout << "dis = " << glm::length(cube[i-4] - v1) << std::endl;
		cube[i].x = cube[i - 4].x + dir1.x*max_dis;
		cube[i].y = cube[i - 4].y + dir1.y*max_dis;
		cube[i].z = cube[i - 4].z + dir1.z*max_dis;
	}
	MeshFileReader cubeReader("cube.vtk");
	Mesh cubeMesh(cubeReader.GetMesh());
	for (int i = 0; i < 8; i++)
	{
		cubeMesh.V.at(i).x = cube[i].x;
		cubeMesh.V.at(i).y = cube[i].y;
		cubeMesh.V.at(i).z = cube[i].z;
	}
//	MeshFileWriter cubeWriter(cubeMesh, "/home/cotrik/Downloads/software/CGAL-4.7/build/examples/Surface_mesh_skeletonization/data/cube2.vtk", MESH_TYPE_HEXAHEDRON_VTK);
//	cubeWriter.WriteMeshFile("/home/cotrik/Downloads/software/CGAL-4.7/build/examples/Surface_mesh_skeletonization/data/cube2.vtk", MESH_TYPE_HEXAHEDRON_VTK);
	///////////////////////////////////////////////
	// magnify cube *2
	glm::vec3 v1m(v_mid.x + dir2.x*max_dis2, v_mid.y + dir2.y*max_dis2, v_mid.z + dir2.z*max_dis2);
	glm::vec3 v2m(v_mid.x - dir2.x*max_dis2, v_mid.y - dir2.y*max_dis2, v_mid.z - dir2.z*max_dis2);
	glm::vec3 cube2[8] =
	{
			glm::vec3(v1m.x + dir3.x*max_dis2, v1m.y + dir3.y*max_dis2, v1m.z + dir3.z*max_dis2),
			glm::vec3(v2m.x + dir3.x*max_dis2, v2m.y + dir3.y*max_dis2, v2m.z + dir3.z*max_dis2),
			glm::vec3(v2m.x - dir3.x*max_dis2, v2m.y - dir3.y*max_dis2, v2m.z - dir3.z*max_dis2),
			glm::vec3(v1m.x - dir3.x*max_dis2, v1m.y - dir3.y*max_dis2, v1m.z - dir3.z*max_dis2),
	};
	for (int i = 4; i < 8; i++)
	{
		cube2[i].x = cube2[i - 4].x + dir1.x*max_dis*2;
		cube2[i].y = cube2[i - 4].y + dir1.y*max_dis*2;
		cube2[i].z = cube2[i - 4].z + dir1.z*max_dis*2;
	}
	///////////////////////////////////////////////
	// parametrize vertex
	std::vector<int> tetIds;
	std::vector<glm::vec3> w;
	set_handle_points(triMesh, point_segment0);
	Parametrize(triMesh, point_segment0, cube, tetIds, w);
	MapToCube(triMesh, point_segment0, cube2, tetIds, w);
	update_vertices(triMesh);
	///////////////////////////////////////////////
	MeshFileWriter meshWriter(triMesh, argv[2]);
	meshWriter.WriteFile();

	std::cout << "Job Finished!" << std::endl;
	return 0;
}

void Parametrize(const Mesh& mesh, const std::vector<unsigned long>& pointIds, const glm::vec3 cube[8],
		std::vector<int>& tetIds, std::vector<glm::vec3>& w)
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

			if (rambda.x > -1e-3 && rambda.y > -1e-3 && rambda.z>-1e-3 && (rambda.x + rambda.y + rambda.z) < 1.001){
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

void MapToCube(Mesh& mesh, const std::vector<unsigned long>& pointIds, const glm::vec3 cube[8],
		const std::vector<int>& tetId, const std::vector<glm::vec3>& w)
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


typedef std::vector<Eigen::Quaterniond, Eigen::aligned_allocator<Eigen::Quaterniond> > RotationList;

const Eigen::RowVector3d sea_green(70. / 255., 252. / 255., 167. / 255.);
Eigen::MatrixXd V, U;
Eigen::MatrixXi F;
Eigen::VectorXi S, b;
Eigen::RowVector3d mid;
double anim_t = 0.0;
double anim_t_dir = 0.03;
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
    //igl::readDMAT("decimated-knight-selection.dmat", S);
    S = Eigen::VectorXi::Zero(mesh.V.size());
    for(int i = 0; i < V.rows(); i++) S(i) = -1;
    for(int i = 0; i < selectedPointIds.size(); i++) S(selectedPointIds[i]) = 1;
    // vertices in selection
    igl::colon<int>(0, V.rows() - 1, b);  // b is 0, 1, 2, ..., V.rows() - 1
//    int b1[502];
//    int s[502];
//    for (int i = 0; i < 502; i++)
//    {
//        b1[i] = b(i);
//        s[i] = S(i);
//    }
    b.conservativeResize(stable_partition(b.data(), b.data() + b.size(), [](int i)->bool
    {   return S(i)>=0;}) - b.data());  // get the indices of selected points

//    int b2[73];
//    for (int i = 0; i < 73; i++) b2[i] = b(i);
    // Centroid
    mid = 0.5 * (V.colwise().maxCoeff() + V.colwise().minCoeff());
    // Precomputation
    arap_data.max_iter = 100;
    igl::arap_precomputation(V, F, V.cols(), b, arap_data);

    // Set color based on selection
    MatrixXd C(F.rows(), 3);
    RowVector3d purple(80.0 / 255.0, 64.0 / 255.0, 255.0 / 255.0);
    RowVector3d gold(255.0 / 255.0, 228.0 / 255.0, 58.0 / 255.0);
    for (int f = 0; f < F.rows(); f++)
    {
        if (S(F(f, 0)) >= 0 && S(F(f, 1)) >= 0 && S(F(f, 2)) >= 0)
        {
            C.row(f) = purple;
        }
        else
        {
            C.row(f) = gold;
        }
    }
}
