/*
 * GeoUtil.cpp
 *
 *  Created on: Feb 6, 2015
 *      Author: cotrik
 */

#include "GeoUtil.h"
#include <iostream>
#include "Eigen/Core"
#include "Eigen/Eigen"
#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "Eigen/Cholesky"
#include "Eigen/LU"
#include "Eigen/SVD"
#include "Eigen/SparseCore"
GeoUtil::GeoUtil()
{
	// TODO Auto-generated constructor stub

}

GeoUtil::~GeoUtil()
{
	// TODO Auto-generated destructor stub
}

void GeoUtil::GetMaxCoordinateValueOfVertex(const Vertex& vertex, Vertex& vertexMax)
{
	if (vertex.x > vertexMax.x)
	{
		vertexMax.x = vertex.x;
	}
	if (vertex.y > vertexMax.y)
	{
		vertexMax.y = vertex.y;
	}
	if (vertex.z > vertexMax.z)
	{
		vertexMax.z = vertex.z;
	}
}

void GeoUtil::GetMinCoordinateValueOfVertex(const Vertex& vertex, Vertex& vertexMin)
{
	if (vertex.x < vertexMin.x)
	{
		vertexMin.x = vertex.x;
	}
	if (vertex.y < vertexMin.y)
	{
		vertexMin.y = vertex.y;
	}
	if (vertex.z < vertexMin.z)
	{
		vertexMin.z = vertex.z;
	}
}

void GeoUtil::GetMaxMinCoordinateValueOfVertex(const Vertex& vertex, Vertex& vertexMax, Vertex& vertexMin)
{
	GetMaxCoordinateValueOfVertex(vertex, vertexMax);
	GetMinCoordinateValueOfVertex(vertex, vertexMin);
}

void GeoUtil::GetMaxMinCoordinateValueOfTriangle(const std::vector<Vertex*>& m_vecVertex, const Triangle& triangle, Vertex& vertexMax, Vertex& vertexMin)
{
	GetMaxMinCoordinateValueOfVertex(*m_vecVertex[triangle.p0], vertexMax, vertexMin);
	GetMaxMinCoordinateValueOfVertex(*m_vecVertex[triangle.p1], vertexMax, vertexMin);
	GetMaxMinCoordinateValueOfVertex(*m_vecVertex[triangle.p2], vertexMax, vertexMin);
}

bool GeoUtil::IsPointInTriangle(const glm::vec3& P, const glm::vec3& A, const glm::vec3& B, const glm::vec3& C, const int i, const int j)
{
	// Compute vectors
	const glm::vec3 v0 = C - A;
	const glm::vec3 v1 = B - A;
	const glm::vec3 v2 = P - A;

	// Compute dot products
	const double dot00 = glm::dot(v0, v0);
	const double dot01 = glm::dot(v0, v1);
	const double dot02 = glm::dot(v0, v2);
	const double dot11 = glm::dot(v1, v1);
	const double dot12 = glm::dot(v1, v2);

	// Compute barycentric coordinates
	const double invDenom = 1.0 / (dot00 * dot11 - dot01 * dot01);
	const double u = (dot11 * dot02 - dot01 * dot12) * invDenom;
	const double v = (dot00 * dot12 - dot01 * dot02) * invDenom;
	if (i == 327 || i == 10515 || i == 14266)
	{
		if(fabs(u) + fabs(v) < 1.3)
		std::cout << "i = " << i << " j = " << j << " u = " << u << ", v = " << v  << std::endl;
	}
	// Check if point is in triangle
	return (u >= 0) && (v >= 0) && (u + v < 1.0);
}

bool GeoUtil::IsPointCanProjectToTriangle(const glm::vec3& orig, const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2, glm::vec3& p)
{
	glm::vec3 n = glm::cross(v1 - v0, v2 - v0);
    double d = -glm::dot(n, v0);
    double length = 1.0 / n.length();
	n *= length;
	d *= length;

	double a = n.x;
	double b = n.y;
	double c = n.z;

	double t = (glm::dot(orig, n) + d)/(a*a + b*b + c*c);

	//glm::vec3 p(orig.x - a*t, orig.y - b*t, orig.z - c*t);
	p.x = orig.x - a*t;
	p.y = orig.y - a*t;
	p.z = orig.z - a*t;

	return IsPointInTriangle(p, v0, v1, v2);
}

bool GeoUtil::IsPointInTetrahedron(const Vertex& p, const Vertex& p0, const Vertex& p1, const Vertex& p2, const Vertex& p3)
{
	glm::vec3 v03((p0.x - p3.x), (p0.y - p3.y), (p0.z - p3.z));
	glm::vec3 v13((p1.x - p3.x), (p1.y - p3.y), (p1.z - p3.z));
	glm::vec3 v23((p2.x - p3.x), (p2.y - p3.y), (p2.z - p3.z));

	glm::vec3 r(p.x, p.y, p.z);
	glm::vec3 r4(p3.x, p3.y, p3.z);

	glm::vec3 rr4 = r - r4;

	glm::mat3x3 T(v03, v13, v23);
	glm::mat3x3 T_inverse = glm::inverse(T);
	glm::vec3 lambda = T_inverse * (r - r4);

	if (lambda.x > -1e-3 && lambda.y > -1e-3 && lambda.z > -1e-3 && (lambda.x + lambda.y + lambda.z) < 1.001)
	{
		return true;
    }

    return false;
}

bool GeoUtil::IsPointInTetrahedron(const Vertex& p, const Vertex& p0, const Vertex& p1, const Vertex& p2, const Vertex& p3, glm::vec3& uvw)
{
    glm::vec3 v03((p0.x - p3.x), (p0.y - p3.y), (p0.z - p3.z));
    glm::vec3 v13((p1.x - p3.x), (p1.y - p3.y), (p1.z - p3.z));
    glm::vec3 v23((p2.x - p3.x), (p2.y - p3.y), (p2.z - p3.z));

    glm::vec3 r(p.x, p.y, p.z);
    glm::vec3 r4(p3.x, p3.y, p3.z);

    glm::vec3 rr4 = r - r4;

    glm::mat3x3 T(v03, v13, v23);
    glm::mat3x3 T_inverse = glm::inverse(T);
    glm::vec3 lambda = T_inverse * (r - r4);

    if (lambda.x > -1e-3 && lambda.y > -1e-3 && lambda.z > -1e-3 && (lambda.x + lambda.y + lambda.z) < 1.001)
    {
        uvw.x = lambda.x;
        uvw.y = lambda.y;
        uvw.z = lambda.z;
        return true;
    }

    return false;
}

bool GeoUtil::IsPointInTetrahedron_Robust(const Vertex& p, const Vertex& p0, const Vertex& p1, const Vertex& p2, const Vertex& p3, glm::vec4& uvwx)
{
    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(4, 4);
    A(0,0) = p0.x; A(0,1) = p1.x; A(0,2) = p2.x; A(0,3) = p3.x;
    A(1,0) = p0.y; A(1,1) = p1.y; A(1,2) = p2.y; A(1,3) = p3.y;
    A(2,0) = p0.z; A(2,1) = p1.z; A(2,2) = p2.z; A(2,3) = p3.z;
    A(3,0) = 1.0;  A(3,1) = 1.0;  A(3,2) = 1.0;  A(3,3) = 1.0;

    Eigen::VectorXd b = Eigen::VectorXd::Zero(4);
    b(0) = p.x; b(1) = p.y; b(2) = p.z; b(3) = 1.0;
    Eigen::VectorXd lambda = A.llt().solve(b);
    if (lambda(0) > 0 && lambda(1) > 0 && lambda(2) > 0 && lambda(3) > 0
            && lambda(0) + lambda(1) + lambda(2) + lambda(3) <= 1.000001
            && lambda(0) + lambda(1) + lambda(2) + lambda(3) >= 0.999999)
    {
        uvwx.x = lambda(0);
        uvwx.y = lambda(1);
        uvwx.z = lambda(2);
        uvwx.w = lambda(3);
        return true;
    }

    return false;
}

bool GeoUtil::IsPointInTetrahedron_Robust_test(const Vertex& p, const Vertex& p0, const Vertex& p1, const Vertex& p2, const Vertex& p3, glm::vec4& uvwx)
{
    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(3, 3);
    A(0,0) = p0.x - p3.x; A(0,1) = p1.x - p3.x; A(0,2) = p2.x - p3.x;
    A(1,0) = p0.y - p3.y; A(1,1) = p1.y - p3.y; A(1,2) = p2.y - p3.y;
    A(2,0) = p0.z - p3.z; A(2,1) = p1.z - p3.z; A(2,2) = p2.z - p3.z;

    Eigen::VectorXd b = Eigen::VectorXd::Zero(3);
    b(0) = p.x - p3.x;
    b(1) = p.y - p3.y;
    b(2) = p.z - p3.z;
    Eigen::VectorXd lambda = A.transpose()*b;
    if (lambda(0) > 0 && lambda(1) > 0 && lambda(2) > 0
            && lambda(0) + lambda(1) + lambda(2) <= 1.000001)
    {
        uvwx.x = lambda(0);
        uvwx.y = lambda(1);
        uvwx.z = lambda(2);
        uvwx.w = 1 - lambda(0) - lambda(1) - lambda(2);
        return true;
    }

    return false;
}

bool GeoUtil::IsVertexInsideTetrahedron(const Vertex& vertex, const Cell& tet, const std::vector<Vertex>& vecVertex, glm::vec3& lambda)
{
    bool bInside = false;
    glm::vec4 v1(vecVertex[tet.at(0)].x, vecVertex[tet.at(0)].y, vecVertex[tet.at(0)].z, 1.0);
    glm::vec4 v2(vecVertex[tet.at(1)].x, vecVertex[tet.at(1)].y, vecVertex[tet.at(1)].z, 1.0);
    glm::vec4 v3(vecVertex[tet.at(2)].x, vecVertex[tet.at(2)].y, vecVertex[tet.at(2)].z, 1.0);
    glm::vec4 v4(vecVertex[tet.at(3)].x, vecVertex[tet.at(3)].y, vecVertex[tet.at(3)].z, 1.0);

    glm::vec4 v(vertex.x, vertex.y, vertex.z, 1.0);

    glm::mat4x4 m0(v1, v2, v3, v4);
    glm::mat4x4 m1(v, v2, v3, v4);
    glm::mat4x4 m2(v1, v, v3, v4);
    glm::mat4x4 m3(v1, v2, v, v4);
    glm::mat4x4 m4(v1, v2, v3, v);

    double d0 = glm::determinant(m0);
    double d1 = glm::determinant(m1);
    double d2 = glm::determinant(m2);
    double d3 = glm::determinant(m3);
    double d4 = glm::determinant(m4);

//    if (d0 < 1e-7)
//        std::cout << "Tet degenerated!" << std::endl;
//
//    if (fabs(d1) < 1e-7 && fabs(d2) < 1e-7 && fabs(d3) < 1e-7 && fabs(d3) < 1e-7){
//        lambda = glm::vec3(0.25, 0.25, 0.25);
//        return true;
//    }

    if (((d0 > 0 && d1 > 0 && d2 > 0 && d3 > 0 && d4 > 0) || (d0 < 0 && d1 < 0 && d2 < 0 && d3 < 0 && d4 < 0))
        && fabs(d0 - d1 - d2 - d3 - d4) < 1e-7 && fabs(d0) > 1e-7)
    {
        lambda = glm::vec3(d1/d0, d2/d0, d3/d0);
        bInside = true;
    }

    return bInside;
}

bool GeoUtil::IsVertexInsideTetrahedron_Robust(const Vertex& p, const Vertex& p0, const Vertex& p1, const Vertex& p2, const Vertex& p3, glm::vec3& lambda)
{
    bool bInside = false;
    glm::vec4 v1(p0.x, p0.y, p0.z, 1.0);
    glm::vec4 v2(p1.x, p1.y, p1.z, 1.0);
    glm::vec4 v3(p2.x, p2.y, p2.z, 1.0);
    glm::vec4 v4(p3.x, p3.y, p3.z, 1.0);

    glm::vec4 v(p.x, p.y, p.z, 1.0);

    glm::mat4x4 m0(v1, v2, v3, v4);
    glm::mat4x4 m1(v, v2, v3, v4);
    glm::mat4x4 m2(v1, v, v3, v4);
    glm::mat4x4 m3(v1, v2, v, v4);
    glm::mat4x4 m4(v1, v2, v3, v);

    std::vector<double> d(5);
    d[0] = glm::determinant(m0);
    d[1] = glm::determinant(m1);
    d[2] = glm::determinant(m2);
    d[3] = glm::determinant(m3);
    d[4] = glm::determinant(m4);

    int hasPostive = 0;
    int hasNegtive = 0;
    for (size_t i = 0; i < 5; i++)
    {
        if (d[i] > 0)
            hasPostive = 1;
        else if (d[i] < 0)
            hasNegtive = -1;
    }

    bool isSameSign = (hasPostive == 1) ^ (hasNegtive == -1);
    if (isSameSign && fabs(d[0]) > 1e-8 && fabs(d[0] - d[1] - d[2] - d[3] - d[4]) < 1e-9)
    {
        lambda = glm::vec3(d[1]/d[0], d[2]/d[0], d[3]/d[0]);
        bInside = true;
    }

    return bInside;
}

bool GeoUtil::IsTetrahedronDegenerated(const Cell& tet, const std::vector<Vertex>& vecVertex, const double eps/* = 1e-9*/)
{
    glm::vec4 v1(vecVertex[tet.at(0)].x, vecVertex[tet.at(0)].y, vecVertex[tet.at(0)].z, 1.0);
    glm::vec4 v2(vecVertex[tet.at(1)].x, vecVertex[tet.at(1)].y, vecVertex[tet.at(1)].z, 1.0);
    glm::vec4 v3(vecVertex[tet.at(2)].x, vecVertex[tet.at(2)].y, vecVertex[tet.at(2)].z, 1.0);
    glm::vec4 v4(vecVertex[tet.at(3)].x, vecVertex[tet.at(3)].y, vecVertex[tet.at(3)].z, 1.0);
    glm::mat4x4 m0(v1, v2, v3, v4);
    double d0 = glm::determinant(m0);
    if (fabs(d0) < eps)
    {
//        std::cout << "volume: " << d0;
        return true;
    }
    else
        return false;
}

FACE_TYPE GeoUtil::GetFaceType(const glm::vec3& normal)
{
	const float l = glm::length(normal);
	if (fabs(normal.x/l) > 0.8)
	{
		return FACE_X;
	}
	else if (fabs(normal.y/l) > 0.8)
	{
		return FACE_Y;
	}
	else if (fabs(normal.z/l) > 0.8)
	{
		return FACE_Z;
	}

	FACE_TYPE faceType = FACE_X;
    double max = fabs(normal.x/l);
    if (fabs(normal.y/l) > max)
    {
    	max = fabs(normal.y/l);
    	faceType = FACE_Y;

    }
    if (fabs(normal.z/l) > max)
    {
    	max = fabs(normal.z/l);
    	faceType = FACE_Z;
    }
	return faceType;
}
