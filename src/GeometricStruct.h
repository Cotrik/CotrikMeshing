#ifndef __GEOMETRIC_STRUCT_H__
#define __GEOMETRIC_STRUCT_H__

#include "string.h"
#include <math.h>
#include <cmath>
#include <string>
#include <vector>
#include <set>
#include <algorithm>
#include "glm/glm.hpp"
#include "glm/gtc/matrix_access.hpp"

enum MESH_TYPE
{
	MESH_TYPE_HEXAHEDRON = 0,
	MESH_TYPE_TETRAHEDRON_VTK,
	MESH_TYPE_HEXAHEDRON_VTK,
	MESH_TYPE_TRIANGLE_OFF,
	MESH_TYPE_QUAD_OFF,
	MESH_TYPE_TETRAHEDRON_OFF,
	MESH_TYPE_TRIANGLE_OBJ,
	MESH_TYPE_TRIANGLE_VTK,
	MESH_TYPE_HEXAHEDRON_OFF,
	MESH_TYPE_TRIANGLE_MESH,
	MESH_TYPE_TETRAHEDRON_MESH
};
enum FileType
{
    FILE_VTK,
    FILE_OFF,
    FILE_MESH,
    FILE_OBJ,
    FILE_STL,
    FILE_PLY
};

enum ElementType
{
    POLYGON,
    TRIANGLE,
    QUAD,
    TETRAHEDRA,
    HEXAHEDRA
};
extern std::string strAjustedPolycubeTetMeshFilename;
enum DIRECTION
{
    X_AXIS = 0,    // face normal direction
    Y_AXIS,
    Z_AXIS,

    X_LINE,    // boundary line direction
    Y_LINE,
    Z_LINE,

    CORNER,   // corners

    X_AXIS_PLUS,    // face normal direction
    Y_AXIS_PLUS,
    Z_AXIS_PLUS,

    X_AXIS_MINUS,    // face normal direction
    Y_AXIS_MINUS,
    Z_AXIS_MINUS,

    UNKNOWN_AXIS,
};
struct VertexInfo
{
    VertexInfo()
    : bCorner(false)
    , bBoundaryLine(false)
    , bSurface(false)
    , bBoundaryCell(false)
    , bSingurality(false)
    , bValid(true)
    , bNeedSmoothing(true)
    , bAjusted(false)
    , bUsed(false)
    , curvature(0.0f)
    , dir(UNKNOWN_AXIS)
    , id(0)
    , duplicatedVertexIndex(0)
    {}
    VertexInfo(const VertexInfo& v)
    : bCorner(v.bCorner)
    , bBoundaryLine(v.bBoundaryLine)
    , bSurface(v.bSurface)
    , bBoundaryCell(v.bBoundaryCell)
    , bSingurality(v.bSingurality)
    , bValid(v.bValid)
    , bNeedSmoothing(v.bNeedSmoothing)
    , bAjusted(v.bAjusted)
    , bUsed(v.bUsed)
    , curvature(v.curvature)
    , dir(v.dir)
    , id(v.id)
    , duplicatedVertexIndex(v.duplicatedVertexIndex)
    {}

    bool bCorner;
    bool bBoundaryLine;
    bool bSurface;
    bool bBoundaryCell;
    bool bSingurality;
    bool bValid;            // if false, the vertex is not belong to the mesh
    bool bNeedSmoothing;
    bool bAjusted;
    bool bUsed;
    float curvature;
    DIRECTION dir;
    unsigned long id;
    unsigned long duplicatedVertexIndex;

    glm::vec3 compensation;
    std::vector<float> scalars;
    std::vector<glm::vec3> vectors;
    std::vector<glm::mat3x3> tensors;
};

struct CellInfo
{
    CellInfo()
    : bCorner(false)
    , bBoundary(false)
    , bSurface(false)
    , dir(UNKNOWN_AXIS)
    , id(0)
    {}
    CellInfo(const CellInfo& v)
    : bCorner(v.bCorner)
    , bBoundary(v.bBoundary)
    , bSurface(v.bSurface)
    , dir(v.dir)
    , id(v.id)
    {}

    bool bCorner;
    bool bBoundary;
    bool bSurface;
    DIRECTION dir;
    unsigned long id;

    std::vector<float> scalars;
    std::vector<glm::vec3> vectors;
    std::vector<glm::mat3x3> tensors;
};

struct VertexDouble
{
	VertexDouble()
		: x(0), y(0), z(0), a(0), bSurface(false), bBoundary(false), index(0)
	{

	}
	VertexDouble(const VertexDouble& vertex)
		: x(vertex.x), y(vertex.y), z(vertex.z), a(vertex.a),
		  bSurface(vertex.bSurface), bBoundary(vertex.bBoundary), bUsed(vertex.bUsed), index(vertex.index)
	{

	}
	VertexDouble& operator=(const VertexDouble& vertex)
	{
		x = vertex.x;
		y = vertex.y;
		z = vertex.z;
		a = vertex.a;
		bSurface = vertex.bSurface;
		bBoundary = vertex.bBoundary;
		index = vertex.index;
		bUsed = vertex.bUsed;
		// VertexDouble(vertex);
		return *this;
	}
	VertexDouble(double x, double y, double z)
		: x(x), y(y), z(z), a(0), bSurface(false), bBoundary(false), bUsed(false), index(0)
	{

	}
	// Subtract
	VertexDouble operator - (const VertexDouble& v) const
	{
		return VertexDouble(x - v.x, y - v.y, z - v.z);
	}

	const VertexDouble& operator + (const VertexDouble& v)
	{
		x += v.x;
		y += v.y;
		z += v.z;

		return *this;
	}

	// Dot product
	double Dot(const VertexDouble& v) const
	{
		return x*v.x +y*v.y + z*v.z;
	}
	// Cross product
	VertexDouble Cross(const VertexDouble& v) const
	{
		return VertexDouble(y*v.z - z*v.y, z*v.x - x*v.z, x*v.y - y*v.x);
	}
	double x;
	double y;
	double z;
	int a;
	bool bSurface;
	bool bBoundary;
	bool bUsed;
	int index;
};

struct VertexFloat
{
	VertexFloat()
		: x(0), y(0), z(0)
	{
	}
	VertexFloat(const VertexFloat& vertex)
		: x(vertex.x), y(vertex.y), z(vertex.z), vinfo(vertex.vinfo)
	{

	}
	const VertexFloat& operator=(const VertexFloat& vertex)
	{
		x = vertex.x;
		y = vertex.y;
		z = vertex.z;
		return *this;
	}
	VertexFloat(const float x, const float y, const float z)
		: x(x), y(y), z(z)
	{
	}
	const VertexFloat& operator=(const glm::vec3& v)
	{
		x = v.x;
		y = v.y;
		z = v.z;
		return *this;
	}
	const bool operator == (const VertexFloat& v) const
	{
		if (fabs(x - v.x) < 1e-6 && fabs(y - v.y) < 1e-6 && fabs(z - v.z) < 1e-6)
			return true;
		return false;
	}
	const bool operator < (const VertexFloat& v) const
	{
		if (x < v.x)
			return true;
		else if (fabs(x - v.x) < 1e-6)
		{
			if (y < v.y)
				return true;
			else if (fabs(y - v.y) < 1e-6)
			{
				if (z < v.z)
					return true;
			}
		}
		return false;
	}
	VertexFloat(const glm::vec3& v)
		: x(v.x), y(v.y), z(v.z)
	{
	}
	// Subtract
	VertexFloat operator - (const VertexFloat& v)
	{
		return VertexFloat(x - v.x, y - v.y, z - v.z);
	}

	const VertexFloat& operator + (const VertexFloat& v)
	{
		x += v.x;
		y += v.y;
		z += v.z;

		return *this;
	}
	// Dot product
	double Dot(const VertexFloat& v) const
	{
		return x*v.x +y*v.y + z*v.z;
	}
	// Cross product
	VertexFloat Cross(const VertexFloat& v) const
	{
		return VertexFloat(y*v.z - z*v.y, z*v.x - x*v.z, x*v.y - y*v.x);
	}

	double Length() const
	{
		return 0.5 * sqrt(x*x + y*y + z*z);
	}

	double Distance(const VertexFloat& v) const
	{
		return sqrt((x - v.x)*(x - v.x) + (y - v.y)*(y - v.y) + (z - v.z)*(z - v.z));
	}

	float x;
	float y;
	float z;

	VertexInfo vinfo;
};

typedef VertexFloat Vertex;

class _Cell : public std::vector<unsigned long>
{
public:
    CellInfo cinfo;
};
//typedef Cell Face;
typedef std::vector<unsigned long> Face;
typedef std::vector<unsigned long> Cell;

class HexahedronVertices
{
public:
	Vertex p0;
	Vertex p1;
	Vertex p2;
	Vertex p3;
	Vertex p4;
	Vertex p5;
	Vertex p6;
	Vertex p7;
};

class TetrahedronVertices
{
public:
	Vertex p0;
	Vertex p1;
	Vertex p2;
	Vertex p3;
};

class TriangleVertices
{
public:
	Vertex p0;
	Vertex p1;
	Vertex p2;
};

class LineVertices
{
public:
	LineVertices(const Vertex& p0, const Vertex& p1)
		: p0(p0), p1(p1)
	{

	}
	Vertex p0;
	Vertex p1;
};

class Line
{
public:
	Line(const unsigned long p0, const unsigned long p1)
		: p0(p0), p1(p1)
	{

	}
	Line(const Line& l)
		: p0(l.p0), p1(l.p1)
	{

	}
	Line()
		: p0(0), p1(0)
	{

	}
	virtual ~Line()
	{

	}
	unsigned long p0;
	unsigned long p1;

	virtual const bool operator == (const Line& right) const
	{
		return (((p0 == right.p0) && (p1 == right.p1)) || ((p0 == right.p1) && (p1 == right.p0)));
	}
	virtual const bool operator < (const Line& right) const
	{
		//return ((p0 < right.p0)/* || (p1 < right.p1) || (p0 < right.p1) || (p1 < right.p0)*/);
		if (p0 < right.p0)
		{
			return true;
		}
		else if (p0 == right.p0)
		{
			if (p1 < right.p1)
			{
				return true;
			}
		}
		return false;
	}
};
typedef Line Edge;

class DirectedEdge : public Line
{
public:
	DirectedEdge(const unsigned long p0, const unsigned long p1, const unsigned long cellIndex = 0)
		: Line(p0, p1), cellIndex(cellIndex)
	{

	}
	DirectedEdge(const DirectedEdge& l)
		: Line(l.p0, l.p1), cellIndex(l.cellIndex)
	{

	}
	DirectedEdge()
		: Line(0, 0), cellIndex(0)
	{

	}
	virtual ~DirectedEdge()
	{

	}
	virtual const bool operator == (const DirectedEdge& right) const
	{
		return ((p0 == right.p0) && (p1 == right.p1))/* || ((p0 == right.p1) && (p1 == right.p0))*/;
	}
	virtual const bool operator < (const DirectedEdge& right) const
	{
		bool bRet = false;
		if (p0 < right.p0)
			bRet = true;
		else if (p0 == right.p0)
			if (p1 < right.p1)
				bRet = true;
		return bRet;
	}
	unsigned long cellIndex;
};


class Quad
{
public:
	Quad(const unsigned long p0, const unsigned long p1, const unsigned long p2, const unsigned long p3)
		: p0(p0), p1(p1), p2(p2), p3(p3)
	{

	}
	Quad(const Quad& quad)
		: p0(quad.p0), p1(quad.p1), p2(quad.p2), p3(quad.p3)
	{

	}
	const bool operator == (const Quad& right) const
	{
		return ((p0 == right.p0 && p1 == right.p1 && p2 == right.p2 && p3 == right.p3));
	}
	unsigned long p0;
	unsigned long p1;
	unsigned long p2;
	unsigned long p3;
};

class QuadVertices
{
public:
	Vertex p0;
	Vertex p1;
	Vertex p2;
	Vertex p3;
};

/*
         3____________________2
         /|                 /|
        / |                / |
       /  |               /  |
   0  ____|_______________ 1 |
      |   |              |   |
      |   |              |   |
      |   |              |   |
      |   ____________________               |
      |   / 7            |  / 6
      |  /               | /
      | /                |/
      ___________________
    4                    5
*/
extern const unsigned int HexFaces[6][4];
extern const unsigned int TetFaces[4][3];

extern const unsigned int HexPoint_Points[8][3];
extern const unsigned int HexPoint_Faces[8][3][4];

extern const unsigned int TetPoint_Points[4][3];
extern const unsigned int TetPoint_Faces[4][3][3];

extern const unsigned int QuadPoint_Points[4][2];
extern const unsigned int TriPoint_Points[3][2];

extern const unsigned int HexEdge[12][2];
extern const unsigned int TetEdge[6][2];
extern const unsigned int QuadEdge[4][2];
extern const unsigned int TriEdge[3][2];

extern const unsigned int HexTrippleEdge[12][6];
struct Hexahedron
{
	Hexahedron()
		: p0(0), p1(0), p2(0), p3(0), p4(0), p5(0), p6(0), p7(0), a(0), bUsed(false)
	{
//		std::vector<int> f0;
//		f0.push_back(0);face.
//		f0.push_back(1);
//		f0.push_back(2);
//		f0.push_back(3);
//		std::vector<int> f1;
//		f1.push_back(1);
//		f1.push_back(5);
//		f1.push_back(6);
//		f1.push_back(2);
//		std::vector<int> f2;
//		f2.push_back(2);
//		f2.push_back(6);
//		f2.push_back(7);
//		f2.push_back(3);
//		std::vector<int> f3;
//		f3.push_back(3);
//		f3.push_back(7);
//		f3.push_back(4);
//		f3.push_back(0);
//		std::vector<int> f4;
//		f4.push_back(0);
//		f4.push_back(4);
//		f4.push_back(5);
//		f4.push_back(1);
//		std::vector<int> f5;
//		f5.push_back(7);
//		f5.push_back(6);
//		f5.push_back(5);
//		f5.push_back(4);
//		face.push_back(f0);
//		face.push_back(f1);
//		face.push_back(f2);
//		face.push_back(f3);
//		face.push_back(f4);
//		face.push_back(f5);
	}

	unsigned long p0;
	unsigned long p1;
	unsigned long p2;
	unsigned long p3;
	unsigned long p4;
	unsigned long p5;
	unsigned long p6;
	unsigned long p7;
	unsigned int a;
	bool bUsed;
//	std::vector<std::vector<int> > face;
};

struct Tetrahedron
{
	Tetrahedron()
		: p0(0), p1(0), p2(0), p3(0), a(0), bBoundary(false)
	{

	}
	Tetrahedron(const unsigned long p0, const unsigned long p1, const unsigned long p2, const unsigned long p3)
		: p0(p0), p1(p1), p2(p2), p3(p3), a(0), bBoundary(false)
	{

	}

	Tetrahedron& operator=(const Tetrahedron& right)
	{
		p0 = right.p0;
		p1 = right.p1;
		p2 = right.p2;
		p3 = right.p3;
		a = right.a;
		bBoundary = right.bBoundary;

		return (*this);
	}
	unsigned long p0;
	unsigned long p1;
	unsigned long p2;
	unsigned long p3;
	unsigned int a;
	bool bBoundary;
};

class Triangle
{
public:
	Triangle(const unsigned long p0, const unsigned long p1, const unsigned long p2)
		: p0(p0), p1(p1), p2(p2), a(0)
	{

	}
	Triangle(const Triangle& triangle)
		: p0(triangle.p0), p1(triangle.p1), p2(triangle.p2), a(0)
	{

	}

	Triangle()
		: p0(0), p1(0), p2(0), a(0)
	{

	}

	const bool operator == (const Triangle& right) const
	{
		return (
			(p0 == right.p0 && p1 == right.p1 && p2 == right.p2)
// 			|| (p0 == right.p1 && p1 == right.p2 && p2 == right.p0)
// 			|| (p0 == right.p2 && p1 == right.p0 && p2 == right.p1)
             );
	}
	const bool operator < (const Triangle& right) const
	{
		bool bRet = false;
		if (!bRet && p0 < right.p0)
		{
			bRet = true;
			return bRet;
		}
		if (!bRet && p1 < right.p1)
		{
			bRet = true;
			return bRet;
		}
		if (!bRet && p2 < right.p2)
		{
			bRet = true;
			return bRet;
		}
		return bRet;
	}

	const bool IsInTetrahedron(const Tetrahedron& tetrahedron)
	{
		std::vector<unsigned long> setTriangle;
		setTriangle.push_back(p0);
		setTriangle.push_back(p1);
		setTriangle.push_back(p2);

		std::vector<unsigned long> setTetrahedron;
		setTetrahedron.push_back(tetrahedron.p0);
		setTetrahedron.push_back(tetrahedron.p1);
		setTetrahedron.push_back(tetrahedron.p2);
		setTetrahedron.push_back(tetrahedron.p3);

		std::sort (setTriangle.begin(), setTriangle.end());
		std::sort (setTetrahedron.begin(), setTetrahedron.end());

        std::vector<unsigned long> v(7);
		std::vector<unsigned long>::iterator it = std::set_intersection (setTriangle.begin(), setTriangle.end(), setTetrahedron.begin(), setTetrahedron.end(), v.begin());
		// 10 20 0  0  0  0  0  0  0  0
		v.resize(it - v.begin());                      // 10 20

		if (v.size() == 3)
		{
			return true;
		}

		return false;
	}

	unsigned long p0;
	unsigned long p1;
	unsigned long p2;
	unsigned int a;
};

struct Keywords
{
	const char* pointKeyword;
	const char* cellKeyword;
};

struct MeshTypeKeywords
{
	MESH_TYPE meshType;
	Keywords keywords;
};


struct Point
{
	Point()
		: x(0), y(0), z(0)
	{

	}
	Point(double x, double y, double z)
		: x(x), y(y), z(z)
	{

	}
	double x;
	double y;
	double z;
};

typedef Point RotationAngle;

class Vector
{
public:
	Vector()
		: i(0), j(0), k(0)
	{

	}
	Vector(double i, double j, double k)
		: i(i), j(j), k(k)
	{

	}

	Vector(const Vector& v)
		: i(v.i), j(v.j), k(v.k)
	{

	}
	const bool operator < (const Vector& right) const
	{
		return i < right.i;
	}
	double i;
	double j;
	double k;
};

// ax + by + cz + d = 0;
class Plane
{
public:
	Plane(const double a, const double b, const double c, const double d)
		: a(a), b(b), c(c), d(d)
	{

	}

	Plane(const Vertex& p1, const Vertex& p2, const Vertex& p3)
	{
		a = ((p2.y - p1.y) * (p3.z - p1.z) - (p2.z - p1.z) * (p3.y - p1.y));
		b = ((p2.z - p1.z) * (p3.x - p1.x) - (p2.x - p1.x) * (p3.z - p1.z));
		c = ((p2.x - p1.x) * (p3.y - p1.y) - (p2.y - p1.y) * (p3.x - p1.x));
		d = (0 - (a * p1.x + b * p1.y + c * p1.z));
	}

	Plane(const Point& point, const Vector& vec)
		: a(vec.i), b(vec.j), c(vec.k), d(-vec.i*point.x - vec.j*point.y - vec.k*point.y)
	{
		double regularScale = sqrt(a*a + b*b + c*c);
		a /= regularScale;
		b /= regularScale;
		c /= regularScale;
		d /= regularScale;
	}

	Plane(const Plane& plane)
	{
		memcpy(this, &plane, sizeof(plane));
	}

	double IntersectionAngle(const Plane& p)
	{
		double cosangle = fabs(a*p.a + b*p.b + c*p.c) / (sqrt(a*a + b*b + c*c)* sqrt(p.a*p.a + p.b*p.b + p.c*p.c));
		return cosangle;
	}

	~Plane(){}
	virtual bool operator < (const Plane& right) const
	{
		return d < right.d;
	}

public:
	double a;
	double b;
	double c;
	double d;
};


class DPlane : public Plane
{
public:
	virtual ~DPlane();
	virtual bool operator < (const Plane& right) const
	{
		return d < right.d;
	}
};


struct Matrix4X4
{
	double q[4][4];
};


class Quaternions
{
public:
	Quaternions(Matrix4X4& matrix)
	{
		GenerateQuaternions();
	}
private:
	// Generate Coherently-orienting Quaternions
	// When constructing a quaternion from a rotation matrix R, we first reorder and reflect the columns of R to maximize its trace. This process guarantees that the scalar
	// components qw of the quaternions are strictly positive in all cases, which ensures that ||q|| > 0 and thus the quaternion is non-degenerate.
	void GenerateQuaternions(void)
	{
		// Get the trace of Matrix R
		double trace = 0;
	}
	double q[4][4];
};



struct Position_Distance
{
	Position_Distance()
	: distance(0.0f)
	{}

	bool operator < (const Position_Distance& right) const
	{
		return distance < right.distance;
	}

	glm::vec3 position;
	float distance;

};

struct TetHexParas
{
	std::vector<glm::vec3> paras;
	std::vector<unsigned long> tetIndex;
	std::vector<bool> flag;
	std::vector<unsigned long> count;
	std::vector<unsigned long> closestTetVertexIndex;
	std::vector<bool> visited;
};

struct TriQuadParas
{
	std::vector<glm::vec2> paras;
	std::vector<unsigned long> triIndex;
	std::vector<bool> flag;
	std::vector<unsigned long> count;
	std::vector<unsigned long> closestTriVertexIndex;
	std::vector<bool> visited;
};

struct OverlappingParas
{
	unsigned long hexIndex;
	std::vector<unsigned long> tetIndex;
	std::vector<glm::vec3> paras;
	std::vector<unsigned long> closestTetVertexIndex;
};

extern const unsigned long INVALID_NUM;

extern const double X_RM[35][3][3];
extern const double Y_RM[35][3][3];
extern const double Z_RM[35][3][3];


enum FACE_TYPE
{
	FACE_UNKNOWN = 0,
	FACE_X = 1,
	FACE_Y = 2,
	FACE_Z = 3,
	FACE_SHARED_2,
	FACE_SHARED_3,
	FACE_OTHER
};


struct Tripe
{
	unsigned long a;
	unsigned long b;
	unsigned long c;
};

extern const Tripe hexTripe[8][3];
extern const unsigned long hexTet[5][4];
#endif // __GEOMETRIC_STRUCT_H__
