#ifndef __Mesh_File_Writer_H_
#define __Mesh_File_Writer_H_

#include <vector>
#include <string>
#include <fstream>

#include "GeometricStruct.h"
#include "Iterator.h"
#include "Mesh.h"

typedef ConcreteAggregate<Vertex> Vertices;
typedef ConcreteAggregate<Vertex> Points;
typedef ConcreteAggregate<Hexahedron> Hexahedrons;
typedef ConcreteAggregate<Tetrahedron> Tetrahedrons;
typedef ConcreteIterator<Vertex> VertexIterator;
typedef ConcreteIterator<Hexahedron> HexahedronIterator;
typedef ConcreteIterator<Tetrahedron> TetrahedronIterator;

class MeshFileWriter
{
public:
	MeshFileWriter(const Mesh& mesh, const char* pFileName = "tempfile");
	MeshFileWriter(const std::vector<Vertex>& v, const std::vector<Cell>& c, const char* pFileName, const ElementType cellType = HEXAHEDRA);
	MeshFileWriter();
	~MeshFileWriter();

	void WriteTriangleSurfaceFile(const char* pSurfaceFilename = NULL);
	void WriteBoundaryTetrahedronFile(std::vector<Vertex*>& vecVertex, std::vector<Tetrahedron>& vecBoundaryTetrahedron,
			const char* pBoundaryFilename = NULL);
	void WriteQuadSurfaceFile(std::vector<Vertex*>& vecVertex, std::vector<Quad>& vecQuadSurfaceIndex,
			const char* pSurfaceFilename = NULL);
	void WriteTetrahedralSurfaceCurvature(const char* curvatureFilename, const char* tetFilename);
	void WriteMeshDesityFieldFile(const char* pFileName = NULL, const MESH_TYPE meshType = MESH_TYPE_TETRAHEDRON_VTK);
	void WriteMeshBoundaryInfo(const char* pFileName = NULL);
	void WriteCellData(const std::vector<int>& cellData);
	void WritePointData(const std::vector<int>& pointData);
public:
	void WriteFile();
	void WriteMeshFile();
	void WriteVtkFile();
	void WriteVtkPolyDataFile();
	void WriteOffFile();
	void WriteObjFile();
	void WriteStlFile();
	void WritePlyFile();

	void WriteSingularitiesVtkFile();

    void SetVertexInfo(const std::vector<VertexInfo>& vi);
    void SetFixFlag(bool bFixed = true);
	void FixMesh();
private:
	std::string m_strFileName;
	Mesh m_mesh;
	bool m_bFixed;
};

#endif // __Mesh_File_Writer_H_
