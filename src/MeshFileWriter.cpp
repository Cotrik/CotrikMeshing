#include "MeshFileWriter.h"
#include <iostream>
#include <iomanip>
#include <memory>
#include "Log.h"
#include "stdio.h"
using namespace std;
#include <boost/smart_ptr.hpp>

#define EPSINON 0.00001

MeshFileWriter::MeshFileWriter()
{

}
MeshFileWriter::~MeshFileWriter()
{

}

const MeshTypeKeywords meshTypeKeywords[] =
{
	{MESH_TYPE_HEXAHEDRON, "Vertices", "Hexahedra"},
	{MESH_TYPE_TETRAHEDRON_VTK, "POINTS", "CELLS"},
	{MESH_TYPE_HEXAHEDRON_VTK, "POINTS", "CELLS"},
	{MESH_TYPE_TRIANGLE_OFF, "", ""},
	{MESH_TYPE_QUAD_OFF, "", ""},
	{MESH_TYPE_TETRAHEDRON_OFF, "", ""},
	{MESH_TYPE_TRIANGLE_OBJ, "", ""},
	{MESH_TYPE_TRIANGLE_VTK, "POINTS", "CELLS"},
	{MESH_TYPE_HEXAHEDRON_OFF, "", ""},
	{MESH_TYPE_TRIANGLE_MESH, "Vertices", "Triangles"},
	{MESH_TYPE_TETRAHEDRON_MESH, "Vertices", "Tetrahedra"}
};

MeshFileWriter::MeshFileWriter(const Mesh& mesh, const char* pFileName)
: m_strFileName(pFileName)
, m_mesh(mesh)
, m_bFixed(false)
{

}

MeshFileWriter::MeshFileWriter(const std::vector<Vertex>& v, const std::vector<Cell>& c,
    const char* pFileName, const ElementType cellType/* = HEXAHEDRA*/)
: m_strFileName(pFileName)
, m_mesh(v, c, cellType)
, m_bFixed(false)
{

}

void MeshFileWriter::SetVertexInfo(const std::vector<VertexInfo>& vi)
{
    for (size_t i = 0; i < vi.size(); i++)
	m_mesh.V.at(i).vinfo = vi.at(i);
}
void MeshFileWriter::WriteFile()
{
    if (m_bFixed) FixMesh();
    if (m_strFileName.find(".vtk") != m_strFileName.npos)       WriteVtkFile();
    else if (m_strFileName.find(".off") != m_strFileName.npos)  WriteOffFile();
    else if (m_strFileName.find(".mesh") != m_strFileName.npos) WriteMeshFile();
    else if (m_strFileName.find(".obj") != m_strFileName.npos)  WriteObjFile();
    else if (m_strFileName.find(".stl") != m_strFileName.npos)  WriteStlFile();
}
void MeshFileWriter::WriteMeshFile()
{
    const std::vector<Vertex>& V = m_mesh.V;
    const std::vector<Cell>& C = m_mesh.C;
    const size_t vnum = V.size();
    const size_t cnum = C.size();

    std::ofstream ofs(m_strFileName.c_str());
    ofs << "MeshVersionFormatted 2" << endl;
    ofs << "Dimension 3" << endl;
    ofs << "Vertices " << vnum << endl;

    for (size_t i = 0; i < vnum; i++)
        ofs << std::fixed << setprecision(7) << V.at(i).x << " " << V.at(i).y << " " << V.at(i).z << " 0" << endl;

    if (m_mesh.m_cellType == TRIANGLE) ofs << "Triangles ";
    else if (m_mesh.m_cellType == QUAD) ofs << "Quadrilaterals ";
    else if (m_mesh.m_cellType == TETRAHEDRA) ofs << "Tetrahedra ";
    else if (m_mesh.m_cellType == HEXAHEDRA) ofs << "Hexahedra ";
    ofs << cnum << std::endl;

    for (size_t i = 0; i < cnum; i++){
        for (size_t j = 0; j < C.at(i).size(); j++)
            ofs << C.at(i).at(j) + 1 << " ";
        ofs << "0" << std::endl;
    }
    ofs << "End" << std::endl;
}
void MeshFileWriter::WriteSingularitiesVtkFile()
{
    const std::vector<Vertex>& V = m_mesh.V;
    const std::vector<Cell>& C = m_mesh.C;
    const size_t vnum = V.size();
    const size_t cnum = C.size();

    std::ofstream ofs(m_strFileName.c_str());
    ofs << "# vtk DataFile Version 2.0" << endl
        << m_strFileName.c_str() << endl
        << "ASCII" << endl << endl
        << "DATASET POLYDATA" << endl;
    ofs << "POINTS " << vnum << " float" << endl;
    for (size_t i = 0; i < vnum; i++)
        ofs << std::fixed << setprecision(7) << V.at(i).x << " " << V.at(i).y << " " << V.at(i).z << endl;
    size_t vnum_singularities = 0;

    for (size_t i = 0; i < vnum; i++)
        if (V.at(i).vinfo.bSingurality)
            vnum_singularities++;
    ofs << "VERTICES " << vnum_singularities << " " << 2 * vnum_singularities<< endl;
    for (size_t i = 0; i < vnum; i++)
        if (V.at(i).vinfo.bSingurality)
            ofs << "1 " << i << std::endl;

    const size_t snum = m_mesh.singularities.size();
    const std::vector<Edge>& E = m_mesh.E;
    const std::vector<unsigned long>& S = m_mesh.singularities;
    ofs << "LINES " << snum << " " << 3 * snum << endl;
    for (size_t i = 0; i < snum; i++)
        ofs << "2 " << E.at(S.at(i)).p0 << " " << E.at(S.at(i)).p1 << std::endl;
}
void MeshFileWriter::WriteVtkFile()
{
    const std::vector<Vertex>& V = m_mesh.V;
    const std::vector<Cell>& C = m_mesh.C;
    const size_t vnum = V.size();
    const size_t cnum = C.size();

    std::ofstream ofs(m_strFileName.c_str());
    ofs << "# vtk DataFile Version 2.0" << endl
        << m_strFileName.c_str() << endl
        << "ASCII" << endl << endl
        << "DATASET UNSTRUCTURED_GRID" << endl;
    ofs << "POINTS " << vnum << " float" << endl;
    for (size_t i = 0; i < vnum; i++)
        ofs << std::fixed << setprecision(7) << V.at(i).x << " " << V.at(i).y << " " << V.at(i).z << endl;
    ofs << "CELLS " << cnum << " ";

    vtkIdType idType = VTK_TRIANGLE;
    if (m_mesh.m_cellType == TRIANGLE) ofs << 4*cnum << std::endl;
    else if (m_mesh.m_cellType == QUAD) {idType = VTK_QUAD;  ofs << 5*cnum << std::endl;}
    else if (m_mesh.m_cellType == TETRAHEDRA) {idType = VTK_TETRA; ofs << 5*cnum << std::endl;}
    else if (m_mesh.m_cellType == HEXAHEDRA) {idType = VTK_HEXAHEDRON; ofs << 9*cnum << std::endl;}

    for (size_t i = 0; i < cnum; i++){
        ofs << C.at(i).size();
        for (size_t j = 0; j < C.at(i).size(); j++)
            ofs << " " << C.at(i).at(j);
        ofs << std::endl;
    }
    ofs << "CELL_TYPES " << cnum << endl;
    for (size_t i = 0; i < cnum; i++)
        ofs << idType << std::endl;
}

void MeshFileWriter::WriteVtkPolyDataFile()
{
    const std::vector<Vertex>& V = m_mesh.V;
    const std::vector<Cell>& C = m_mesh.C;
    const size_t vnum = V.size();
    const size_t cnum = C.size();

    std::ofstream ofs(m_strFileName.c_str());
    ofs << "# vtk DataFile Version 2.0" << endl
        << m_strFileName.c_str() << endl
        << "ASCII" << endl << endl
        << "DATASET POLYDATA" << endl;
    ofs << "POINTS " << vnum << " float" << endl;
    for (size_t i = 0; i < vnum; i++)
        ofs << std::fixed << setprecision(7) << V.at(i).x << " " << V.at(i).y << " " << V.at(i).z << endl;
    ofs << "POLYGONS " << cnum << " ";

    vtkIdType idType = VTK_TRIANGLE;
    if (m_mesh.m_cellType == TRIANGLE) ofs << 4*cnum << std::endl;
    else if (m_mesh.m_cellType == QUAD) {idType = VTK_QUAD;  ofs << 5*cnum << std::endl;}
    else if (m_mesh.m_cellType == TETRAHEDRA) {idType = VTK_TETRA; ofs << 5*cnum << std::endl;}
    else if (m_mesh.m_cellType == HEXAHEDRA) {idType = VTK_HEXAHEDRON; ofs << 9*cnum << std::endl;}

    for (size_t i = 0; i < cnum; i++){
        ofs << C.at(i).size();
        for (size_t j = 0; j < C.at(i).size(); j++)
            ofs << " " << C.at(i).at(j);
        ofs << std::endl;
    }
//    ofs << "CELL_TYPES " << cnum << endl;
//    for (size_t i = 0; i < cnum; i++)
//        ofs << idType << std::endl;
}
void MeshFileWriter::WriteOffFile()
{
    const std::vector<Vertex>& V = m_mesh.V;
    const std::vector<Cell>& C = m_mesh.C;
    const size_t vnum = V.size();
    const size_t cnum = C.size();

    std::ofstream ofs(m_strFileName.c_str());
    ofs << "OFF\n";
    ofs << vnum << " " << cnum << " 0\n";
    for (size_t i = 0; i < vnum; i++)
        ofs << std::fixed << setprecision(7) << V.at(i).x << " " << V.at(i).y << " " << V.at(i).z << endl;

    for (size_t i = 0; i < cnum; i++){
        if (m_mesh.m_cellType == TRIANGLE) ofs << 3;
        else if (m_mesh.m_cellType == QUAD) ofs << 4;
        else if (m_mesh.m_cellType == TETRAHEDRA) ofs << 5;
        else if (m_mesh.m_cellType == HEXAHEDRA) ofs << 10;

        for (size_t j = 0; j < C.at(i).size(); j++)
            ofs << " " << C.at(i).at(j);
        if (m_mesh.m_cellType == TETRAHEDRA) ofs << " 0\n";
        else if (m_mesh.m_cellType == HEXAHEDRA) ofs << " 0 0\n";
        else ofs << "\n";
    }
}

void MeshFileWriter::WriteObjFile()
{

}

void MeshFileWriter::WriteStlFile()
{

}

void MeshFileWriter::WriteTriangleSurfaceFile(const char* pSurfaceFilename)
{
	const std::vector<Cell>::iterator iterSurfaceBegin = m_mesh.surface.begin();
	const std::vector<Cell>::iterator iterSurfaceEnd = m_mesh.surface.end();

	const std::vector<Vertex>::iterator iterVertexBegin = m_mesh.V.begin();
	const std::vector<Vertex>::iterator iterVertexEnd = m_mesh.V.end();
	if (pSurfaceFilename == NULL)
	{
		cout << "surface filename is null" << endl;
		return;
	}
	std::ofstream surface(pSurfaceFilename);
	surface << "OFF" << endl;
	surface << m_mesh.V.size() << " " << m_mesh.surface.size() << " " << 0 << endl;

	// Write Vertices Coordinates
	for (std::vector<Vertex>::iterator iterVertex = iterVertexBegin; iterVertex != iterVertexEnd; ++iterVertex)
	{
		const Vertex& vertex = *iterVertex;
		surface << std::fixed << setprecision(5) << vertex.x << " " << vertex.y << " " << vertex.z << endl;
	}

	// Write Triangle Vertices Index
	for (std::vector<Cell>::iterator iterSurface = iterSurfaceBegin; iterSurface != iterSurfaceEnd; ++iterSurface)
	{
		surface << "3 " << iterSurface->at(0) << " " << iterSurface->at(1) << " " << iterSurface->at(2) << endl;
	}
	surface.close();
}

void MeshFileWriter::WriteTetrahedralSurfaceCurvature(const char* curvatureFilename, const char* tetFilename)
{
	std::ifstream curv_f(curvatureFilename);
	std::string strCurvatureFilename = std::string(tetFilename) + std::string(".curvature.vtk");
	std::ofstream curv_of(strCurvatureFilename.c_str(), std::ios_base::app);
	curv_of << "POINT_DATA " << m_mesh.V.size() << std::endl
			<< "SCALARS scalars float 1" << std::endl
            << "LOOKUP_TABLE default" << std::endl;
	float cur = 0;
	int i = 0;
	while (curv_f >> cur)
	{
		while(true)
		{
			if (m_mesh.V[i].vinfo.bSurface)
			{
				m_mesh.V[i].vinfo.curvature = cur;
				curv_of<< std::fixed << std::setprecision(5) << cur << std::endl;
				i++;
				break;
			}
			curv_of << "0.0" << std::endl;

			i++;
		}
	}
	curv_f.close();
	curv_of.close();
}
void MeshFileWriter::WriteMeshDesityFieldFile(const char* pFileName, const MESH_TYPE meshType)
{
	WriteFile();
	std::ofstream density(pFileName, std::ios_base::app);
	density << "POINT_DATA " << m_mesh.V.size() << std::endl
			<< "SCALARS scalars float 1" << std::endl
			<< "LOOKUP_TABLE default" << std::endl;
	for (int i = 0; i < m_mesh.vec_densityFiled.size(); i++)
	{
		density << m_mesh.vec_densityFiled[i] << std::endl;
	}
	density.close();
}

void MeshFileWriter::WriteCellData(const std::vector<int>& cellData)
{
    std::ofstream ofs(m_strFileName.c_str(), std::ios_base::app);
    ofs << "CELL_DATA " << cellData.size() << std::endl
            << "SCALARS scalars int 1" << std::endl
            << "LOOKUP_TABLE default" << std::endl;
    for (int i = 0; i < cellData.size(); i++)
    {
        ofs << cellData[i] << std::endl;
    }
    ofs.close();
}

void MeshFileWriter::WritePointData(const std::vector<int>& pointData)
{
    std::ofstream ofs(m_strFileName.c_str(), std::ios_base::app);
    ofs << "POINT_DATA " << pointData.size() << std::endl
            << "SCALARS point_scalars int 1" << std::endl
            << "LOOKUP_TABLE default" << std::endl;
    for (int i = 0; i < pointData.size(); i++)
    {
        ofs << pointData[i] << std::endl;
    }
    ofs.close();
}

void MeshFileWriter::WriteMeshBoundaryInfo(const char* pFileName)
{
	std::ofstream boundery(pFileName, std::ios_base::app);
	boundery << "POINT_DATA " << m_mesh.V.size() << std::endl
			<< "SCALARS scalars float 1" << std::endl
			<< "LOOKUP_TABLE default" << std::endl;
	for (int i = 0; i < m_mesh.V.size(); i++)
	{
		boundery << (m_mesh.V[i].vinfo.bSurface ? "0" : "1") << std::endl;
	}
	boundery.close();
}
void MeshFileWriter::WriteBoundaryTetrahedronFile(std::vector<Vertex*>& vecVertex, std::vector<Tetrahedron>& vecBoundaryTetrahedron, const char* pBoundaryFilename)
{
	std::vector<unsigned long> vecVertexSurfaceTriangle;
	const std::vector<Tetrahedron>::iterator iterTetrahedronBegin = vecBoundaryTetrahedron.begin();
	const std::vector<Tetrahedron>::iterator iterTetrahedronEnd = vecBoundaryTetrahedron.end();

	const std::vector<Vertex*>::iterator iterVertexBegin = vecVertex.begin();
	const std::vector<Vertex*>::iterator iterVertexEnd = vecVertex.end();
	if (pBoundaryFilename == NULL)
	{
		cout << "Boundary filename is null" << endl;
		return;
	}

	std::ofstream surface(pBoundaryFilename);
	surface << "OFF" << endl;
	surface << vecVertex.size() << " " << vecBoundaryTetrahedron.size() << " " << 0 << endl;
	for (std::vector<Vertex*>::iterator iterVertex = iterVertexBegin; iterVertex != iterVertexEnd; ++iterVertex)
	{
		Vertex* pVertex = *iterVertex;
		surface << std::fixed << setprecision(5) << pVertex->x << " " << pVertex->y << " " << pVertex->z << endl;
	}

	for (std::vector<Tetrahedron>::iterator iterTetrahedron = iterTetrahedronBegin; iterTetrahedron != iterTetrahedronEnd; ++iterTetrahedron)
	{
		surface << "4 " << iterTetrahedron->p0 << " " << iterTetrahedron->p1 << " " << iterTetrahedron->p2 << " " << iterTetrahedron->p3 << endl;
	}

	surface.close();
}

void MeshFileWriter::WriteQuadSurfaceFile(std::vector<Vertex*>& vecVertex, std::vector<Quad>& vecQuadSurfaceIndex, const char* pSurfaceFilename)
{
	std::vector<unsigned long> vecVertexSurfaceQuad;
	const std::vector<Quad>::iterator iterQuadBegin = vecQuadSurfaceIndex.begin();
	const std::vector<Quad>::iterator iterQuadEnd = vecQuadSurfaceIndex.end();

	const std::vector<Vertex*>::iterator iterVertexBegin = vecVertex.begin();
	const std::vector<Vertex*>::iterator iterVertexEnd = vecVertex.end();
	if (pSurfaceFilename == NULL)
	{
		cout << "Quad surface filename is null" << endl;
		return;
	}
	std::ofstream surface(pSurfaceFilename);
	surface << "OFF" << endl;
	surface << vecVertex.size() << " " << vecQuadSurfaceIndex.size() << " " << 0 << endl;
	for (std::vector<Vertex*>::iterator iterVertex = iterVertexBegin; iterVertex != iterVertexEnd; ++iterVertex)
	{
		Vertex* pVertex = *iterVertex;
		surface << std::fixed << setprecision(5) << pVertex->x << " " << pVertex->y << " " << pVertex->z << endl;
	}

	for (std::vector<Quad>::iterator iterQuad = iterQuadBegin; iterQuad != iterQuadEnd; ++iterQuad)
	{
		surface << "4 " << iterQuad->p0 << " " << iterQuad->p3 << " " << iterQuad->p2 << " " << iterQuad->p1 << endl;
	}
	surface.close();
}
void MeshFileWriter::SetFixFlag(bool bFixed)
{
	m_bFixed = bFixed;
}
void MeshFileWriter::FixMesh()
{
	std::vector<unsigned long> v_real_index;
	size_t c_size = 0;
	for (unsigned long i = 0; i < m_mesh.C.size(); i++)
	{
		if (m_mesh.C.at(i).size() > 8)
			continue;
		else
		{
			for (int j = 0; j < m_mesh.C.at(i).size(); j++)
			{
				v_real_index.push_back(m_mesh.C.at(i).at(j));
			}
			c_size++;
		}
	}

	std::sort(v_real_index.begin(), v_real_index.end());
	std::vector<unsigned long>::iterator iter = std::unique(v_real_index.begin(), v_real_index.end());
	v_real_index.resize(std::distance(v_real_index.begin(), iter));

	std::vector<Vertex> V(v_real_index.size());
	for (unsigned long i = 0; i < v_real_index.size(); i++)
	{
		const Vertex& v = m_mesh.V.at(v_real_index.at(i));
		V.at(i) = v;
	}
	m_mesh.V = V;
	//////////////////////////////////////////////////////
	std::map<unsigned long, unsigned long> v_v;
	unsigned long index = 0;
	for (unsigned long i = 0; i < v_real_index.size(); i++)
	{
		v_v[v_real_index.at(i)] = i;
	}

	std::vector<Cell> C;
	for (unsigned long i = 0; i < m_mesh.C.size(); i++)
	{
		Cell c;
		const size_t cellSize = m_mesh.C.at(i).size();
		if (cellSize == 8 || cellSize == 3 || cellSize == 4)
		{
			for (int j = 0; j < m_mesh.C.at(i).size(); j++)
			{
				const Cell& cell = m_mesh.C.at(i);
				c.push_back(v_v[cell.at(j)]);
			}
			C.push_back(c);
		}
	}
	C.resize(C.size());
	m_mesh.C = C;
}
