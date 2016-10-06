/*
 * Mesh.cpp
 *
 *  Created on: Dec 26, 2014
 *      Author: cotrik
 */
#include <iostream>
#include <fstream>
#include <algorithm>
#include "Mesh.h"
#include "GeoUtil.h"
#include "MeshFileWriter.h"
#include "include.h"
#include "glm/glm.hpp"
#include "glm/gtx/intersect.hpp"
double MAX_COORDINATE_VALUE = -1;

#ifdef _WIN32
HANDLE hMutex;
#else
pthread_mutex_t lock = PTHREAD_MUTEX_INITIALIZER;
#endif

int g_i;
std::vector<Cell>::iterator g_iterHexahedron;
std::vector<Cell>::iterator g_iterHexahedronBegin;
std::vector<Cell>::iterator g_iterHexahedronEnd;

int vecHexCellSize;
unsigned int THREAD_NUM = 4;

const glm::vec3 dirz(0.0, 0.0, 1.0);
const glm::vec3 diry(0.0, 1.0, 0.0);
const glm::vec3 dirx(1.0, 0.0, 0.0);

Mesh::Mesh()
: m_maxVertex(-1.0, -1.0, -1.0)
, m_minVertex(1.0, 1.0, 1.0)
, m_patchNumber(0)
, pUnstructuredGrid(NULL)
, m_genus(0)
{
	// TODO Auto-generated constructor stub

}

Mesh::~Mesh()
{
	// TODO Auto-generated destructor stub
}

Mesh::Mesh(const std::vector<Vertex>& v, const std::vector<Cell>& c, const ElementType cellType/* = HEXAHEDRA*/)
: V(v)
, C(c)
, m_cellType(cellType)
, m_genus(0)
{
	Init();
}

Mesh::Mesh(const Mesh& mesh)
: V(mesh.V)
, E(mesh.E)
, F(mesh.F)
, C(mesh.C)
, VI_VI(mesh.VI_VI)
, VI_EI(mesh.VI_EI)
, VI_FI(mesh.VI_FI)
, VI_CI(mesh.VI_CI)
, EI_FI(mesh.EI_FI)
, EI_CI(mesh.EI_CI)
, FI_CI(mesh.FI_CI)
, SF(mesh.SF)
, uniqueE(mesh.uniqueE)
, uniqueF(mesh.uniqueF)
, V_E(mesh.V_E)
, V_F(mesh.V_F)
//, V_C(mesh.V_C)
, surface(mesh.surface)
, vec_densityFiled(mesh.vec_densityFiled)
, m_meshType(mesh.m_meshType)
, m_cellType(mesh.m_cellType)
, pUnstructuredGrid(mesh.pUnstructuredGrid)
, m_genus(mesh.m_genus)
, scalarFields(mesh.scalarFields)
, vectorFields(mesh.vectorFields)
, singularities(mesh.singularities)
{
	//InitVertexInfo();
	//InitFace();
	Init();
}

bool IsPairExisted(const std::multimap<unsigned long, unsigned long>& V_V,
		const std::pair<unsigned long, unsigned long>& pair)
{
	const std::pair<std::multimap<unsigned long, unsigned long>::const_iterator,
			std::multimap<unsigned long, unsigned long>::const_iterator> ret =
			V_V.equal_range(pair.first);

    const std::multimap<unsigned long, unsigned long>::const_iterator iterLowerBound = ret.first;
    const std::multimap<unsigned long, unsigned long>::const_iterator iterUpperBound = ret.second;
    if (iterLowerBound == V_V.end() || iterUpperBound == V_V.end())
        return false;
    std::multimap<unsigned long, unsigned long>::const_iterator iterBound = iterLowerBound;

	// if the v_e is not existed, insert it to V_E
	for (; iterBound != iterUpperBound; ++iterBound)
	{
		const unsigned long& e = iterBound->second;
		if (e == pair.second)
		{
			return true;
			//break;
		}
	}

	return false;
}

bool IsPairExisted(const std::multimap<unsigned long, Edge>& V_E, const std::pair<unsigned long, Edge>& pair)
{
    const std::pair<std::multimap<unsigned long, Edge>::const_iterator,
    std::multimap<unsigned long, Edge>::const_iterator> ret = V_E.equal_range(pair.first);

    const std::multimap<unsigned long, Edge>::const_iterator iterLowerBound = ret.first;
    const std::multimap<unsigned long, Edge>::const_iterator iterUpperBound = ret.second;
    if (iterLowerBound == V_E.end() || iterUpperBound == V_E.end())
        return false;
    std::multimap<unsigned long, Edge>::const_iterator iterBound = iterLowerBound;

    // if the v_e is not existed, insert it to V_E
    for (; iterBound != iterUpperBound; ++iterBound)
    {
        const Edge& e = iterBound->second;
        if (e == pair.second)
        {
            return true;
            //break;
        }
    }

    return false;
}

bool IsPairExisted(const std::multimap<unsigned long, Cell>& V_C, const std::pair<unsigned long, Cell>& pair)
{
    const std::pair<std::multimap<unsigned long, Cell>::const_iterator,
    std::multimap<unsigned long, Cell>::const_iterator> ret = V_C.equal_range(pair.first);

    const std::multimap<unsigned long, Cell>::const_iterator iterLowerBound = ret.first;
    const std::multimap<unsigned long, Cell>::const_iterator iterUpperBound = ret.second;
    if (iterLowerBound == V_C.end() || iterUpperBound == V_C.end())
        return false;
    std::multimap<unsigned long, Cell>::const_iterator iterBound = iterLowerBound;

    // if the v_c is not existed, insert it to V_C
    for (; iterBound != iterUpperBound; ++iterBound)
    {
        const Cell& c = iterBound->second;
        Cell sortedc = c;
        Cell sortedp = pair.second;
        std::sort(sortedc.begin(), sortedc.end());
        std::sort(sortedp.begin(), sortedp.end());
        if (sortedc == sortedp)
        {
            return true;
            //break;
        }
    }

    return false;
}

bool IsPairExisted(const std::multimap<unsigned long, SFace>& V_F, const std::pair<unsigned long, SFace>& pair)
{
    const std::pair<std::multimap<unsigned long, SFace>::const_iterator,
    std::multimap<unsigned long, SFace>::const_iterator> ret = V_F.equal_range(pair.first);

    const std::multimap<unsigned long, SFace>::const_iterator iterLowerBound = ret.first;
    const std::multimap<unsigned long, SFace>::const_iterator iterUpperBound = ret.second;
    if (iterLowerBound == V_F.end() || iterUpperBound == V_F.end())
        return false;
    std::multimap<unsigned long, SFace>::const_iterator iterBound = iterLowerBound;

    // if the v_c is not existed, insert it to V_C
    for (; iterBound != iterUpperBound; ++iterBound)
    {
        const SFace& f = iterBound->second;
        if (f == pair.second)
        {
            return true;
        }
    }

    return false;
}

bool IsPairExisted(const std::multimap<Edge, unsigned long>& E_F, const std::pair<Edge, unsigned long>& pair)
{
    const std::pair<std::multimap<Edge, unsigned long>::const_iterator,
    std::multimap<Edge, unsigned long>::const_iterator> ret = E_F.equal_range(pair.first);

    const std::multimap<Edge, unsigned long>::const_iterator iterLowerBound = ret.first;
    const std::multimap<Edge, unsigned long>::const_iterator iterUpperBound = ret.second;
    if (iterLowerBound == E_F.end() || iterUpperBound == E_F.end())
        return false;
    std::multimap<Edge, unsigned long>::const_iterator iterBound = iterLowerBound;

    // if the e_f is not existed, insert it to E_F
    for (; iterBound != iterUpperBound; ++iterBound)
    {
        const unsigned long& f = iterBound->second;
        if (f == pair.second)
        {
            return true;
        }
    }

    return false;
}
void Mesh::Init()
{
    const std::vector<Vertex>::iterator iterVertexBegin = V.begin();
    const std::vector<Vertex>::iterator iterVertexEnd = V.end();
    std::vector<Vertex>::iterator iterV = iterVertexBegin;

    const std::vector<Cell>::iterator iterCellBegin = C.begin();
    const std::vector<Cell>::iterator iterCellEnd = C.end();
    std::vector<Cell>::iterator iterCell = iterCellBegin;
    GetF();
//    GetE();
//    GetVI_VI();
//    GetVI_EI();
//    GetVI_FI();
//    GetVI_CI();
//    GetEI_FI();
//    GetEI_CI();
//    GetFI_CI();
}
void Mesh::GetF()
{
    if (!F.empty())
        return;
    {
        if (m_cellType == TRIANGLE || m_cellType == QUAD)
        {
            F = C;
            surface = F;
        }
        else
        {
            const std::vector<Cell>::iterator iterCellBegin = C.begin();
            const std::vector<Cell>::iterator iterCellEnd = C.end();
            std::vector<Cell>::iterator iterCell = iterCellBegin;
            for (iterCell = iterCellBegin; iterCell != iterCellEnd; ++iterCell)
            {
                if (m_cellType == HEXAHEDRA)
                {
                    for (int i = 0; i < 6; i++)
                    {
                        Face face(4, 0);
                        for (int j = 0; j < 4; j++)
                            //face.push_back(iterCell->at(HexFaces[i][j]));
                            face.at(j) = iterCell->at(HexFaces[i][j]);
                        F.push_back(face);
                    }
                }
                else if (m_cellType == TETRAHEDRA)
                {
                    for (int i = 0; i < 4; i++)
                    {
                        Face face;
                        for (int j = 0; j < 3; j++)
                            face.push_back(iterCell->at(TetFaces[i][j]));
                        F.push_back(face);
                    }
                }
            }
            F.resize(F.size());
        }
    }
//    if (uniqueF.empty())
//    {
//        SF.resize(F.size());
//        std::copy(F.begin(), F.end(), SF.begin());
//        for (int i = 0; i < SF.size(); i++)
//        {
//            const SFace& f = SF.at(i);
//            std::vector<SFace>::iterator iterF = std::find(uniqueF.begin(), uniqueF.end(), f);
//            if (iterF == uniqueF.end())
//            {
//                uniqueF.push_back(f);
//            }
//        }
//    }

}
void Mesh::GetE()
{
    if (!E.empty())
        return;
    const std::vector<Cell>::iterator iterCellBegin = C.begin();
    const std::vector<Cell>::iterator iterCellEnd = C.end();
    std::vector<Cell>::iterator iterCell = iterCellBegin;

    {
        for (iterCell = iterCellBegin; iterCell != iterCellEnd; ++iterCell)
        {
            if (m_cellType == HEXAHEDRA)
            {
                Edge e;
                for (int i = 0; i < 12; i++)
                {
                    e.p0 = iterCell->at(HexEdge[i][0]);
                    e.p1 = iterCell->at(HexEdge[i][1]);
                    E.push_back(e);
                }
            }
            else if (m_cellType == TETRAHEDRA)
            {
                Edge e;
                for (int i = 0; i < 6; i++)
                {
                    e.p0 = iterCell->at(TetEdge[i][0]);
                    e.p1 = iterCell->at(TetEdge[i][1]);
                    E.push_back(e);
                }
            }
            else if (m_cellType == QUAD)
            {
                Edge e;
                for (int i = 0; i < 4; i++)
                {
                    e.p0 = iterCell->at(QuadEdge[i][0]);
                    e.p1 = iterCell->at(QuadEdge[i][1]);
                    E.push_back(e);
                }
            }
            else if (m_cellType == TRIANGLE)
            {
                Edge e;
                for (int i = 0; i < 3; i++)
                {
                    e.p0 = iterCell->at(TriEdge[i][0]);
                    e.p1 = iterCell->at(TriEdge[i][1]);
                    E.push_back(e);
                }
            }
        }
//        std::cout << "before ||E|| = " <<  E.size() << std::endl;
//        std::sort(E.begin(), E.end());
//        std::vector<Edge>::iterator iterE = std::unique(E.begin(), E.end());
//        E.resize(std::distance(E.begin(), iterE));
//        std::cout << "after ||E|| = " <<  E.size() << std::endl;
        if (uniqueE.empty())
            for (int i = 0; i < E.size(); i++)
            {
                const Edge& e = E.at(i);
                std::vector<Edge>::iterator iterE = std::find(uniqueE.begin(), uniqueE.end(), e);
                if (iterE == uniqueE.end())
                {
                    uniqueE.push_back(e);
                }
            }
        E = uniqueE;
    }
}

void Mesh::GetVI_VI()
{
    if (!VI_VI.empty())
        return;
    const std::vector<Cell>::iterator iterCellBegin = C.begin();
    const std::vector<Cell>::iterator iterCellEnd = C.end();
    std::vector<Cell>::iterator iterCell = iterCellBegin;
    for (iterCell = iterCellBegin; iterCell != iterCellEnd; ++iterCell)
    {
        const size_t cellSize = iterCell->size();
        for (unsigned long j = 0; j < cellSize; j++)
        {
            const unsigned long& v = iterCell->at(j);
            if (m_cellType == HEXAHEDRA) {
                for (int k = 0; k < 3; k++) {
                    const std::pair<unsigned long, unsigned long> v_v(v, iterCell->at(HexPoint_Points[j][k]));
                    if (!IsPairExisted(VI_VI, v_v))
                        VI_VI.insert(v_v);
                }
            }
            else if (m_cellType == TETRAHEDRA) {
                for (int k = 0; k < 3; k++) {
                    const std::pair<unsigned long, unsigned long> v_v(v, iterCell->at(TetPoint_Points[j][k]));
                    if (!IsPairExisted(VI_VI, v_v))
                        VI_VI.insert(v_v);
                }
            }
            else if (m_cellType == QUAD) {
                for (int k = 0; k < 2; k++) {
                    const std::pair<unsigned long, unsigned long> v_v(v, iterCell->at(QuadPoint_Points[j][k]));
                    if (!IsPairExisted(VI_VI, v_v))
                        VI_VI.insert(v_v);
                }
            }
            if (m_cellType == TRIANGLE)
            {
                for (int k = 0; k < 2; k++) {
                    const std::pair<unsigned long, unsigned long> v_v(v, iterCell->at(TriPoint_Points[j][k]));
                    if (!IsPairExisted(VI_VI, v_v))
                        VI_VI.insert(v_v);
                }
            }
        }
    }
    V_V = VI_VI;
}
void Mesh::GetVI_EI()
{
    if (!VI_EI.empty())
        return;
    const std::vector<Cell>::iterator iterCellBegin = C.begin();
    const std::vector<Cell>::iterator iterCellEnd = C.end();
    std::vector<Cell>::iterator iterCell = iterCellBegin;
    for (iterCell = iterCellBegin; iterCell != iterCellEnd; ++iterCell)
    {
        const size_t cellSize = iterCell->size();
        for (unsigned long j = 0; j < cellSize; j++)
        {
            const unsigned long& v = iterCell->at(j);
            if (m_cellType == HEXAHEDRA) {
                for (int k = 0; k < 3; k++) {
                    const Edge e(v, iterCell->at(HexPoint_Points[j][k]));
                    unsigned long ei = GetEI(e);
                    const std::pair<unsigned long, unsigned long> vi_ei(v, ei);
                    if (!IsPairExisted(VI_EI, vi_ei))
                        VI_EI.insert(vi_ei);
                }
            }
            else if (m_cellType == TETRAHEDRA) {
                for (int k = 0; k < 3; k++) {
                    const Edge e(v, iterCell->at(TetPoint_Points[j][k]));
                    unsigned long ei = GetEI(e);
                    const std::pair<unsigned long, unsigned long> vi_ei(v, ei);
                    if (!IsPairExisted(VI_EI, vi_ei))
                        VI_EI.insert(vi_ei);
                }
            }
            else if (m_cellType == QUAD) {
                for (int k = 0; k < 2; k++) {
                    const Edge e(v, iterCell->at(QuadPoint_Points[j][k]));
                    unsigned long ei = GetEI(e);
                    const std::pair<unsigned long, unsigned long> vi_ei(v, ei);
                    if (!IsPairExisted(VI_EI, vi_ei))
                        VI_EI.insert(vi_ei);
                }
            }
            else if (m_cellType == TRIANGLE) {
                for (int k = 0; k < 2; k++) {
                    const Edge e(v, iterCell->at(TriPoint_Points[j][k]));
                    unsigned long ei = GetEI(e);
                    const std::pair<unsigned long, unsigned long> vi_ei(v, ei);
                    if (!IsPairExisted(VI_EI, vi_ei))
                        VI_EI.insert(vi_ei);
                }
            }
        }
    }
}
void Mesh::GetVI_FI()
{
    if (!VI_FI.empty())
        return;
    const std::vector<Cell>::iterator iterCellBegin = C.begin();
    const std::vector<Cell>::iterator iterCellEnd = C.end();
    std::vector<Cell>::iterator iterCell = iterCellBegin;
    for (iterCell = iterCellBegin; iterCell != iterCellEnd; ++iterCell)
    {
        const size_t cellSize = iterCell->size();
        for (unsigned long j = 0; j < cellSize; j++)
        {
            const unsigned long& v = iterCell->at(j);
            //-------------------------------------------------
            // construct V_F
            if (m_cellType == HEXAHEDRA)
            {
                SFace f;
                f.resize(4);
                for (int k = 0; k < 3; k++) {
                    for (int m = 0; m < 4; m++)
                        f.at(m) = iterCell->at(HexPoint_Faces[j][k][m]);
                    const std::pair<unsigned long, SFace> v_f(v, f);
                    if (!IsPairExisted(V_F, v_f))
                        V_F.insert(v_f);
                }
            }
            else if (m_cellType == TETRAHEDRA) {
                for (int k = 0; k < 3; k++) {
                    SFace f;
                    for (int m = 0; m < 3; m++)
                        f.push_back(iterCell->at(TetPoint_Faces[j][k][m]));
                    const std::pair<unsigned long, SFace> v_f(v, f);
                    if (!IsPairExisted(V_F, v_f))
                        V_F.insert(v_f);
                }
            }
        }
    }

}
void Mesh::GetVI_CI()
{
    if (!VI_CI.empty())
        return;
    for (int i = 0; i < C.size(); i++) {
        for (unsigned long j = 0; j < C.at(i).size(); j++) {
            const std::pair<unsigned long, unsigned long> v_ci(C.at(i).at(j), i);
            if (!IsPairExisted(VI_CI, v_ci))
                VI_CI.insert(v_ci);
        }
    }
}
void Mesh::GetEI_FI()
{

}
void Mesh::GetEI_CI()
{
    if (!EI_CI.empty())
        return;
    const std::vector<Cell>::iterator iterCellBegin = C.begin();
    const std::vector<Cell>::iterator iterCellEnd = C.end();
    std::vector<Cell>::iterator iterCell = iterCellBegin;
    for (iterCell = iterCellBegin; iterCell != iterCellEnd; ++iterCell)
    {
        if (m_cellType == HEXAHEDRA)
        {
            Edge e;
            for (int i = 0; i < 12; i++)
            {
                e.p0 = iterCell->at(HexEdge[i][0]);
                e.p1 = iterCell->at(HexEdge[i][1]);
                unsigned long ei = GetEI(e);
                unsigned long ci = std::distance(iterCellBegin, iterCell);
                std::pair<unsigned long, unsigned long> p(ei, ci);
                EI_CI.insert(p);
            }
        }
        else if (m_cellType == TETRAHEDRA)
        {
            Edge e;
            for (int i = 0; i < 6; i++)
            {
                e.p0 = iterCell->at(TetEdge[i][0]);
                e.p1 = iterCell->at(TetEdge[i][1]);
                unsigned long ei = GetEI(e);
                unsigned long ci = std::distance(iterCellBegin, iterCell);
                std::pair<unsigned long, unsigned long> p(ei, ci);
                EI_CI.insert(p);
            }
        }
        else if (m_cellType == QUAD)
        {
            Edge e;
            for (int i = 0; i < 4; i++)
            {
                e.p0 = iterCell->at(QuadEdge[i][0]);
                e.p1 = iterCell->at(QuadEdge[i][1]);
                unsigned long ei = GetEI(e);
                unsigned long ci = std::distance(iterCellBegin, iterCell);
                std::pair<unsigned long, unsigned long> p(ei, ci);
                EI_CI.insert(p);
            }
        }
        else if (m_cellType == TRIANGLE)
        {
            Edge e;
            for (int i = 0; i < 3; i++)
            {
                e.p0 = iterCell->at(TriEdge[i][0]);
                e.p1 = iterCell->at(TriEdge[i][1]);
                unsigned long ei = GetEI(e);
                unsigned long ci = std::distance(iterCellBegin, iterCell);
                std::pair<unsigned long, unsigned long> p(ei, ci);
                EI_CI.insert(p);
            }
        }
    }
}
void Mesh::GetFI_CI()
{

}

unsigned long Mesh::GetEI(const Edge& e)
{
    std::vector<Edge>::iterator iter = std::find(E.begin(), E.end(), e);
    return std::distance(E.begin(), iter);
}
void Mesh::GetNeighboringInfo()
{
	const std::vector<Vertex>::iterator iterVertexBegin = V.begin();
	const std::vector<Vertex>::iterator iterVertexEnd = V.end();
	std::vector<Vertex>::iterator iterV = iterVertexBegin;

	const std::vector<Cell>::iterator iterCellBegin = C.begin();
	const std::vector<Cell>::iterator iterCellEnd = C.end();
	std::vector<Cell>::iterator iterCell = iterCellBegin;

	// Get All Faces
	if (F.empty()) {
		if (m_cellType == TRIANGLE || m_cellType == QUAD) {
			F = C;
			surface = F;
		}
		else {
			for (iterCell = iterCellBegin; iterCell != iterCellEnd; ++iterCell)	{
				if (m_cellType == HEXAHEDRA) {
					for (int i = 0; i < 6; i++) {
						Face face(4, 0);
						for (int j = 0; j < 4; j++)
							//face.push_back(iterCell->at(HexFaces[i][j]));
							face.at(j) = iterCell->at(HexFaces[i][j]);
						F.push_back(face);
					}
				}
				else if (m_cellType == TETRAHEDRA)
				{
					for (int i = 0; i < 4; i++)
					{
						Face face;
						for (int j = 0; j < 3; j++)
							face.push_back(iterCell->at(TetFaces[i][j]));
						F.push_back(face);
					}
				}
			}
			F.resize(F.size());
		}
	}
	const size_t VertexNum = Size();

	if (VI_CI.empty())
	{
		for (int i = 0; i < C.size(); i++)
		{
			for (unsigned long j = 0; j < C.at(i).size(); j++)
			{
				const std::pair<unsigned long, unsigned long> v_ci(C.at(i).at(j), i);
				if (!IsPairExisted(VI_CI, v_ci))
					VI_CI.insert(v_ci);
			}
		}
	}
	// construct V_E
	if (V_E.empty())
		for (iterCell = iterCellBegin; iterCell != iterCellEnd; ++iterCell)
		{
			const size_t cellSize = iterCell->size();
			for (unsigned long j = 0; j < cellSize; j++)
			{
				const unsigned long& v = iterCell->at(j);
				//-------------------------------------------------
				// construct V_F
				if (m_cellType == HEXAHEDRA)
				{
					SFace f;
					f.resize(4);
					for (int k = 0; k < 3; k++)
					{
						for (int m = 0; m < 4; m++)
							//f.push_back(iterCell->at(HexPoint_Faces[j][k][m]));
							f.at(m) = iterCell->at(HexPoint_Faces[j][k][m]);
						const std::pair<unsigned long, SFace> v_f(v, f);
						if (!IsPairExisted(V_F, v_f))
						{
							V_F.insert(v_f);
						}
					}
				}
				if (m_cellType == TETRAHEDRA)
				{
					for (int k = 0; k < 3; k++)
					{
						SFace f;
						for (int m = 0; m < 3; m++)
							f.push_back(iterCell->at(TetPoint_Faces[j][k][m]));
						const std::pair<unsigned long, SFace> v_f(v, f);
						if (!IsPairExisted(V_F, v_f))
						{
							V_F.insert(v_f);
						}
					}
				}
				//-------------------------------------------------
				// construct V_E
				if (m_cellType == HEXAHEDRA)
				{
					for (int k = 0; k < 3; k++)
					{
						const Edge e(v, iterCell->at(HexPoint_Points[j][k]));
						const std::pair<unsigned long, Edge> v_e(v, e);
						if (!IsPairExisted(V_E, v_e))
						{
							V_E.insert(v_e);
						}
						if (e.p0 == v)
						{
							const std::pair<unsigned long, unsigned long> v_v(v, e.p1);
							if (!IsPairExisted(V_V, v_v))
								V_V.insert(v_v);
						}
						else if (e.p1 == v)
						{
							const std::pair<unsigned long, unsigned long> v_v(v, e.p0);
							if (!IsPairExisted(V_V, v_v))
								V_V.insert(v_v);
						}
					}
				}

				if (m_cellType == TETRAHEDRA)
				{
					for (int k = 0; k < 3; k++)
					{
						const Edge e(v, iterCell->at(TetPoint_Points[j][k]));
						const std::pair<unsigned long, Edge> v_e(v, e);
						if (!IsPairExisted(V_E, v_e))
						{
							V_E.insert(v_e);
						}
						if (e.p0 == v)
						{
							const std::pair<unsigned long, unsigned long> v_v(v, e.p1);
							if (!IsPairExisted(V_V, v_v))
								V_V.insert(v_v);
						}
						else if (e.p1 == v)
						{
							const std::pair<unsigned long, unsigned long> v_v(v, e.p0);
							if (!IsPairExisted(V_V, v_v))
								V_V.insert(v_v);
						}
					}
				}
				if (m_cellType == TRIANGLE)
				{
					for (int k = 0; k < 2; k++)
					{
						const Edge e(v, iterCell->at(TriPoint_Points[j][k]));
						const std::pair<unsigned long, Edge> v_e(v, e);
						if (!IsPairExisted(V_E, v_e))
						{
							V_E.insert(v_e);
						}
						if (e.p0 == v)
						{
							const std::pair<unsigned long, unsigned long> v_v(v, e.p1);
							if (!IsPairExisted(V_V, v_v))
								V_V.insert(v_v);
						}
						else if (e.p1 == v)
						{
							const std::pair<unsigned long, unsigned long> v_v(v, e.p0);
							if (!IsPairExisted(V_V, v_v))
								V_V.insert(v_v);
						}
					}
				}
			}
		}
}


void Mesh::ExtractSurfaceStep(const Face& face)
{
	Face sortedFace(face.begin(), face.end());
	std::sort(sortedFace.begin(), sortedFace.end());
	// Look up the triangle in the set
	const std::vector<Face>::iterator iterBegin = sortedSurface.begin();
	const std::vector<Face>::iterator iterEnd = sortedSurface.end();
	std::vector<Face>::iterator iter = iterBegin;//std::find(iterBegin, iterEnd, sortedFace);

	for (; iter != iterEnd; ++iter)
	{
		bool bFlag = true;
		for (unsigned int i = 0; i < iter->size(); i++)
			bFlag &= (iter->at(i) == sortedFace.at(i));
		if (bFlag)
			break;
	}
	if (iter != iterEnd) // If the triangle is in the set, delete it
	{
		int index = std::distance(iterBegin, iter);
		iter = sortedSurface.erase(iter);
		surface.erase(surface.begin() + index);
	}
	else // If the triangle is not in the set, insert it
	{
		sortedSurface.push_back(sortedFace);
		surface.push_back(face);
	}
}

void Mesh::ExtractSurface()
{
	if (!surface.empty() || F.empty())
	{
		return;
	}

	std::vector<FacePair> vecFacePair;
	const std::vector<Face>::iterator iterBegin = F.begin();
	const std::vector<Face>::iterator iterEnd = F.end();
	for (std::vector<Face>::iterator iter = iterBegin; iter != iterEnd; ++iter)
	{
		//ExtractSurfaceStep(*iter);
		Face sortedFace(iter->begin(), iter->end());
		std::sort(sortedFace.begin(), sortedFace.end());

		FacePair facePair(sortedFace, *iter);
		vecFacePair.push_back(facePair);
	}
	std::sort(vecFacePair.begin(), vecFacePair.end());

	for (int i = 0; i < vecFacePair.size() - 1; i++)
	{
		if (vecFacePair.at(i) == vecFacePair.at(i + 1))
		{
			i++;
		}
		else
		{
			surface.push_back(vecFacePair.at(i).originFace);
		}
	}
	if (vecFacePair.at(vecFacePair.size() - 2)
			== vecFacePair.at(vecFacePair.size() - 1))
	{
		;
	}
	else
	{
		surface.push_back(vecFacePair.at(vecFacePair.size() - 1).originFace);
	}

//	std::for_each(F.begin(), F.end(), ExtractSurfaceStep);
	for (std::vector<Face>::const_iterator iter = surface.begin();
			iter != surface.end(); ++iter)
	{
		for (int i = 0; i < iter->size(); i++)
		{
			V[iter->at(i)].vinfo.bSurface = true;
		}
	}
	for (std::multimap<unsigned long, Edge>::const_iterator iter = V_E.begin();
			iter != V_E.end(); ++iter)
	{
		if (V[iter->first].vinfo.bSurface)
		{
			const Edge &e = iter->second;
			V[e.p0].vinfo.bBoundaryCell = true;
			V[e.p1].vinfo.bBoundaryCell = true;
		}
	}
	std::cout << "Done extracting surface. surface faces size = " << surface.size() << std::endl;
}

void Mesh::GetSingularities()
{
    if (!singularities.empty())
        return;
    if (m_cellType == HEXAHEDRA || m_cellType == TETRAHEDRA)
    {
        Init();
        ExtractSurface();
        for (size_t i = 0; i < E.size(); i++)
        {
            std::pair<std::multimap<unsigned long, unsigned long>::iterator,
                      std::multimap<unsigned long, unsigned long>::iterator> ret = EI_CI.equal_range(i);
            const unsigned long valence = std::distance(ret.first, ret.second);
            const Edge& e = E.at(i);
            if (  (valence != 4 && (!V.at(e.p0).vinfo.bSurface || !V.at(e.p1).vinfo.bSurface)) // internal Edge
                ||(valence != 2 && V.at(e.p0).vinfo.bSurface && V.at(e.p1).vinfo.bSurface) ) // surface Edge
            {
                V.at(e.p0).vinfo.bSingurality = true;
                V.at(e.p1).vinfo.bSingurality = true;
                singularities.push_back(i);
                //std::cout << "valence = " << valence << std::endl;
            }
        }
    }
}
void Mesh::NormalizeCoordinateValue()
{
	m_maxVertex = Vertex(-1.0, -1.0, -1.0);
	m_minVertex = Vertex(1.0, 1.0, 1.0);
	GetMaxMinCoordinates();
	float boxMaxLength = m_maxVertex.x - m_minVertex.x;
	if (m_maxVertex.y - m_minVertex.y > boxMaxLength)
	{
		boxMaxLength = m_maxVertex.y - m_minVertex.y;
	}
	if (m_maxVertex.z - m_minVertex.z > boxMaxLength)
	{
		boxMaxLength = m_maxVertex.z - m_minVertex.z;
	}

	glm::vec3 vecDiag(m_maxVertex.x - m_minVertex.x,
			m_maxVertex.y - m_minVertex.y, m_maxVertex.z - m_minVertex.z);
	float radius = glm::length(vecDiag) / 2;
	float radius_r = 2.0 / boxMaxLength;
	glm::vec3 centerPoint(m_minVertex.x + vecDiag.x / 2,
			m_minVertex.y + vecDiag.y / 2, m_minVertex.z + vecDiag.z / 2);

	std::cout << "centerPoint is (" << centerPoint.x << ", " << centerPoint.y
			<< ", " << centerPoint.z << ")" << std::endl;
	const std::vector<Vertex>::iterator iterBegin = V.begin();
	const std::vector<Vertex>::iterator iterEnd = V.end();
	std::vector<Vertex>::iterator iter = iterBegin;
	while (iter != iterEnd)
	{
		iter->x = (iter->x - centerPoint.x);
		iter->y = (iter->y - centerPoint.y);
		iter->z = (iter->z - centerPoint.z);
		++iter;
	}

	float maxCoordinateValue = 0;
	GetMaxMinCoordinates();
	iter = iterBegin;
	while (iter != iterEnd)
	{
		if (fabs(iter->x) > maxCoordinateValue)
		{
			maxCoordinateValue = fabs(iter->x);
		}
		if (fabs(iter->y) > maxCoordinateValue)
		{
			maxCoordinateValue = fabs(iter->y);
		}
		if (fabs(iter->z) > maxCoordinateValue)
		{
			maxCoordinateValue = fabs(iter->z);
		}
		++iter;
	}
	MAX_COORDINATE_VALUE = 2 * maxCoordinateValue;

	iter = iterBegin;
	while (iter != iterEnd)
	{
		iter->x = iter->x / MAX_COORDINATE_VALUE;
		iter->y = iter->y / MAX_COORDINATE_VALUE;
		iter->z = iter->z / MAX_COORDINATE_VALUE;
		++iter;
	}
}
void Mesh::GetCenterPointOfEdge(const Edge& e, glm::vec3& centerPoint)
{
	centerPoint.x = (V[e.p0].x + V[e.p1].x) / 2;
	centerPoint.y = (V[e.p0].y + V[e.p1].y) / 2;
	centerPoint.z = (V[e.p0].z + V[e.p1].z) / 2;
}

void Mesh::GetCurvatureCenterPointOfEdge(const Edge& e, glm::vec3& centerPoint)
{
	double w = (1.0 / V[e.p0].vinfo.curvature + 1.0 / V[e.p1].vinfo.curvature);
	centerPoint.x = (V[e.p0].x / V[e.p0].vinfo.curvature + V[e.p1].x / V[e.p1].vinfo.curvature) / w;
	centerPoint.y = (V[e.p0].y / V[e.p0].vinfo.curvature + V[e.p1].y / V[e.p1].vinfo.curvature) / w;
	centerPoint.z = (V[e.p0].z / V[e.p0].vinfo.curvature + V[e.p1].z / V[e.p1].vinfo.curvature) / w;
}

void Mesh::GetCenterPointOfCell(const Cell& c, glm::vec3& centerPoint)
{
	size_t cellSize = c.size();
	for (unsigned int i = 0; i < cellSize; i++)
	{
		centerPoint.x += V[c.at(i)].x;
		centerPoint.y += V[c.at(i)].y;
		centerPoint.z += V[c.at(i)].z;
	}
	centerPoint.x = centerPoint.x / cellSize;
	centerPoint.y = centerPoint.y / cellSize;
	centerPoint.z = centerPoint.z / cellSize;
}

double Mesh::GetCenterCurvatureOfCell(const Cell& c)
{
	const size_t cellSize = c.size();
	double curvature = 0.0;
	for (unsigned int i = 0; i < cellSize; i++)
	{
		curvature += V[c.at(i)].vinfo.curvature;
	}

	return curvature / cellSize;
}

void Mesh::GetCurvatureCenterPointOfCell(const Cell& c, glm::vec3& centerPoint)
{
	const size_t cellSize = c.size();
	double w = 0;
	for (unsigned int i = 0; i < cellSize; i++)
	{
		centerPoint.x += V[c.at(i)].x / V[c.at(i)].vinfo.curvature;
		centerPoint.y += V[c.at(i)].y / V[c.at(i)].vinfo.curvature;
		centerPoint.z += V[c.at(i)].z / V[c.at(i)].vinfo.curvature;
		w += 1.0 / V[c.at(i)].vinfo.curvature;
	}
	centerPoint.x = centerPoint.x / w;
	centerPoint.y = centerPoint.y / w;
	centerPoint.z = centerPoint.z / w;
}

void Mesh::GetCurvatureCenterPointOfCell_P(const Cell& c,
		glm::vec3& centerPoint)
{
	const size_t cellSize = c.size();
	double w = 0;
	for (unsigned int i = 0; i < cellSize; i++)
	{
		centerPoint.x += V[c.at(i)].x * V[c.at(i)].vinfo.curvature;
		centerPoint.y += V[c.at(i)].y * V[c.at(i)].vinfo.curvature;
		centerPoint.z += V[c.at(i)].z * V[c.at(i)].vinfo.curvature;
		w += V[c.at(i)].vinfo.curvature;
	}
	centerPoint.x = centerPoint.x / w;
	centerPoint.y = centerPoint.y / w;
	centerPoint.z = centerPoint.z / w;
}

bool Mesh::IsEdgeOnBondaryLine(const Edge& e)
{
	if (!V[e.p0].vinfo.bBoundaryLine || !V[e.p1].vinfo.bBoundaryLine)
	{
		return false;
	}

	return true;
}

bool Mesh::IsEdgeOnPatch(const Edge& e, const int patchLabel)
{
	const std::vector<int>& labels0 = vertexLabel.at(e.p0);
	std::vector<int>::const_iterator iter0 = std::find(labels0.begin(),
			labels0.end(), patchLabel);
	if (iter0 == labels0.end())
	{
		return false;
	}

	const std::vector<int>& labels1 = vertexLabel.at(e.p1);
	std::vector<int>::const_iterator iter1 = std::find(labels1.begin(),
			labels1.end(), patchLabel);
	if (iter1 == labels1.end())
	{
		return false;
	}

	return true;
}

bool Mesh::IsFaceOnSurface(const Face& f)
{
	size_t cellSize = f.size();
	for (unsigned int i = 0; i < cellSize; i++)
	{
		if (!V[f.at(i)].vinfo.bSurface)
		{
			return false;
		}
	}
	return true;
}

bool Mesh::IsFaceOnPatch(const Face& f, const int patchLabel)
{
	size_t cellSize = f.size();
	for (unsigned int i = 0; i < cellSize; i++)
	{
		const std::vector<int>& labels = vertexLabel.at(f.at(i));
		std::vector<int>::const_iterator iter = std::find(labels.begin(),
				labels.end(), patchLabel);
		if (iter == labels.end())
		{
			return false;
		}
	}
	return true;
}
const DIRECTION Mesh::GetFaceAxis(const Face& f)
{
	const glm::vec3 v0(V[f.at(0)].x, V[f.at(0)].y, V[f.at(0)].z);
	const glm::vec3 v1(V[f.at(1)].x, V[f.at(1)].y, V[f.at(1)].z);
	const glm::vec3 v2(V[f.at(2)].x, V[f.at(2)].y, V[f.at(2)].z);

	const glm::vec3 v10 = v0 - v1;
	const glm::vec3 v12 = v2 - v1;
	const glm::vec3 normal = glm::cross(v12, v10);

	DIRECTION dir = X_AXIS;
	float dir_value = 0;
	if (fabs(normal.x) > fabs(dir_value))
	{
		dir_value = fabs(normal.x);
		dir = X_AXIS;
	}
	if (fabs(normal.y) > fabs(dir_value))
	{
		dir_value = fabs(normal.y);
		dir = Y_AXIS;
	}
	if (fabs(normal.z) > fabs(dir_value))
	{
		dir_value = fabs(normal.z);
		dir = Z_AXIS;
	}

	return dir;
}

void Mesh::GetCurvature(const char* surfaceCurvatureFilename/* = NULL*/)
{
	if (surfaceCurvatureFilename != NULL)
	{
		std::ifstream curv_f(surfaceCurvatureFilename);
		float cur = 0;
		int i = 0;
		while (curv_f >> cur)
		{
			while (true)
			{
				if (V[i].vinfo.bSurface)
				{
					V[i].vinfo.curvature = cur;
					i++;
					break;
				}
				i++;
			}
		}
		curv_f.close();
	}
}

void Mesh::CurvatureWeightingSmoothSurface(const unsigned int time)
{
//	GetVertexInfo();
	int count = time;
	while (count-- != 0)
	{
		std::map<unsigned long, glm::vec3> map_point_meanCoordinates;
		std::multimap<unsigned long, SFace>::iterator iterBegin = V_F.begin();
		std::multimap<unsigned long, SFace>::iterator iterEnd = V_F.end();
		for (std::multimap<unsigned long, SFace>::iterator iter = iterBegin; iter != iterEnd;/* ++iter*/)
		{
			if (V[iter->first].vinfo.bSurface
				&& !V[iter->first].vinfo.bBoundaryLine)
			{
				int pointCount = 0;
				const unsigned long vertexIndex = iter->first;
				glm::vec3 mean(0, 0, 0);
				std::pair<std::multimap<unsigned long, SFace>::iterator,
						std::multimap<unsigned long, SFace>::iterator> ret =
						V_F.equal_range(iter->first);
				for (iter = ret.first; iter != ret.second; ++iter)
				{
					const Face& f = iter->second;
					if (IsFaceOnSurface(f))
					{
						glm::vec3 m(0.0f, 0.0f, 0.0f);
						GetCurvatureCenterPointOfCell_P(f, m);
						mean += m;
						pointCount++;
					}
				}
				mean.x = mean.x / pointCount;
				mean.y = mean.y / pointCount;
				mean.z = mean.z / pointCount;
				map_point_meanCoordinates.insert(
						std::pair<unsigned long, glm::vec3>(vertexIndex, mean));
			}
			else
			{
				++iter;
			}
		}

		std::map<unsigned long, glm::vec3>::iterator iterB =
				map_point_meanCoordinates.begin();
		std::map<unsigned long, glm::vec3>::iterator iterE =
				map_point_meanCoordinates.end();
		for (std::map<unsigned long, glm::vec3>::iterator iter = iterB;
				iter != iterE; ++iter)
		{
			Vertex& v = V[iter->first];
			v.x = 1.05 * iter->second.x;
			v.y = 1.05 * iter->second.y;
			v.z = 1.05 * iter->second.z;
		}
	}
}

void Mesh::GetMaxMinCoordinates()
{
	m_maxVertex.x = -1;
	m_maxVertex.y = -1;
	m_maxVertex.z = -1;

	m_minVertex.x = 1;
	m_minVertex.y = 1;
	m_minVertex.z = 1;
	for (std::vector<Vertex>::iterator iter = V.begin(); iter != V.end(); ++iter)
	{
		if (V.at(std::distance(V.begin(), iter)).vinfo.bSurface)
		GeoUtil::GetMaxMinCoordinateValueOfVertex(*iter, m_maxVertex, m_minVertex);
	}
}

void Mesh::GetMaxMinCoordinates_Patch()
{
	m_maxVertex.x = -1;
	m_maxVertex.y = -1;
	m_maxVertex.z = -1;

	m_minVertex.x = 1;
	m_minVertex.y = 1;
	m_minVertex.z = 1;
	for (std::vector<Patch>::iterator iter = m_patches.begin(); iter != m_patches.end(); ++iter)
	{
		if (iter->dir == X_AXIS)
		{
			if (iter->center.x > m_maxVertex.x) m_maxVertex.x = iter->center.x;
			if (iter->center.x < m_minVertex.x) m_minVertex.x = iter->center.x;
		}
		else if (iter->dir == Y_AXIS)
		{
			if (iter->center.y > m_maxVertex.y) m_maxVertex.y = iter->center.y;
			if (iter->center.y < m_minVertex.y) m_minVertex.y = iter->center.y;
		}
		else if (iter->dir == Z_AXIS)
		{
			if (iter->center.z > m_maxVertex.z) m_maxVertex.z = iter->center.z;
			if (iter->center.z < m_minVertex.z) m_minVertex.z = iter->center.z;
		}
//		std::cout << "Dir " << iter->dir << "Pos(" << iter->center.x  << ", " << iter->center.y << ", " << iter->center.z
//				<< " ) max( " << m_maxVertex.x << ", " << m_maxVertex.y << ", " << m_maxVertex.z
//				<< " ) min(" << m_minVertex.x << ", " << m_minVertex.y << ", " << m_minVertex.z << " )" << std::endl;
	}
}

void Mesh::GenerateSmallCubes_Magifier(const float& delta)
{
	for (std::vector<Vertex>::iterator iterPoint = V.begin(); iterPoint != V.end(); ++iterPoint)
	{
		GeoUtil::GetMaxMinCoordinateValueOfVertex(*iterPoint, m_maxVertex,	m_minVertex);
	}

	Vector axis(1.0, 1.0, 1.0);
	Vertex orig(-0.5, -0.5, -0.5);
	DivideHexahedronIntoSmallerOnes(orig, axis, delta, delta, delta,
			m_vecPolycubeHexVertex, m_vecPolycubeHexCell);
}

void Mesh::GenerateSmallCubes(const float& delta, const char* outputFilename)
{
	for (std::vector<Vertex>::iterator iterPoint = V.begin();
			iterPoint != V.end(); ++iterPoint)
	{
		GeoUtil::GetMaxMinCoordinateValueOfVertex(*iterPoint, m_maxVertex,	m_minVertex);
	}
	Vector axis(m_maxVertex.x - m_minVertex.x, m_maxVertex.y - m_minVertex.y, m_maxVertex.z - m_minVertex.z);
	DivideHexahedronIntoSmallerOnes(m_minVertex, axis, delta, delta, delta,
			m_vecPolycubeHexVertex, m_vecPolycubeHexCell);
}

void Mesh::GenerateSmallCubes_plus(const float& delta,
		const char* outputFilename, const float plus)
{
	for (std::vector<Vertex>::iterator iterPoint = V.begin();
			iterPoint != V.end(); ++iterPoint)
	{
		GeoUtil::GetMaxMinCoordinateValueOfVertex(*iterPoint, m_maxVertex,	m_minVertex);
	}
	Vector axis(m_maxVertex.x - m_minVertex.x, m_maxVertex.y - m_minVertex.y,
			m_maxVertex.z - m_minVertex.z);
	m_minVertex.x += plus;
	m_minVertex.y += plus;
	m_minVertex.z += plus;
	DivideHexahedronIntoSmallerOnes(m_minVertex, axis, delta, delta, delta,
			m_vecPolycubeHexVertex, m_vecPolycubeHexCell);
}

//bool IsVertexInsideTetrahedron(const Vertex& vertex, const Cell& tet, const std::vector<Vertex>& vecVertex)
//{
//	bool bInside = false;
//	glm::vec4 v1(vecVertex[tet.at(0)].x, vecVertex[tet.at(0)].y, vecVertex[tet.at(0)].z, 1.0);
//	glm::vec4 v2(vecVertex[tet.at(1)].x, vecVertex[tet.at(1)].y, vecVertex[tet.at(1)].z, 1.0);
//	glm::vec4 v3(vecVertex[tet.at(2)].x, vecVertex[tet.at(2)].y, vecVertex[tet.at(2)].z, 1.0);
//	glm::vec4 v4(vecVertex[tet.at(3)].x, vecVertex[tet.at(3)].y, vecVertex[tet.at(3)].z, 1.0);
//
//	glm::vec4 v(vertex.x, vertex.y, vertex.z, 1.0);
//
//	glm::mat4x4 m0(v1, v2, v3, v4);
//	glm::mat4x4 m1(v, v2, v3, v4);
//	glm::mat4x4 m2(v1, v, v3, v4);
//	glm::mat4x4 m3(v1, v2, v, v4);
//	glm::mat4x4 m4(v1, v2, v3, v);
//
//	float d0 = glm::determinant(m0);
//	float d1 = glm::determinant(m1);
//	float d2 = glm::determinant(m2);
//	float d3 = glm::determinant(m3);
//	float d4 = glm::determinant(m4);
//
//	if ((d0 >= 0 && d1 >= 0 && d2 >= 0 && d3 >= 0 && d4 >= 0)
//		|| (d0 <= 0 && d1 <= -0 && d2 <= 0 && d3 <= 0 && d4 <= 0))
//	{
//		bInside = true;
//	}
//
//	return bInside;
//}

bool IsVertexInsideTetrahedron(const Vertex& vertex, const Cell& tet,
		const std::vector<Vertex>& vecVertex)
{
	const Vertex& p0 = vecVertex[tet.at(0)];
	const Vertex& p1 = vecVertex[tet.at(1)];
	const Vertex& p2 = vecVertex[tet.at(2)];
	const Vertex& p3 = vecVertex[tet.at(3)];

	glm::vec3 v03((p0.x - p3.x), (p0.y - p3.y), (p0.z - p3.z));
	glm::vec3 v13((p1.x - p3.x), (p1.y - p3.y), (p1.z - p3.z));
	glm::vec3 v23((p2.x - p3.x), (p2.y - p3.y), (p2.z - p3.z));

	glm::vec3 r(vertex.x, vertex.y, vertex.z);
	glm::vec3 r4(p3.x, p3.y, p3.z);

	glm::vec3 rr4 = r - r4;

	glm::mat3x3 T(v03, v13, v23);
	glm::mat3x3 T_inverse = glm::inverse(T);
	glm::vec3 rambda = T_inverse * (r - r4);

	if (rambda.x > -1e-3 && rambda.y > -1e-3 && rambda.z > -1e-3
			&& (rambda.x + rambda.y + rambda.z) < 1.01)
	{
		return true;
	}
	else
	{
		return false;
	}
}

void Mesh::AddVertex(const std::vector<Vertex>& vecOrgTetVertex,
		Vertex& hexVertex, const Cell& tet,
		const std::vector<Vertex>& vecPolycubeVertex,
		std::vector<Vertex>& vecHexVertex)
{
	if (!hexVertex.vinfo.bUsed)
	{
		if (IsVertexInsideTetrahedron(hexVertex, tet, vecPolycubeVertex))
		{
			const Vertex& p0 = vecPolycubeVertex[tet.at(0)];
			const Vertex& p1 = vecPolycubeVertex[tet.at(1)];
			const Vertex& p2 = vecPolycubeVertex[tet.at(2)];
			const Vertex& p3 = vecPolycubeVertex[tet.at(3)];

			glm::vec3 v03((p0.x - p3.x), (p0.y - p3.y), (p0.z - p3.z));
			glm::vec3 v13((p1.x - p3.x), (p1.y - p3.y), (p1.z - p3.z));
			glm::vec3 v23((p2.x - p3.x), (p2.y - p3.y), (p2.z - p3.z));

			glm::mat3x3 T(v03, v13, v23);
			glm::vec3 r(hexVertex.x, hexVertex.y, hexVertex.z);
			glm::vec3 r4(p3.x, p3.y, p3.z);
			glm::mat3x3 T_inverse = glm::inverse(T);
			glm::vec3 rambda = T_inverse * (r - r4);
			//////////////////////////////////////////
			const Vertex& origT0 = vecOrgTetVertex[tet.at(0)];
			const Vertex& origT1 = vecOrgTetVertex[tet.at(1)];
			const Vertex& origT2 = vecOrgTetVertex[tet.at(2)];
			const Vertex& origT3 = vecOrgTetVertex[tet.at(3)];

			glm::vec3 p03((origT0.x - origT3.x), (origT0.y - origT3.y),
					(origT0.z - origT3.z));
			glm::vec3 p13((origT1.x - origT3.x), (origT1.y - origT3.y),
					(origT1.z - origT3.z));
			glm::vec3 p23((origT2.x - origT3.x), (origT2.y - origT3.y),
					(origT2.z - origT3.z));
			glm::mat3x3 OrigT(p03, p13, p23);
			glm::vec3 new_r4(origT3.x, origT3.y, origT3.z);
			glm::vec3 new_r = OrigT * rambda;
			new_r += new_r4;
			Vertex baryCenter(new_r.x, new_r.y, new_r.z);
			if (baryCenter.x != baryCenter.x || baryCenter.y != baryCenter.y
					|| baryCenter.z != baryCenter.z)
				cout << " Fail in parameterization!" << std::endl << "lambda("
						<< rambda.x << ", " << rambda.y << ", " << rambda.z
						<< ")" << std::endl;
			vecHexVertex.push_back(baryCenter);

			hexVertex.vinfo.id = vecHexVertex.size() - 1;
			hexVertex.vinfo.bUsed = true;
		}
	}
}
static long zerohexnum = 2;
void Mesh::CreateHexMesh(const Mesh& origTetMesh,
		std::vector<Vertex>& vecHexVertex, std::vector<Cell>& vecHexCell)
{
	const std::vector<Cell>::iterator iterHexBegin =
			m_vecPolycubeHexCell.begin();
	const std::vector<Cell>::iterator iterHexEnd = m_vecPolycubeHexCell.end();
//#pragma omp parallel for
	for (std::vector<Cell>::iterator iterHex = iterHexBegin;
			iterHex != iterHexEnd; ++iterHex)
	{
		Cell& hex = *iterHex;
		if (hex.size() > 8)
			continue;
		Vertex* p[8];
		for (int i = 0; i < 8; i++)
			p[i] = &m_vecPolycubeHexVertex[hex.at(i)];
		const std::vector<Cell>::iterator iterTetBegin = C.begin();
		const std::vector<Cell>::iterator iterTetEnd = C.end();
		for (std::vector<Cell>::iterator iterTet = iterTetBegin;
				iterTet != iterTetEnd; ++iterTet)
		{
			for (int i = 0; i < 8; i++)
				AddVertex(origTetMesh.V, *p[i], *iterTet, V, vecHexVertex);
			bool bUsed = true;
			for (int i = 0; i < 8; i++)
				bUsed &= p[i]->vinfo.bUsed;
			if (bUsed)
				break;
		}
		Cell hexCell;
		int zeroNum = 0;
		for (int i = 0; i < 8; i++)
		{
			hexCell.push_back(p[i]->vinfo.id);
			if (p[i]->vinfo.id == 0)
				zeroNum++;
		}
		if (zeroNum > 1/* && zerohexnum-- <= 0*/)
			continue;
		vecHexCell.push_back(hexCell);
		int index = std::distance(iterHexBegin, iterHex);
		if (index % 10000 == 0)
		{
			std::cout << "index = " << index << ". Hex Cell size = "
					<< m_vecPolycubeHexCell.size() << std::endl;
			std::cout << ((float) index) * 100 / m_vecPolycubeHexCell.size()
					<< "% Generated" << std::endl;
		}
	}
}

struct ThreadStructTriangle
{
	ThreadStructTriangle(
	/*const*/std::vector<Vertex>& vecPolycubeHexVertex,
	/*const*/std::vector<Vertex>& vecPolycubeTetVertex,
	/*const*/std::vector<Face>& polycubeTetSurface, const float delta = 0.008f,
			const DIRECTION dir = Z_AXIS) :
			vecPolycubeHexVertex(vecPolycubeHexVertex), vecPolycubeTetVertex(
					vecPolycubeTetVertex), polycubeTetSurface(
					polycubeTetSurface), delta(delta), dir(dir)
	{

	}
	std::vector<Vertex>& vecPolycubeHexVertex;
	std::vector<Vertex>& vecPolycubeTetVertex;
	std::vector<Face>& polycubeTetSurface;
	float delta;
	DIRECTION dir;
};

bool IsPointInsidePolyhedron(const std::vector<Vertex>& vecVertex,
		const std::vector<Face>& vecSurfaceTriangle, const glm::vec3& orig,
		const glm::vec3& dir)
{
	bool bInside = false;
	const std::vector<Face>::const_iterator iterBegin =
			vecSurfaceTriangle.begin();
	const std::vector<Face>::const_iterator iterEnd = vecSurfaceTriangle.end();
	for (std::vector<Face>::const_iterator iter = iterBegin; iter != iterEnd;
			++iter)
	{
		const Face& tri = *iter;
		const Vertex& vertex0 = vecVertex[tri.at(0)];
		const Vertex& vertex1 = vecVertex[tri.at(1)];
		const Vertex& vertex2 = vecVertex[tri.at(2)];
		const glm::vec3 v0(vertex0.x, vertex0.y, vertex0.z);
		glm::vec3 v1(vertex1.x, vertex1.y, vertex1.z);
		glm::vec3 v2(vertex2.x, vertex2.y, vertex2.z);
		glm::vec3 position;
		if (glm::intersectRayTriangle(orig, dir, v0, v1, v2, position)
				|| glm::intersectRayTriangle(orig, dir, v0, v2, v1, position))
		{
			bInside = !bInside;
		}
	}
	return bInside;
}

bool Mesh::IsPointInsideMesh(const Vertex& p, const Mesh& mesh)
{
	const unsigned long cubeIndex = mesh.GetCubeIndex(p, CUBESIZE);
	std::vector<unsigned long> vecTetIndex;
	for (unsigned long i = 0; i < mesh.V.size(); i++)
	{
		const unsigned long cube_index = mesh.GetCubeIndex(mesh.V.at(i),
				CUBESIZE);
		if (cubeIndex == cubeIndex)
		{
			const unsigned long v_idex = i;
			std::pair<
					std::multimap<unsigned long, unsigned long>::const_iterator,
					std::multimap<unsigned long, unsigned long>::const_iterator> iterBound =
					mesh.VI_CI.equal_range(v_idex);
			for (std::multimap<const unsigned long, unsigned long>::const_iterator iter =
					iterBound.first; iter != iterBound.second; ++iter)
			{
				vecTetIndex.push_back(iter->second);
			}
		}
	}

	std::sort(vecTetIndex.begin(), vecTetIndex.end());
	std::vector<unsigned long>::iterator iter = std::unique(vecTetIndex.begin(),
			vecTetIndex.end());
	vecTetIndex.resize(std::distance(vecTetIndex.begin(), iter));

	for (int i = 0; i < vecTetIndex.size(); i++)
	{
		const unsigned long tetIndex = vecTetIndex.at(i);
		const Cell& tet = mesh.C.at(tetIndex);

		const Vertex& p0 = mesh.V.at(tet.at(0));
		const Vertex& p1 = mesh.V.at(tet.at(1));
		const Vertex& p2 = mesh.V.at(tet.at(2));
		const Vertex& p3 = mesh.V.at(tet.at(3));

		if (GeoUtil::IsPointInTetrahedron(p, p0, p1, p2, p3))
			return true;
	}

	return false;
}

#ifdef _WIN32
DWORD WINAPI ProcTriangle_Magnifier(LPVOID lpParameter)
#else
void* ProcTriangle_Magnifier(void* lpParameter)
#endif
{
	const ThreadStructTriangle& p = *(ThreadStructTriangle*) lpParameter;
	glm::vec3 ray_direction = dirz;
	if (p.dir == X_AXIS)
	{
		ray_direction = dirx;
	}
	else if (p.dir == Y_AXIS)
	{
		ray_direction = diry;
	}

	while (g_iterHexahedron < g_iterHexahedronEnd)
	{
#ifdef _WIN32
		WaitForSingleObject(hMutex, INFINITE);
#else
		pthread_mutex_lock(&lock);
#endif
		const std::vector<Cell>::iterator iterHexahedron = g_iterHexahedron;
		Cell& hex = *iterHexahedron;
		++g_iterHexahedron;
		g_i++;
#ifdef _WIN32
		ReleaseMutex(hMutex);
#else
		pthread_mutex_unlock(&lock);
#endif
		if (hex.empty())
			continue;
// 		const Vertex& v2 = (p.vecPolycubeHexVertex.at(hex.at(2)));
//		const glm::vec3 origin(v2.x + p.delta/2, v2.y + p.delta/2, v2.z + p.delta/2);
//
//		if (!IsPointInsidePolyhedron(p.vecPolycubeTetVertex, p.polycubeTetSurface, origin, ray_direction))
//		{
//			hex.push_back(INVALID_NUM);
//		}
		bool bOneCubeVertexInsideTetMesh = false;
		for (unsigned int i = 0; i < hex.size(); i++)
		{
			const Vertex& v = p.vecPolycubeHexVertex.at(hex.at(i));
			const glm::vec3 origin(v.x, v.y, v.z);

			if (IsPointInsidePolyhedron(p.vecPolycubeTetVertex,
					p.polycubeTetSurface, origin, ray_direction))
			{
				bOneCubeVertexInsideTetMesh = true;
				break;
			}
		}
		if (!bOneCubeVertexInsideTetMesh)
			hex.push_back(INVALID_NUM);
		if (g_i % 10000 == 0)
		{
			std::cout << (((float) g_i * 100) / vecHexCellSize) << "% is done"
					<< std::endl;
		}
	}
#ifdef _WIN32
	return 0;
#else
	return NULL;
#endif
}

#ifdef _WIN32
DWORD WINAPI ProcTriangle(LPVOID lpParameter)
#else
void* ProcTriangle(void* lpParameter)
#endif
{
	const ThreadStructTriangle& p = *(ThreadStructTriangle*) lpParameter;
	glm::vec3 ray_direction = dirz;
	if (p.dir == X_AXIS)
	{
		ray_direction = dirx;
	}
	else if (p.dir == Y_AXIS)
	{
		ray_direction = diry;
	}

	while (g_iterHexahedron < g_iterHexahedronEnd)
	{
#ifdef _WIN32
		WaitForSingleObject(hMutex, INFINITE);
#else
		pthread_mutex_lock(&lock);
#endif
		const std::vector<Cell>::iterator iterHexahedron = g_iterHexahedron;
		Cell& hex = *iterHexahedron;
		++g_iterHexahedron;
		g_i++;
#ifdef _WIN32
		ReleaseMutex(hMutex);
#else
		pthread_mutex_unlock(&lock);
#endif
		if (hex.empty())
			continue;
		const Vertex& v2 = (p.vecPolycubeHexVertex.at(hex.at(2)));
		const Vertex& v4 = (p.vecPolycubeHexVertex.at(hex.at(4)));
		//const glm::vec3 origin(v2.x + p.delta/2, v2.y + p.delta/2, v2.z + p.delta/2);
		const glm::vec3 origin(0.5 * (v2.x + v4.x), 0.5 * (v2.y + v4.y),
				0.5 * (v2.z + v4.z));

		if (!IsPointInsidePolyhedron(p.vecPolycubeTetVertex,
				p.polycubeTetSurface, origin, ray_direction))
		{
			hex.push_back(INVALID_NUM);
		}
		if (g_i % 100000 == 0)
		{
			std::cout << (((float) g_i * 100) / vecHexCellSize) << "% is done"
					<< std::endl;
		}
	}
#ifdef _WIN32
	return 0;
#else
	return NULL;
#endif
}

struct ThreadRemoveOutsidePointsParameters
{
	ThreadRemoveOutsidePointsParameters(Mesh& mesh, const float delta = 0.008f) :
			mesh(mesh), delta(delta)
	{

	}

	Mesh& mesh;
	float delta;
};

#ifdef _WIN32
DWORD WINAPI ThreadRemoveOutsidePoints(LPVOID lpParameter)
#else
void* ThreadRemoveOutsidePoints(void* lpParameter)
#endif
{
	const ThreadRemoveOutsidePointsParameters& p =
			*(ThreadRemoveOutsidePointsParameters*) lpParameter;

	while (g_iterHexahedron < g_iterHexahedronEnd)
	{
#ifdef _WIN32
		WaitForSingleObject(hMutex, INFINITE);
#else
		pthread_mutex_lock(&lock);
#endif
		const std::vector<Cell>::iterator iterHexahedron = g_iterHexahedron;
		Cell& hex = *iterHexahedron;
		++g_iterHexahedron;
		g_i++;
#ifdef _WIN32
		ReleaseMutex(hMutex);
#else
		pthread_mutex_unlock(&lock);
#endif
		if (hex.empty())
			continue;
		const Vertex& v2 = (p.mesh.m_vecPolycubeHexVertex.at(hex.at(2)));
		//const glm::vec3 origin(v2.x + p.delta/2, v2.y + p.delta/2, v2.z + p.delta/2);
		const Vertex point(v2.x + p.delta / 2, v2.y + p.delta / 2,
				v2.z + p.delta / 2);

		if (!Mesh::IsPointInsideMesh(point, p.mesh))
		{
			hex.push_back(INVALID_NUM);
		}
		if (g_i % 10000 == 0)
		{
			std::cout << (((float) g_i * 100) / vecHexCellSize) << "% is done"
					<< std::endl;
		}
	}
#ifdef _WIN32
	return 0;
#else
	return NULL;
#endif
}

#ifdef _WIN32
DWORD WINAPI ProcTriangle_new(LPVOID lpParameter)
#else
void* ProcTriangle_new(void* lpParameter)
#endif
{
	const ThreadStructTriangle& p = *(ThreadStructTriangle*) lpParameter;
	glm::vec3 ray_direction = dirz;
	if (p.dir == X_AXIS)
	{
		ray_direction = dirx;
	}
	else if (p.dir == Y_AXIS)
	{
		ray_direction = diry;
	}

	while (g_iterHexahedron < g_iterHexahedronEnd)
	{
#ifdef _WIN32
		WaitForSingleObject(hMutex, INFINITE);
#else
		pthread_mutex_lock(&lock);
#endif
		const std::vector<Cell>::iterator iterHexahedron = g_iterHexahedron;
		Cell& hex = *iterHexahedron;
		++g_iterHexahedron;
		g_i++;
#ifdef _WIN32
		ReleaseMutex(hMutex);
#else
		pthread_mutex_unlock(&lock);
#endif
		if (hex.empty())
			continue;
		for (int i = 0; i < 8; i++)
		{
			const Vertex& v = (p.vecPolycubeHexVertex.at(hex.at(i)));
			const glm::vec3 origin(v.x, v.y, v.z);

			if (!IsPointInsidePolyhedron(p.vecPolycubeTetVertex,
					p.polycubeTetSurface, origin, ray_direction))
			{
				hex.push_back(INVALID_NUM);
				break;
			}
		}
		if (g_i % 10000 == 0)
		{
			std::cout << (((float) g_i * 100) / vecHexCellSize) << "% is done"
					<< std::endl;
		}
	}
#ifdef _WIN32
	return 0;
#else
	return NULL;
#endif
}


#ifdef _WIN32
DWORD WINAPI ProcTriangle_ref(LPVOID lpParameter)
#else
void* ProcTriangle_ref(void* lpParameter)
#endif
{
    const ThreadStructTriangle& p = *(ThreadStructTriangle*) lpParameter;
    glm::vec3 ray_direction = dirz;
    if (p.dir == X_AXIS)
    {
        ray_direction = dirx;
    }
    else if (p.dir == Y_AXIS)
    {
        ray_direction = diry;
    }

    while (g_iterHexahedron < g_iterHexahedronEnd)
    {
#ifdef _WIN32
        WaitForSingleObject(hMutex, INFINITE);
#else
        pthread_mutex_lock(&lock);
#endif
        const std::vector<Cell>::iterator iterHexahedron = g_iterHexahedron;
        Cell& hex = *iterHexahedron;
        ++g_iterHexahedron;
        g_i++;
#ifdef _WIN32
        ReleaseMutex(hMutex);
#else
        pthread_mutex_unlock(&lock);
#endif
        if (hex.empty())
            continue;
        for (int i = 0; i < 8; i++)
        {
            const Vertex& v = (p.vecPolycubeHexVertex.at(hex.at(i)));
            const glm::vec3 origin(v.x, v.y, v.z);

            if (!IsPointInsidePolyhedron(p.vecPolycubeTetVertex,
                    p.polycubeTetSurface, origin, ray_direction))
            {
                hex.push_back(INVALID_NUM);
                break;
            }
        }

        const Vertex& v2 = (p.vecPolycubeHexVertex.at(hex.at(2)));
        const Vertex& v4 = (p.vecPolycubeHexVertex.at(hex.at(4)));
        //const glm::vec3 origin(v2.x + p.delta/2, v2.y + p.delta/2, v2.z + p.delta/2);
        const glm::vec3 origin(0.5 * (v2.x + v4.x), 0.5 * (v2.y + v4.y), 0.5 * (v2.z + v4.z));

        if (!IsPointInsidePolyhedron(p.vecPolycubeTetVertex,
                p.polycubeTetSurface, origin, ray_direction))
        {
            hex.push_back(INVALID_NUM);
        }

        if (g_i % 10000 == 0)
        {
            std::cout << (((float) g_i * 100) / vecHexCellSize) << "% is done"
                    << std::endl;
        }
    }
#ifdef _WIN32
    return 0;
#else
    return NULL;
#endif
}


void Mesh::removeInvalidHexahedronsInBox(const float delta, const DIRECTION dir)
{
	const std::vector<Cell>::iterator iterHexBegin =
			m_vecPolycubeHexCell.begin();
	const std::vector<Cell>::iterator iterHexEnd = m_vecPolycubeHexCell.end();
	g_iterHexahedronBegin = iterHexBegin;
	g_iterHexahedronEnd = iterHexEnd;
	g_i = 0;
	vecHexCellSize = m_vecPolycubeHexCell.size();
	g_iterHexahedron = iterHexBegin;
#ifdef _WIN32
	hMutex = CreateMutex(NULL, false, NULL);
#else
	if (pthread_mutex_init(&lock, NULL) != 0)
	{
		std::cerr << "\n mutex init failed\n" << std::endl;
	}
#endif

	ThreadStructTriangle p(m_vecPolycubeHexVertex, V, surface, delta, dir);
#ifdef _WIN32
	HANDLE hThread[32];
	for (int i = 0; i < THREAD_NUM; i++)
	{
		hThread[i] = CreateThread(NULL, 0, ProcTriangle, &p, 0, NULL);
	}
	for (int i = 0; i < THREAD_NUM; i++)
	{
		WaitForSingleObject (hThread[i], INFINITE);
		CloseHandle (hThread[i]);
	}
#else
	pthread_t tid[32];
	for (unsigned int i = 0; i < THREAD_NUM; i++)
	{
		/*int err = */pthread_create(&tid[i], NULL, &ProcTriangle, &p);
	}
	for (unsigned int i = 0; i < THREAD_NUM; i++)
	{
		pthread_join(tid[i], NULL);
	}
	pthread_mutex_destroy(&lock);
#endif
}

void Mesh::removeInvalidHexahedronsInBox_ref(std::vector<Vertex>& Vs, const float delta,
        const DIRECTION dir)
{
    const std::vector<Cell>::iterator iterHexBegin =
            m_vecPolycubeHexCell.begin();
    const std::vector<Cell>::iterator iterHexEnd = m_vecPolycubeHexCell.end();
    g_iterHexahedronBegin = iterHexBegin;
    g_iterHexahedronEnd = iterHexEnd;
    g_i = 0;
    vecHexCellSize = m_vecPolycubeHexCell.size();
    g_iterHexahedron = iterHexBegin;
#ifdef _WIN32
    hMutex = CreateMutex(NULL, false, NULL);
#else
    if (pthread_mutex_init(&lock, NULL) != 0)
    {
        std::cerr << "\n mutex init failed\n" << std::endl;
    }
#endif

    ThreadStructTriangle p(m_vecPolycubeHexVertex, Vs, surface, delta, dir);
#ifdef _WIN32
    HANDLE hThread[32];
    for (int i = 0; i < THREAD_NUM; i++)
    {
        hThread[i] = CreateThread(NULL, 0, ProcTriangle_ref, &p, 0, NULL);
    }
    for (int i = 0; i < THREAD_NUM; i++)
    {
        WaitForSingleObject (hThread[i], INFINITE);
        CloseHandle (hThread[i]);
    }
#else
    pthread_t tid[32];
    for (unsigned int i = 0; i < THREAD_NUM; i++)
    {
        /*int err = */pthread_create(&tid[i], NULL, &ProcTriangle_ref, &p);
    }
    for (unsigned int i = 0; i < THREAD_NUM; i++)
    {
        pthread_join(tid[i], NULL);
    }
    pthread_mutex_destroy(&lock);
#endif
}

void Mesh::removeInvalidHexahedronsInBox_new(const float delta,
		const DIRECTION dir)
{
	const std::vector<Cell>::iterator iterHexBegin =
			m_vecPolycubeHexCell.begin();
	const std::vector<Cell>::iterator iterHexEnd = m_vecPolycubeHexCell.end();
	g_iterHexahedronBegin = iterHexBegin;
	g_iterHexahedronEnd = iterHexEnd;
	g_i = 0;
	vecHexCellSize = m_vecPolycubeHexCell.size();
	g_iterHexahedron = iterHexBegin;
#ifdef _WIN32
	hMutex = CreateMutex(NULL, false, NULL);
#else
	if (pthread_mutex_init(&lock, NULL) != 0)
	{
		std::cerr << "\n mutex init failed\n" << std::endl;
	}
#endif

	ThreadStructTriangle p(m_vecPolycubeHexVertex, V, surface, delta, dir);
#ifdef _WIN32
	HANDLE hThread[32];
	for (int i = 0; i < THREAD_NUM; i++)
	{
		hThread[i] = CreateThread(NULL, 0, ProcTriangle_new, &p, 0, NULL);
	}
	for (int i = 0; i < THREAD_NUM; i++)
	{
		WaitForSingleObject (hThread[i], INFINITE);
		CloseHandle (hThread[i]);
	}
#else
	pthread_t tid[32];
	for (unsigned int i = 0; i < THREAD_NUM; i++)
	{
		/*int err = */pthread_create(&tid[i], NULL, &ProcTriangle_new, &p);
	}
	for (unsigned int i = 0; i < THREAD_NUM; i++)
	{
		pthread_join(tid[i], NULL);
	}
	pthread_mutex_destroy(&lock);
#endif
}

void Mesh::removeInvalidHexahedronsInBox_fast(const float delta, const DIRECTION dir)
{
	const std::vector<Cell>::iterator iterHexBegin = m_vecPolycubeHexCell.begin();
	const std::vector<Cell>::iterator iterHexEnd = m_vecPolycubeHexCell.end();
	g_iterHexahedronBegin = iterHexBegin;
	g_iterHexahedronEnd = iterHexEnd;
	g_i = 0;
	vecHexCellSize = m_vecPolycubeHexCell.size();
	g_iterHexahedron = iterHexBegin;
#ifdef _WIN32
	hMutex = CreateMutex(NULL, false, NULL);
#else
	if (pthread_mutex_init(&lock, NULL) != 0)
	{
		std::cerr << "\n mutex init failed\n" << std::endl;
	}
#endif

	ThreadRemoveOutsidePointsParameters p(*this, delta);
#ifdef _WIN32
	HANDLE hThread[32];
	for (int i = 0; i < THREAD_NUM; i++)
	{
		hThread[i] = CreateThread(NULL, 0, ThreadRemoveOutsidePoints, &p, 0, NULL);
	}
	for (int i = 0; i < THREAD_NUM; i++)
	{
		WaitForSingleObject (hThread[i], INFINITE);
		CloseHandle (hThread[i]);
	}
#else
	pthread_t tid[32];
	for (unsigned int i = 0; i < THREAD_NUM; i++)
	{
		/*int err = */pthread_create(&tid[i], NULL, &ThreadRemoveOutsidePoints,
				&p);
	}
	for (unsigned int i = 0; i < THREAD_NUM; i++)
	{
		pthread_join(tid[i], NULL);
	}
	pthread_mutex_destroy(&lock);
#endif
}

void Mesh::removeInvalidHexahedronsInBox_Magnifier(const float delta,
		const DIRECTION dir)
{
	std::vector<unsigned long> remainCube;
	for (int i = 0; i < V.size(); i++)
	{
		remainCube.push_back(GetCubeIndex(V[i], delta));
	}
	std::sort(remainCube.begin(), remainCube.end());
	std::vector<unsigned long>::iterator iter = std::unique(remainCube.begin(),
			remainCube.end());
	remainCube.resize(std::distance(remainCube.begin(), iter));

	for (unsigned long i = 0; i < m_vecPolycubeHexCell.size(); i++)
	{
		if (std::find(remainCube.begin(), remainCube.end(), i)
				!= remainCube.end())
		{
			continue;
		}
		else
		{
			Cell& hex = m_vecPolycubeHexCell[i];
			hex.push_back(INVALID_NUM);
		}
	}
}

void Mesh::WriteHexahedralmesh(const char* outputFilename)
{
	MeshFileWriter meshFileWriter(m_vecPolycubeHexVertex, m_vecPolycubeHexCell,	outputFilename, HEXAHEDRA);
	meshFileWriter.WriteFile();
	std::cout << "poly hex is generated" << std::endl;
}
float __plus = 0.001f;
void Mesh::GenerateHexahedralMesh_plus(const Mesh& orgTet,
		const char* outputPolycubeHexFilename, const char* outputHexFilename,
		const float delta, const DIRECTION dir)
{
	GenerateSmallCubes_plus(delta, NULL, __plus);
	removeInvalidHexahedronsInBox(delta, dir);
	WriteHexahedralmesh(outputPolycubeHexFilename);

	std::vector<Vertex> vecHexVertex;
	std::vector<Cell> vecHexCell;
	CreateHexMesh(orgTet, vecHexVertex, vecHexCell);
	MeshFileWriter meshFileWriter(vecHexVertex, vecHexCell, outputHexFilename,	HEXAHEDRA);
	meshFileWriter.WriteFile();
}

void Mesh::GenerateHexahedralMesh_align(const Mesh& orgTet,
		const char* outputPolycubeHexFilename, const char* outputHexFilename,
		const float delta, const DIRECTION dir)
{
	GenerateSmallCubes_plus(delta, NULL, __plus);
	removeInvalidHexahedronsInBox(delta, dir);
	WriteHexahedralmesh(outputPolycubeHexFilename);

	std::vector<Vertex> vecHexVertex;
	std::vector<Cell> vecHexCell;
	CreateHexMesh(orgTet, vecHexVertex, vecHexCell);
	MeshFileWriter meshFileWriter(vecHexVertex, vecHexCell, outputHexFilename,	HEXAHEDRA);
	meshFileWriter.WriteFile();
}

void Mesh::GenerateHexahedralMesh_fast(const Mesh& orgTet,
		const char* outputPolycubeHexFilename, const char* outputHexFilename,
		const float delta, const DIRECTION dir)
{
	GenerateSmallCubes_plus(delta, NULL, __plus);
	removeInvalidHexahedronsInBox_fast(delta, dir);
	WriteHexahedralmesh(outputPolycubeHexFilename);

	std::vector<Vertex> vecHexVertex;
	std::vector<Cell> vecHexCell;
	CreateHexMesh(orgTet, vecHexVertex, vecHexCell);
	MeshFileWriter meshFileWriter(vecHexVertex, vecHexCell, outputHexFilename,	HEXAHEDRA);
	meshFileWriter.WriteFile();
}

void Mesh::GenerateHexahedralMesh_Magnifier(const Mesh& orgTet,
		const char* outputPolycubeHexFilename, const char* outputHexFilename,
		const float delta, const DIRECTION dir)
{
	GenerateSmallCubes_Magifier(delta);
	removeInvalidHexahedronsInBox_Magnifier(delta, dir);
	WriteHexahedralmesh(outputPolycubeHexFilename);

	std::vector<Vertex> vecHexVertex;
	std::vector<Cell> vecHexCell;
	CreateHexMesh(orgTet, vecHexVertex, vecHexCell);
	MeshFileWriter meshFileWriter(vecHexVertex, vecHexCell, outputHexFilename,	HEXAHEDRA);
	meshFileWriter.WriteFile();
}

void Mesh::GenerateHexahedralMesh(const Mesh& orgTet,
		const char* outputPolycubeHexFilename, const char* outputHexFilename,
		const float delta, const DIRECTION dir)
{
	GenerateSmallCubes(delta, outputPolycubeHexFilename);
	removeInvalidHexahedronsInBox(delta, dir);
	WriteHexahedralmesh(outputPolycubeHexFilename);

	std::vector<Vertex> vecHexVertex;
	std::vector<Cell> vecHexCell;
	CreateHexMesh(orgTet, vecHexVertex, vecHexCell);
	MeshFileWriter meshFileWriter(vecHexVertex, vecHexCell, outputHexFilename, HEXAHEDRA);
	meshFileWriter.WriteFile();
}

unsigned long Mesh::GetCubeIndex(const Vertex& v, const float delta) const
{
	const Vertex originPoint(-0.5, -0.5, -0.5);
	Vertex vv;
	vv.x = v.x - originPoint.x;
	vv.y = v.y - originPoint.y;
	vv.z = v.z - originPoint.z;
	const Vector axis(1.0, 1.0, 1.0);
	const float xInterval = delta;
	const float yInterval = delta;
	const float zInterval = delta;

	unsigned int prevPointNumber = 0;
	unsigned int PointNumber = prevPointNumber;

	unsigned int xTotalStep = ceil(fabs(axis.i / xInterval)) + 1;
	unsigned int yTotalStep = ceil(fabs(axis.j / yInterval)) + 1;
	unsigned int zTotalStep = ceil(fabs(axis.k / zInterval)) + 1;
	unsigned int CurStep = (unsigned int) (fabs(axis.k / zInterval) + 0.5 + 1);

	unsigned int yzCubeNumber = (yTotalStep - 1) * (zTotalStep - 1);

	unsigned int xCurStep = (unsigned int) (fabs(vv.x / xInterval));
	unsigned int yCurStep = (unsigned int) (fabs(vv.y / yInterval));
	unsigned int zCurStep = (unsigned int) (fabs(vv.z / zInterval));

	return (xCurStep * yzCubeNumber + yCurStep * (zTotalStep - 1) + zCurStep);
}

unsigned int Mesh::DivideHexahedronIntoSmallerOnes(const Vertex& originPoint,
		const Vector& axis, const float& xInterval, const float& yInterval,
		const float& zInterval, std::vector<Vertex>& vecVertex,	std::vector<Cell>& vecHexahedron)
{
	unsigned int prevPointNumber = vecVertex.size();
	unsigned int PointNumber = prevPointNumber;

	unsigned int xTotalStep = ceil(fabs(axis.i / xInterval)) + 1;
	unsigned int yTotalStep = ceil(fabs(axis.j / yInterval)) + 1;
	unsigned int zTotalStep = ceil(fabs(axis.k / zInterval)) + 1;
	unsigned int yzPlaneNumber = yTotalStep * zTotalStep;

	for (unsigned int xStep = 0; xStep < xTotalStep; xStep++)
	{
		for (unsigned int yStep = 0; yStep < yTotalStep; yStep++)
		{
			for (unsigned int zStep = 0; zStep < zTotalStep; zStep++)
			{
				Vertex v(originPoint.x + xStep * xInterval,
						originPoint.y + yStep * yInterval,
						originPoint.z + zStep * zInterval);
				vecVertex.push_back(v);
				if (xStep > 0 && yStep > 0 && zStep > 0)
				{
					Cell hex;
					hex.push_back(PointNumber - yzPlaneNumber);
					hex.push_back(PointNumber - yzPlaneNumber - zTotalStep);
					hex.push_back(PointNumber - yzPlaneNumber - zTotalStep - 1);
					hex.push_back(PointNumber - yzPlaneNumber - 1);
					hex.push_back(PointNumber);
					hex.push_back(PointNumber - zTotalStep);
					hex.push_back(PointNumber - zTotalStep - 1);
					hex.push_back(PointNumber - 1);
					vecHexahedron.push_back(hex);
				}

				PointNumber++;
			}
		}
	}

	return vecVertex.size();
}

void Mesh::SmoothVolume(const unsigned int time)
{
	//GetVertexInfo();
	int count = time;
	while (count-- != 0)
	{
		std::map<unsigned long, glm::vec3> map_point_meanCoordinates;
		std::multimap<unsigned long, unsigned long>::iterator iterBegin =
				VI_CI.begin();
		std::multimap<unsigned long, unsigned long>::iterator iterEnd =
				VI_CI.end();
		for (std::multimap<unsigned long, unsigned long>::iterator iter =
				iterBegin; iter != iterEnd;/* ++iter*/)
		{
//			if (iter->first == 25905)
//			{
//				std::cout << " volume 25905 connected cell centroid :" << std::endl;
//			}
			if (!V[iter->first].vinfo.bSurface
					&& V[iter->first].vinfo.bBoundaryCell)
			{
				const unsigned long vertexIndex = iter->first;
				int pointCount = 0;
				glm::vec3 mean(0, 0, 0);
				std::pair<std::multimap<unsigned long, unsigned long>::iterator,
						std::multimap<unsigned long, unsigned long>::iterator> ret =
						VI_CI.equal_range(iter->first);
				for (iter = ret.first; iter != ret.second; ++iter)
				{
					const Cell& c = C.at(iter->second);
					glm::vec3 m(0.0f, 0.0f, 0.0f);
					GetCenterPointOfCell(c, m);
					mean += m;
					pointCount++;
//					if (iter->first == 25905)
//					{
//						std::cout << " volume 25905 connected cell centroid :" << m.x << m.y << m.z << std::endl;
//					}
				}
				mean.x = mean.x / pointCount;
				mean.y = mean.y / pointCount;
				mean.z = mean.z / pointCount;
//				if (iter->first == 25905)
//				{
//					std::cout << " volume 25905 mean :" << mean.x << mean.y << mean.z << std::endl;
//				}
				map_point_meanCoordinates.insert(
						std::pair<unsigned long, glm::vec3>(
								vertexIndex/*iter->first*/, mean));
			}
			else
			{
				++iter;
			}
		}

		std::map<unsigned long, glm::vec3>::iterator iterB =
				map_point_meanCoordinates.begin();
		std::map<unsigned long, glm::vec3>::iterator iterE =
				map_point_meanCoordinates.end();
		for (std::map<unsigned long, glm::vec3>::iterator iter = iterB;
				iter != iterE; ++iter)
		{
			if (!V[iter->first].vinfo.bSurface)
			{
				Vertex& v = V[iter->first];
				v.x = iter->second.x;
				v.y = iter->second.y;
				v.z = iter->second.z;
			}
		}
	}
}

//void Mesh::SmoothVolume(const unsigned int time)
//{
//	//GetVertexInfo();
//	int count = time;
//	while(count-- != 0)
//	{
//		std::map<unsigned long, glm::vec3> map_point_meanCoordinates;
//		std::multimap<unsigned long, Cell>::iterator iterBegin = V_C.begin();
//		std::multimap<unsigned long, Cell>::iterator iterEnd = V_C.end();
//		for (std::multimap<unsigned long, Cell>::iterator iter = iterBegin; iter != iterEnd;/* ++iter*/)
//		{
//			if (iter->first == 25905)
//			{
//				std::cout << " volume 25905 connected cell centroid :" << std::endl;
//			}
//			if (!vertexInfo[iter->first].bSurface && vertexInfo[iter->first].bBoundaryCell)
//			{
//				int pointCount = 0;
//				glm::vec3 mean(0, 0, 0);
//				std::pair<std::multimap<unsigned long, Cell>::iterator, std::multimap<unsigned long, Cell>::iterator> ret = V_C.equal_range(iter->first);
//				for (iter = ret.first; iter != ret.second; ++iter)
//				{
//					const Cell& c = iter->second;
//					glm::vec3 m(0.0f, 0.0f, 0.0f);
//					GetCenterPointOfCell(c, m);
//					mean += m;
//					pointCount++;
//					if (iter->first == 25905)
//					{
//						std::cout << " volume 25905 connected cell centroid :" << m.x << m.y << m.z << std::endl;
//					}
//				}
//				mean.x = mean.x/pointCount;
//				mean.y = mean.y/pointCount;
//				mean.z = mean.z/pointCount;
//				if (iter->first == 25905)
//				{
//					std::cout << " volume 25905 mean :" << mean.x << mean.y << mean.z << std::endl;
//				}
//				map_point_meanCoordinates.insert(std::pair<unsigned long, glm::vec3>(iter->first, mean));
//			}
//			else
//			{
//				++iter;
//			}
//		}
//
//		std::map<unsigned long, glm::vec3>::iterator iterB = map_point_meanCoordinates.begin();
//		std::map<unsigned long, glm::vec3>::iterator iterE = map_point_meanCoordinates.end();
//		for (std::map<unsigned long, glm::vec3>::iterator iter = iterB; iter != iterE; ++iter)
//		{
//			if (!vertexInfo[iter->first].bSurface){
//			Vertex& v = V[iter->first];
//			v.x = iter->second.x;
//			v.y = iter->second.y;
//			v.z = iter->second.z;}
//		}
//	}
//}

void Mesh::SmoothSurface(const unsigned int time)
{
	//GetVertexInfo();
	int count = time;
	while (count-- != 0)
	{
		std::map<unsigned long, glm::vec3> map_point_meanCoordinates;
		std::multimap<unsigned long, SFace>::iterator iterBegin = V_F.begin();
		std::multimap<unsigned long, SFace>::iterator iterEnd = V_F.end();
		for (std::multimap<unsigned long, SFace>::iterator iter = iterBegin;
				iter != iterEnd;/* ++iter*/)
		{
			if (V[iter->first].vinfo.bSurface
				&& !V[iter->first].vinfo.bBoundaryLine
				&& V[iter->first].vinfo.bNeedSmoothing)
			{
				int pointCount = 0;
				const unsigned long vertexIndex = iter->first;
				glm::vec3 mean(0, 0, 0);
				std::pair<std::multimap<unsigned long, SFace>::iterator,
						std::multimap<unsigned long, SFace>::iterator> ret =
						V_F.equal_range(iter->first);
				for (iter = ret.first; iter != ret.second; ++iter)
				{
					const Face& f = iter->second;
					if (IsFaceOnSurface(f))
					{
						glm::vec3 m(0.0f, 0.0f, 0.0f);
						GetCenterPointOfCell(f, m);
						mean += m;
						pointCount++;
					}
				}
				mean.x = mean.x / pointCount;
				mean.y = mean.y / pointCount;
				mean.z = mean.z / pointCount;
				map_point_meanCoordinates.insert(
						std::pair<unsigned long, glm::vec3>(vertexIndex, mean));
			}
			else
			{
				++iter;
			}
		}

		std::map<unsigned long, glm::vec3>::iterator iterB =
				map_point_meanCoordinates.begin();
		std::map<unsigned long, glm::vec3>::iterator iterE =
				map_point_meanCoordinates.end();
		for (std::map<unsigned long, glm::vec3>::iterator iter = iterB;
				iter != iterE; ++iter)
		{
			Vertex& v = V[iter->first];
			v.x = iter->second.x;
			v.y = iter->second.y;
			v.z = iter->second.z;
		}
	}
}

void Mesh::GetParameterAndCubeIndicesOfVertexIndices()
{
	tets.clear();
	paras.clear();
	paras_flag.clear();
	vec_cubeIndex.clear();
	std::cout << "cubesize = " << CUBESIZE << std::endl;
	for (int i = 0; i < V.size()/*vertices.size()*/; i++)
	{
		const Vertex& v = V.at(i/*vertices.at(i)*/);
		// 1. judge v is in which cube
		const unsigned long cubeIndex = GetCubeIndex(v, CUBESIZE);
		vec_cubeIndex.push_back(cubeIndex);
		const Cell& cubeCell = m_vecPolycubeHexCell.at(cubeIndex);
		bool flag = false;
		for (int k = 0; k < 5; k++)
		{
			Cell tet;
			for (int j = 0; j < 4; j++)
			{
				tet.push_back(cubeCell.at(hexTet[k][j]));
			}
			if (IsVertexInsideTetrahedron(v, tet, m_vecPolycubeHexVertex))
			{
				tet.resize(4);
				tets.push_back(tet);
				const Vertex& p0 = m_vecPolycubeHexVertex[tet.at(0)];
				const Vertex& p1 = m_vecPolycubeHexVertex[tet.at(1)];
				const Vertex& p2 = m_vecPolycubeHexVertex[tet.at(2)];
				const Vertex& p3 = m_vecPolycubeHexVertex[tet.at(3)];

				glm::vec3 v03((p0.x - p3.x), (p0.y - p3.y), (p0.z - p3.z));
				glm::vec3 v13((p1.x - p3.x), (p1.y - p3.y), (p1.z - p3.z));
				glm::vec3 v23((p2.x - p3.x), (p2.y - p3.y), (p2.z - p3.z));

				glm::mat3x3 T(v03, v13, v23);
				glm::vec3 r(v.x, v.y, v.z);
				glm::vec3 r4(p3.x, p3.y, p3.z);
				glm::mat3x3 T_inverse = glm::inverse(T);
				glm::vec3 rambda = T_inverse * (r - r4);
				paras.push_back(rambda);

				//std::cout << "###" << hexTet[k][0] << " " << hexTet[k][1] << " " << hexTet[k][2] << " " << hexTet[k][3] << " " << std::endl;
				//std::cout << tet[0] << " " << tet[1] << " " << tet[2] << " " << tet[3] << " " << rambda.x << " " << rambda.y << " " << rambda.z << std::endl;
				flag = true;
				break;
			}
		}
		if (!flag)
		{
			Cell tet;
			for (int j = 0; j < 4; j++)
			{
				tet.push_back(cubeCell.at(hexTet[0][j]));
			}
			tets.push_back(tet);
			glm::vec3 rambda(0.25, 0.25, 0.25);
			paras.push_back(rambda);

			std::cout << "i = " << i << " Bad!" << std::endl;
			paras_flag.push_back(false);
		}
		else
			paras_flag.push_back(true);
		//vec_w.push_back(C);
	}
}

void Mesh::MapbackFromDeformedCubes()
{
	//parameterization of new OrigTet to new Cube Mesh
	std::cout << "parameterization of new OrigTet to new Cube Mesh"
			<< std::endl;
	for (int i = 0; i < V.size()/*vertices.size()*/; i++)
	{
		if (!paras_flag[i])
			continue;
		Vertex& v = V.at(i/*vertices.at(i)*/);
		// 1. judge v is in which cube
		const unsigned long cubeIndex = vec_cubeIndex.at(i);
		const Cell& cubeCell = m_vecPolycubeHexCell.at(cubeIndex);
		const Cell& tet = tets.at(i);

		const Vertex& origT0 = m_vecPolycubeHexVertex[tet.at(0)];
		const Vertex& origT1 = m_vecPolycubeHexVertex[tet.at(1)];
		const Vertex& origT2 = m_vecPolycubeHexVertex[tet.at(2)];
		const Vertex& origT3 = m_vecPolycubeHexVertex[tet.at(3)];

		glm::vec3 p03((origT0.x - origT3.x), (origT0.y - origT3.y),
				(origT0.z - origT3.z));
		glm::vec3 p13((origT1.x - origT3.x), (origT1.y - origT3.y),
				(origT1.z - origT3.z));
		glm::vec3 p23((origT2.x - origT3.x), (origT2.y - origT3.y),
				(origT2.z - origT3.z));
		glm::mat3x3 OrigT(p03, p13, p23);
		glm::vec3 new_r4(origT3.x, origT3.y, origT3.z);
		const glm::vec3& rambda = paras.at(i);
		glm::vec3 new_r = OrigT * rambda;
		new_r += new_r4;
		Vertex baryCenter(new_r.x, new_r.y, new_r.z);
		v.x = baryCenter.x;
		v.y = baryCenter.y;
		v.z = baryCenter.z;
		//std::cout << tet[0] << " " << tet[1] << " " << tet[2] << " " << tet[3] << " " << rambda.x << " " << rambda.y << " " << rambda.z << std::endl;
	}
}

void Mesh::GetDensityField()
{
	// 1. Compute the average length, denoted by e, of each hex-element of the hex-mesh (there are totally 12 edges for each element)
	// 2. Assign 1/e to the corresponding hex-element. Compute this for all the elements, which will give rise to a scalar field
	vec_densityFiled.clear();
	vec_averageLength_V.clear();
//	std::multimap<unsigned long, Cell>::iterator iterBegin = V_C.begin();
//	std::multimap<unsigned long, Cell>::iterator iterEnd = V_C.end();
	std::multimap<unsigned long, unsigned long>::iterator iterBegin =
			VI_CI.begin();
	std::multimap<unsigned long, unsigned long>::iterator iterEnd = VI_CI.end();
	float max_density = -10000000.0f;
	float min_density = 10000000.0f;
//	for (std::multimap<unsigned long, Cell>::iterator iter = iterBegin; iter != iterEnd;/*++iter*/)
	for (std::multimap<unsigned long, unsigned long>::iterator iter = iterBegin;
			iter != iterEnd;/*++iter*/)
	{
		std::vector<float> edgeLength;
//		std::pair<std::multimap<unsigned long, Cell>::iterator, std::multimap<unsigned long, Cell>::iterator> ret =
//				V_C.equal_range(iter->first);
		std::pair<std::multimap<unsigned long, unsigned long>::iterator,
				std::multimap<unsigned long, unsigned long>::iterator> ret =
				VI_CI.equal_range(iter->first);
		for (iter = ret.first; iter != ret.second; ++iter)
		{
			//const Cell& c = iter->second;
			const Cell& c = C.at(iter->second);
			if (c.size() >= 8)
			{
				for (int i = 0; i < 12; i++)
				{
					unsigned long e0 = c[HexEdge[i][0]];
					unsigned long e1 = c[HexEdge[i][1]];
					glm::vec3 v0(V[e0].x, V[e0].y, V[e0].z);
					glm::vec3 v1(V[e1].x, V[e1].y, V[e1].z);
					edgeLength.push_back(glm::length(v0 - v1));
				}
			}
			else if (c.size() == 4)
			{
				for (int i = 0; i < 4; i++)
				{
					unsigned long e0 = c[TetEdge[i][0]];
					unsigned long e1 = c[TetEdge[i][1]];
					glm::vec3 v0(V[e0].x, V[e0].y, V[e0].z);
					glm::vec3 v1(V[e1].x, V[e1].y, V[e1].z);
					edgeLength.push_back(glm::length(v0 - v1));
				}
			}
			else if (c.size() == 3)
			{
				for (int i = 0; i < 3; i++)
				{
					unsigned long e0 = c[TriEdge[i][0]];
					unsigned long e1 = c[TriEdge[i][1]];
					glm::vec3 v0(V[e0].x, V[e0].y, V[e0].z);
					glm::vec3 v1(V[e1].x, V[e1].y, V[e1].z);
					edgeLength.push_back(glm::length(v0 - v1));
				}
			}
		}
		float density = 0.0;
		for (int i = 0; i < edgeLength.size(); i++)
		{
			density += edgeLength[i];
		}

		//std::cout << "density = " << density << std::endl;
		float averageLength = density / edgeLength.size();
		vec_averageLength_V.push_back(averageLength);
		density = 1.0 / averageLength;
		max_density = max_density > density ? max_density : density;
		min_density = min_density < density ? min_density : density;
		vec_densityFiled.push_back(density);
	}

	float density_category = 1.0 / (max_density - min_density);
	for (int i = 0; i < vec_densityFiled.size(); i++)
	{
		vec_densityFiled.at(i) = (vec_densityFiled.at(i) - min_density)
				* density_category;
	}
}

void Mesh::GetSurfaceDensityField()
{
	// 1. Compute the average length, denoted by e, of each hex-element of the hex-mesh (there are totally 12 edges for each element)
	// 2. Assign 1/e to the corresponding hex-element. Compute this for all the elements, which will give rise to a scalar field
	vec_densityFiled.clear();
	vec_averageLength_V.clear();
	std::multimap<unsigned long, SFace>::iterator iterBegin = V_F.begin();
	std::multimap<unsigned long, SFace>::iterator iterEnd = V_F.end();
	float max_density = -10000000.0f;
	float min_density = 10000000.0f;
	for (std::multimap<unsigned long, SFace>::iterator iter = iterBegin;
			iter != iterEnd;/*++iter*/)
	{
		if (!V[iter->first].vinfo.bSurface)
		{
			++iter;
			vec_densityFiled.push_back(0);
			continue;
		}
		std::vector<float> edgeLength;
		std::pair<std::multimap<unsigned long, SFace>::iterator,
				std::multimap<unsigned long, SFace>::iterator> ret =
				V_F.equal_range(iter->first);
		for (iter = ret.first; iter != ret.second; ++iter)
		{
			const Cell& c = iter->second;
			if (!IsFaceOnSurface(c))
			{
				continue;
			}
			if (c.size() >= 8)
			{
				for (int i = 0; i < 12; i++)
				{
					unsigned long e0 = c[HexEdge[i][0]];
					unsigned long e1 = c[HexEdge[i][1]];
					glm::vec3 v0(V[e0].x, V[e0].y, V[e0].z);
					glm::vec3 v1(V[e1].x, V[e1].y, V[e1].z);
					edgeLength.push_back(glm::length(v0 - v1));
				}
			}
			else if (c.size() == 4)
			{
				for (int i = 0; i < 4; i++)
				{
					unsigned long e0 = c[TetEdge[i][0]];
					unsigned long e1 = c[TetEdge[i][1]];
					glm::vec3 v0(V[e0].x, V[e0].y, V[e0].z);
					glm::vec3 v1(V[e1].x, V[e1].y, V[e1].z);
					edgeLength.push_back(glm::length(v0 - v1));
				}
			}
			else if (c.size() == 3)
			{
				for (int i = 0; i < 3; i++)
				{
					unsigned long e0 = c[TriEdge[i][0]];
					unsigned long e1 = c[TriEdge[i][1]];
					glm::vec3 v0(V[e0].x, V[e0].y, V[e0].z);
					glm::vec3 v1(V[e1].x, V[e1].y, V[e1].z);
					edgeLength.push_back(glm::length(v0 - v1));
				}
			}
		}
		float density = 0.0;
		for (int i = 0; i < edgeLength.size(); i++)
		{
			density += edgeLength[i];
		}

		//std::cout << "density = " << density << std::endl;
		float averageLength = density / edgeLength.size();
		vec_averageLength_V.push_back(averageLength);
		density = 1.0 / averageLength;
		max_density = max_density > density ? max_density : density;
		min_density = min_density < density ? min_density : density;
		vec_densityFiled.push_back(density);
	}

	float density_category = 1.0 / (max_density - min_density);
	for (int i = 0; i < vec_densityFiled.size(); i++)
	{
		float normal = (vec_densityFiled.at(i) - min_density)
				* density_category;
		if (normal >= 0 && normal <= 1)
		{
			vec_densityFiled.at(i) = normal;
		}
		else
		{
			vec_densityFiled.at(i) = 0;
		}
	}
}

static float _GetScaledJacobian(const glm::vec3& i, const glm::vec3& j,
		const glm::vec3& k)
{
	const glm::mat3x3 m(i, j, k);
	return glm::determinant(m);
}
template<typename T>
void Swap(T& t1, T& t2)
{
	T t = t1;
	t1 = t2;
	t2 = t;
}

float cal_volume_Tet_real(float v0[3],float v1[3],float v2[3],float v3[3])
{
	float v1v0[3],v2v0[3],v3v0[3];
	for(int i=0;i<3;i++)
	{
		v1v0[i]=v1[i]-v0[i];
		v2v0[i]=v2[i]-v0[i];
		v3v0[i]=v3[i]-v0[i];
	}

	float norm1=sqrt(v1v0[0]*v1v0[0]+v1v0[1]*v1v0[1]+v1v0[2]*v1v0[2]);
	float norm2=sqrt(v2v0[0]*v2v0[0]+v2v0[1]*v2v0[1]+v2v0[2]*v2v0[2]);
	float norm3=sqrt(v3v0[0]*v3v0[0]+v3v0[1]*v3v0[1]+v3v0[2]*v3v0[2]);

	float volume=v1v0[0]*(v2v0[1]*v3v0[2]-v2v0[2]*v3v0[1])-v1v0[1]*(v2v0[0]*v3v0[2]-v2v0[2]*v3v0[0])+v1v0[2]*(v2v0[0]*v3v0[1]-v2v0[1]*v3v0[0]);
	return volume;
}
const int V_T[8][4] =
{
    {0, 3, 4, 1},
	{1, 0, 5, 2},
	{2, 1, 6, 3},
	{3, 2, 7, 0},
	{4, 7, 5, 0},
	{5, 4, 6, 1},
	{6, 5, 7, 2},
	{7, 6, 4, 3}
};

bool Mesh::JudgeDirection(const Cell& c)
{
	float v[8][3];
	for (int i = 0; i < c.size(); i++)
	{
		v[i][0] = V.at(c.at(i)).x;
		v[i][1] = V.at(c.at(i)).y;
		v[i][2] = V.at(c.at(i)).z;
	}

	float VL[8];
	for (int i = 0; i < 8; i++)
	{
		const int* p = V_T[i];
		VL[i] = cal_volume_Tet_real(v[p[0]], v[p[1]], v[p[2]], v[p[3]]);
	}

	if (VL[0] + VL[1] + VL[2] + VL[3] + VL[4] + VL[5] + VL[6] + VL[7] < 0)
		return false;
	else
		return true;
}

const float Mesh::GetScaledJacobian(const Cell& c)
{
//	assert(c.size() >= 8);
//	unsigned long ccc[8] = {cc[4], cc[5], cc[6], cc[7], cc[0], cc[1], cc[2], cc[3]};
//	Cell c;
//	for (int i = 0; i < 8; i++)
//		c.push_back(ccc[i]);

	Cell c1(c);
	if (!JudgeDirection(c))
	{
		swap(c1[0], c1[3]);
		swap(c1[1], c1[2]);
		swap(c1[4], c1[7]);
		swap(c1[5], c1[6]);
	}
	float minScaledJacobian = 1;
	for (int n = 0; n < 7; n++)
	{
		const Vertex& o = V.at(c1.at(n));

		const Vertex& i = V.at(c1.at(HexPoint_Points[n][0]));
		const Vertex& j = V.at(c1.at(HexPoint_Points[n][1]));
		const Vertex& k = V.at(c1.at(HexPoint_Points[n][2]));

		const glm::vec3 ei(i.x - o.x, i.y - o.y, i.z - o.z);
		const glm::vec3 ej(j.x - o.x, j.y - o.y, j.z - o.z);
		const glm::vec3 ek(k.x - o.x, k.y - o.y, k.z - o.z);

		const float length_i = glm::length(ei);
		const float length_j = glm::length(ej);
		const float length_k = glm::length(ek);

		const glm::vec3 ni(ei.x / length_i, ei.y / length_i, ei.z / length_i);
		const glm::vec3 nj(ej.x / length_j, ej.y / length_j, ej.z / length_j);
		const glm::vec3 nk(ek.x / length_k, ek.y / length_k, ek.z / length_k);

		float scaledJacobian = _GetScaledJacobian(ni, nj, nk);
		minScaledJacobian =
				minScaledJacobian < scaledJacobian ?
						minScaledJacobian : scaledJacobian;

//		float scaled_jacobian=ei[0]*(ej[1]*ek[2]-ej[2]*ek[1])-ei[1]*(ej[0]*ek[2]-ej[2]*ek[0])+ei[2]*(ej[0]*ek[1]-ej[1]*ek[0]);
//		scaled_jacobian/=length_i*length_j*length_k;
//		minScaledJacobian = minScaledJacobian < scaled_jacobian ? minScaledJacobian : scaled_jacobian;
	}

	return minScaledJacobian;
}

void WriteInvertedHexFile(const Cell& c, const char* filename = "inverted1.vtk")
{
	std::ofstream fileStream(filename, std::ios_base::out);
	fileStream << "OFF" << std::endl;

}
void Mesh::OutputScaledJacobianDataFile(const std::vector<Cell>& C,
		const char* filename)
{
	invertedCellIndex.clear();
	lowQualityCellIndex.clear();
	float minScaledJacobian = 1;
	float maxScaledJacobian = -1;
	double sumScaledJacobian = 0;
	int invertedCellCount = 0;
	int historgramBin[5] =
	{ 0, 0, 0, 0, 0 }; // 0~0.2, 0.2~0.4, 0.4~0.6, 0.6~0.8, 0.8~1.0
	std::ofstream fileStream(filename, std::ios_base::out);
	for (int i = 0; i < C.size(); i++)
	{
		float scaledJacobian = GetScaledJacobian(C.at(i));
		fileStream << scaledJacobian << std::endl;
		minScaledJacobian =
				minScaledJacobian < scaledJacobian ?
						minScaledJacobian : scaledJacobian;
		maxScaledJacobian =
				maxScaledJacobian > scaledJacobian ?
						maxScaledJacobian : scaledJacobian;
		sumScaledJacobian += scaledJacobian;
		if (scaledJacobian < 0)
		{
			invertedCellCount++;
			invertedCellIndex.push_back(i);
			/////////////////////////////////////
			const Cell& c = C.at(i);
			char file_name[64];
			sprintf(file_name, "inverted%d.off", invertedCellCount);
			std::ofstream file_stream(file_name, std::ios_base::out);
			file_stream << "OFF" << std::endl;
			file_stream << "8 6 0" << std::endl;
			for (int j = 0; j < 8; j++)
				file_stream << V[c[j]].x << " " << V[c[j]].y << " " << V[c[j]].z
						<< std::endl;
			for (int j = 0; j < 6; j++)
				file_stream << 4 << " " << HexFaces[j][0] << " "
						<< HexFaces[j][1] << " " << HexFaces[j][2] << " "
						<< HexFaces[j][3] << std::endl;
			file_stream.close();
			/////////////////////////////////////
		}
		else
		{
			int binIndex = int(scaledJacobian / 0.2);
			historgramBin[binIndex]++;
			if (binIndex == 0)
			{
				lowQualityCellIndex.push_back(i);
				/////////////////////////////////////
				const Cell& c = C.at(i);
				char file_name[64];
				sprintf(file_name, "lowQuality%d.off", historgramBin[binIndex]);
				std::ofstream file_stream(file_name, std::ios_base::out);
				file_stream << "OFF" << std::endl;
				file_stream << "8 6 0" << std::endl;
				for (int j = 0; j < 8; j++)
					file_stream << V[c[j]].x << " " << V[c[j]].y << " "
							<< V[c[j]].z << std::endl;
				for (int j = 0; j < 6; j++)
					file_stream << 4 << " " << HexFaces[j][0] << " "
							<< HexFaces[j][1] << " " << HexFaces[j][2] << " "
							<< HexFaces[j][3] << std::endl;
				file_stream.close();
				/////////////////////////////////////
			}
		}
	}

	fileStream << "min = " << minScaledJacobian << std::endl;
	fileStream << "avg = " << sumScaledJacobian / C.size() << std::endl;
	fileStream << "max = " << maxScaledJacobian << std::endl;
	fileStream << std::endl;
	fileStream << "#inverted cells = " << invertedCellCount << std::endl;
	fileStream << std::endl;
	fileStream << "[0.0, 0.2) = "
			<< ((float) (historgramBin[0] * 100)) / C.size() << "%"
			<< std::endl;
	fileStream << "[0.2, 0.4) = "
			<< ((float) (historgramBin[1] * 100)) / C.size() << "%"
			<< std::endl;
	fileStream << "[0.4, 0.6) = "
			<< ((float) (historgramBin[2] * 100)) / C.size() << "%"
			<< std::endl;
	fileStream << "[0.6, 0.8) = "
			<< ((float) (historgramBin[3] * 100)) / C.size() << "%"
			<< std::endl;
	fileStream << "[0.8, 1.0) = "
			<< ((float) (historgramBin[4] * 100)) / C.size() << "%"
			<< std::endl;
	fileStream.close();

	std::cout << " min = " << minScaledJacobian;
	std::cout << " avg = " << sumScaledJacobian / C.size();
	std::cout << " max = " << maxScaledJacobian;
	std::cout << std::endl;
	std::cout << "#inverted cells = " << invertedCellCount << std::endl;
	std::cout << "[0.0, 0.2) = "
			<< ((float) (historgramBin[0] * 100)) / C.size() << "%";
	std::cout << "[0.2, 0.4) = "
			<< ((float) (historgramBin[1] * 100)) / C.size() << "%";
	std::cout << "[0.4, 0.6) = "
			<< ((float) (historgramBin[2] * 100)) / C.size() << "%";
	std::cout << "[0.6, 0.8) = "
			<< ((float) (historgramBin[3] * 100)) / C.size() << "%";
	std::cout << "[0.8, 1.0) = "
			<< ((float) (historgramBin[4] * 100)) / C.size() << "%"
			<< std::endl;
}

glm::vec3 Mesh::GetFaceNormal(const Face& f) const
{
	const Vertex* pVertex0 = &V[f.at(0)];
	const Vertex* pVertex1 = &V[f.at(1)];
	const Vertex* pVertex2 = &V[f.at(2)];

	const glm::vec3 v0(pVertex0->x, pVertex0->y, pVertex0->z);
	const glm::vec3 v1(pVertex1->x, pVertex1->y, pVertex1->z);
	const glm::vec3 v2(pVertex2->x, pVertex2->y, pVertex2->z);

	const glm::vec3 v10 = v0 - v1;
	const glm::vec3 v12 = v2 - v1;
	glm::vec3 normal = glm::cross(v12, v10);
	return normal;
}
void Mesh::GetNormalOfSurfaceFaces()
{
	N_F.clear();
	const std::vector<Face>::const_iterator iterSurfaceBegin = surface.begin();
	const std::vector<Face>::const_iterator iterSurfaceEnd = surface.end();
	for (std::vector<Face>::const_iterator iterSurface = iterSurfaceBegin;
			iterSurface != iterSurfaceEnd; ++iterSurface)
	{
		const Vertex* pVertex0 = &V[iterSurface->at(0)];
		const Vertex* pVertex1 = &V[iterSurface->at(1)];
		const Vertex* pVertex2 = &V[iterSurface->at(2)];

		const glm::vec3 v0(pVertex0->x, pVertex0->y, pVertex0->z);
		const glm::vec3 v1(pVertex1->x, pVertex1->y, pVertex1->z);
		const glm::vec3 v2(pVertex2->x, pVertex2->y, pVertex2->z);

		const glm::vec3 v10 = v0 - v1;
		const glm::vec3 v12 = v2 - v1;
		const glm::vec3 normal = glm::cross(v12, v10);
		N_F.push_back(normal);
	}
	N_F.resize(N_F.size());
}

void Mesh::GetNormalOfSurfaceVertices()
{
	N_V.clear();
	//N_V.resize(V.size());
	glm::vec3 n(0.0f, 0.0f, 0.0f);
	std::multimap<unsigned long, SFace>::iterator iterBegin = V_F.begin();
	std::multimap<unsigned long, SFace>::iterator iterEnd = V_F.end();
	for (std::multimap<unsigned long, SFace>::iterator iter = iterBegin;
			iter != iterEnd;/* ++iter*/)
	{
		if (V[iter->first].vinfo.bSurface)
		{
			int faceCount = 0;
			glm::vec3 normal_v(0.0f, 0.0f, 0.0f);
			const unsigned long vertexIndex = iter->first;
			std::pair<std::multimap<unsigned long, SFace>::iterator,
					std::multimap<unsigned long, SFace>::iterator> ret =
					V_F.equal_range(iter->first);
			for (iter = ret.first; iter != ret.second; ++iter)
			{
				const Face& f = iter->second;
				if (IsFaceOnSurface(f))
				{
					const Vertex* pVertex0 = &V[f.at(0)];
					const Vertex* pVertex1 = &V[f.at(1)];
					const Vertex* pVertex2 = &V[f.at(2)];

					const glm::vec3 v0(pVertex0->x, pVertex0->y, pVertex0->z);
					const glm::vec3 v1(pVertex1->x, pVertex1->y, pVertex1->z);
					const glm::vec3 v2(pVertex2->x, pVertex2->y, pVertex2->z);

					const glm::vec3 v10 = v0 - v1;
					const glm::vec3 v12 = v2 - v1;
					const glm::vec3 normal = glm::cross(v12, v10);
					faceCount++;
					normal_v += normal;
				}
			}
//			normal_v.x = normal_v.x/faceCount;
//			normal_v.y = normal_v.y/faceCount;
//			normal_v.z = normal_v.z/faceCount;
			float normal_len = glm::length(normal_v);
			normal_v.x = normal_v.x / normal_len;
			normal_v.y = normal_v.y / normal_len;
			normal_v.z = normal_v.z / normal_len;
			N_V.push_back(normal_v);
			//N_V.at(vertexIndex) = normal_v;
		}
		else
		{
			N_V.push_back(n);
			++iter;
		}
	}
	N_V.resize(N_V.size());
}

void Mesh::GetSurfaceVertexIndices()
{
	if (!surfaceVertexIndices.empty())
		return;
	for (int i = 0; i < surface.size(); i++)
	{
		const Cell& cell = surface.at(i);
		for (int i = 0; i < cell.size(); i++)
		{
			surfaceVertexIndices.push_back(cell.at(i));
		}
	}
	std::sort(surfaceVertexIndices.begin(), surfaceVertexIndices.end());
	std::vector<unsigned long>::iterator iter = std::unique(
			surfaceVertexIndices.begin(), surfaceVertexIndices.end());
	surfaceVertexIndices.resize(
			std::distance(surfaceVertexIndices.begin(), iter));
}

void Mesh::VTKExtractSurface()
{

}
#include <omp.h>
bool Mesh::IsAllVerticesInsideMesh(const Mesh& mesh)
{
	if (V.size() == 0)
	{
		std::cout << "ERROR testing IsAllVerticesInsideMesh" << std::endl;
	}
#pragma omp for
	for (int i = 0; i < V.size(); i++)
	{
		const Vertex& v = V.at(i);
		const glm::vec3 p(v.x, v.y, v.z);
		if (!IsPointInsidePolyhedron(mesh.V, mesh.surface, p, dirz))
		{
//			std::cout << "Vertex " << i << " is not inside the mesh" << std::endl;
			outsidePoints.push_back(i);
		}
	}

	std::cout << "all " << V.size() << " Vertices are Inside Mesh. Test DONE" << std::endl;
}

void Mesh::FixMesh()
{
	if (m_vecPolycubeHexCell.empty() || m_vecPolycubeHexCell.at(0).size() < 8)
		return;

	std::vector<unsigned long> v_real_index;
	size_t c_size = 0;
	for (unsigned long i = 0; i < m_vecPolycubeHexCell.size(); i++)
	{
		if (m_vecPolycubeHexCell.at(i).size() > 8)
			continue;
		else
		{
			for (int j = 0; j < 8; j++)
			{
				v_real_index.push_back(m_vecPolycubeHexCell.at(i).at(j));
			}
			c_size++;
		}
	}

	std::sort(v_real_index.begin(), v_real_index.end());
	std::vector<unsigned long>::iterator iter = std::unique(
			v_real_index.begin(), v_real_index.end());
	v_real_index.resize(std::distance(v_real_index.begin(), iter));

	std::vector<Vertex> V_(v_real_index.size());
	for (unsigned long i = 0; i < v_real_index.size(); i++)
	{
		const Vertex& v = m_vecPolycubeHexVertex.at(v_real_index.at(i));
		V_.at(i) = v;
	}
	m_vecPolycubeHexVertex = V_;
	//////////////////////////////////////////////////////
	std::map<unsigned long, unsigned long> v_v;
	unsigned long index = 0;
	for (unsigned long i = 0; i < v_real_index.size(); i++)
	{
		v_v[v_real_index.at(i)] = i;
	}

	std::vector<Cell> _C;
	Cell c(8);
	for (unsigned long i = 0; i < m_vecPolycubeHexCell.size(); i++)
	{
		if (m_vecPolycubeHexCell.at(i).size() == 8)
		{
			for (int j = 0; j < 8; j++)
			{
				const Cell& cell = m_vecPolycubeHexCell.at(i);
				c.at(j) = v_v[cell.at(j)];
			}
			_C.push_back(c);
		}
	}
	_C.resize(_C.size());
	m_vecPolycubeHexCell = _C;
}

//void Mesh::InitGeodesicDistance(const char* polydataFilename)
//{
//	reader = vtkPolyDataReader::New();
//	reader->SetFileName(polydataFilename);
//
//	normals = vtkPolyDataNormals::New();
//
//	const int geodesicMethod = 0;
//	const int interpolationOrder = 0;
//	const double distanceOffset = 0;
//
//	// We need to ensure that the dataset has normals if a distance offset was
//	// specified.
//	if (fabs(distanceOffset) > 1e-6)
//	{
//		normals->SetInputConnection(reader->GetOutputPort());
//		normals->SplittingOff();
//
//		// vtkPolygonalSurfacePointPlacer needs cell normals
//		// vtkPolygonalSurfaceContourLineInterpolator needs vertex normals
//		normals->ComputeCellNormalsOn();
//		normals->ComputePointNormalsOn();
//		normals->Update();
//	}
//
//	vtkPolyData *pd =
//			(fabs(distanceOffset) > 1e-6) ?
//					normals->GetOutput() : reader->GetOutput();
//
//	mapper = vtkPolyDataMapper::New();
//	mapper->SetInputConnection(
//			fabs(distanceOffset) > 1e-6 ?
//					normals->GetOutputPort() : reader->GetOutputPort());
//
//	actor = vtkActor::New();
//	actor->SetMapper(mapper);
//
//	ren = vtkRenderer::New();
//	renWin = vtkRenderWindow::New();
//	renWin->AddRenderer(ren);
//	iren = vtkRenderWindowInteractor::New();
//	iren->SetRenderWindow(renWin);
//
//	ren->AddActor(actor);
//	ren->GetActiveCamera()->SetPosition(-3.68, .447, 1.676);
//	ren->GetActiveCamera()->Roll(150);
//	ren->ResetCamera();
//	ren->ResetCameraClippingRange();
//
//	ren->AddActor(actor);
//
//	// Create a contour widget to interactively trace the path
//
//	contourWidget = vtkContourWidget::New();
//	contourWidget->SetInteractor(iren);
//	vtkOrientedGlyphContourRepresentation *rep =
//			vtkOrientedGlyphContourRepresentation::SafeDownCast(
//					contourWidget->GetRepresentation());
//	rep->GetLinesProperty()->SetColor(1, 0.2, 0);
//	rep->GetLinesProperty()->SetLineWidth(5.0);
//
//	pointPlacer = vtkPolygonalSurfacePointPlacer::New();
//	pointPlacer->AddProp(actor);
//	pointPlacer->GetPolys()->AddItem(pd);
//	rep->SetPointPlacer(pointPlacer);
//
//	// Snap the contour nodes to the closest vertices on the mesh
//	pointPlacer->SnapToClosestPointOn();
//
//	interpolator = vtkPolygonalSurfaceContourLineInterpolator2::New();
//	interpolator->GetPolys()->AddItem(pd);
//	interpolator->SetGeodesicMethod(geodesicMethod);
//	interpolator->SetInterpolationOrder(interpolationOrder);
//	rep->SetLineInterpolator(interpolator);
//	if (fabs(distanceOffset) > 1e-6)
//	{
//		pointPlacer->SetDistanceOffset(distanceOffset);
//		interpolator->SetDistanceOffset(distanceOffset);
//	}
//
////	renWin->Render();
////	iren->Initialize();
//	contourWidget->EnabledOn();
//	//interpolator->InterpolateLine(reader->GetOutput(), 14483, 14574);
//}

//double Mesh::ComputeGeodesicDistance(vtkPolyData *polyData,
//		vtkIdType beginVertId, vtkIdType endVertId)
//{
//	std::vector<vtkIdType> ids;
//	interpolator->InterpolateLine(reader->GetOutput(), beginVertId, endVertId,
//			ids);
//
//	double geodesicDistance = 0.0;
//	int lastIndex = ids.at(0);
//	for (int i = 1; i < ids.size(); i++)
//	{
//		int curIndex = ids.at(i);
//		geodesicDistance += V.at(curIndex).Distance(V.at(lastIndex));
//		lastIndex = curIndex;
//	}
//	return geodesicDistance;
//}

bool IsPairExisted(const std::multimap<DirectedEdge, DirectedEdge>& E_E, const std::pair<DirectedEdge, DirectedEdge>& pair)
{
	const std::pair<std::multimap<DirectedEdge, DirectedEdge>::const_iterator, std::multimap<DirectedEdge, DirectedEdge>::const_iterator> ret =
			E_E.equal_range(pair.first);

	const std::multimap<DirectedEdge, DirectedEdge>::const_iterator iterLowerBound = ret.first;
	const std::multimap<DirectedEdge, DirectedEdge>::const_iterator iterUpperBound = ret.second;
	if (iterLowerBound == E_E.end() || iterUpperBound == E_E.end())
		return false;
	std::multimap<DirectedEdge, DirectedEdge>::const_iterator iterBound = iterLowerBound;

	// if the e_e is not existed, insert it to E_E
	for (; iterBound != iterUpperBound; ++iterBound)
	{
		const DirectedEdge& e = iterBound->second;
		if (e == pair.second)
		{
			return true;
			//break;
		}
	}

	return false;
}

void Mesh::GetConeConfiguration()
{
	if (C.empty())
		return;
//	const std::vector<Cell>::const_iterator iterCellBegin = C.begin();
//	const std::vector<Cell>::const_iterator iterCellEnd = C.end();
//	std::vector<Cell>::const_iterator iterCell = iterCellBegin;
//
//	if (m_meshType == MESH_TYPE_HEXAHEDRON
//			|| m_meshType == MESH_TYPE_HEXAHEDRON_VTK)
//	{
//		for (iterCell = iterCellBegin; iterCell != iterCellEnd; ++iterCell)
//		{
//			for (int i = 0; i < 12; i++)
//			{
//				const Edge e(HexEdge[i][0], HexEdge[i][1]);
//				E.push_back(e);
//				const DirectedEdge de1(HexEdge[i][0], HexEdge[i][1]);
//				const DirectedEdge de2(HexEdge[i][1], HexEdge[i][0]);
//				DE.push_back(de1);
//				DE.push_back(de2);
//			}
//		}
//	}
//	std::sort(E.begin(), E.end());
//	std::vector<Edge>::iterator iterE = std::unique(E.begin(), E.end());
//	E.resize(std::distance(E.begin(), iterE));
//
//	std::sort(DE.begin(), DE.end());
//	std::vector<DirectedEdge>::iterator iterDE = std::unique(DE.begin(), DE.end());
//	DE.resize(std::distance(DE.begin(), iterDE));
//
//	std::copy(F.begin(), F.end(), SF.begin());
//	std::sort(SF.begin(), SF.end());
//	std::vector<SFace>::iterator iterSF = std::unique(SF.begin(), SF.end());
//	SF.resize(std::distance(SF.begin(), iterSF));
//
//	for (unsigned long i = 0; i < SF.size(); i++)
//	{
//		const SFace& sf = SF.at(i);
//		const size_t s = sf.size();
//		for (unsigned int j = 0; j < s; j++)
//		{
//			const Edge e(sf.at(j % s), sf.at((j + 1) % s));
//			const std::pair<Edge, unsigned long> e_f(e, i);
//			if (!IsPairExisted(E_F, e_f))
//			{
//				E_F.insert(e_f);
//			}
//		}
//	}
	for (unsigned long i = 0; i < C.size(); i++)
	{
		const Cell& c = C.at(i);
		for (int j = 0; j < 12; j++)
		{
			const DirectedEdge de1(c.at(HexTrippleEdge[j][0]), c.at(HexTrippleEdge[j][1]));
			const DirectedEdge de2(c.at(HexTrippleEdge[j][1]), c.at(HexTrippleEdge[j][0]));
			const DirectedEdge de3(c.at(HexTrippleEdge[j][2]), c.at(HexTrippleEdge[j][3]), i);
			const DirectedEdge de4(c.at(HexTrippleEdge[j][5]), c.at(HexTrippleEdge[j][4]), i);

			const std::pair<DirectedEdge, DirectedEdge> de_de1(de1, de3);
			const std::pair<DirectedEdge, DirectedEdge> de_de2(de2, de4);

			// if (!IsPairExisted(E_E, de_de1));
				E_E.insert(de_de1);
			// if (!IsPairExisted(E_E, de_de2));
				E_E.insert(de_de2);
		}
	}

}

double Mesh::GetSingleEcone(const DirectedEdge& e_ij, const DirectedEdge& u_k_kplus)
{
	const int i = e_ij.p0;
	const int j = e_ij.p1;
	const int k = u_k_kplus.p0;
	const int k1 = u_k_kplus.p1;
	const glm::vec3 vi(V[i].x, V[i].y, V[i].z);
	const glm::vec3 vj(V[j].x, V[j].y, V[j].z);
	const glm::vec3 uk(V[k].x, V[k].y, V[k].z);
	const glm::vec3 uk1(V[k1].x, V[k1].y, V[k1].z);

	glm::vec3 nk = -glm::cross(uk - vi, uk1 - vi);
	const float len_k = glm::length(nk);
	nk.x /= len_k;
	nk.y /= len_k;
	nk.z /= len_k;

	glm::vec3 eij = vj - vi;
	const float len = glm::length(eij);
	eij.x /= len;
	eij.y /= len;
	eij.z /= len;

	const double energy = (eij.x - nk.x)*(eij.x - nk.x) + (eij.y - nk.y)*(eij.y - nk.y) + (eij.z - nk.z)*(eij.z - nk.z);
	return energy;
}

double Mesh::GetSingleEcone(const DirectedEdge& e_ij, const DirectedEdge& u_k_kplus, const glm::vec3& nk_c)
{
	const int i = e_ij.p0;
	const int j = e_ij.p1;
	const int k = u_k_kplus.p0;
	const int k1 = u_k_kplus.p1;
	const glm::vec3 vi(V[i].x, V[i].y, V[i].z);
	const glm::vec3 vj(V[j].x, V[j].y, V[j].z);
	const glm::vec3 uk(V[k].x, V[k].y, V[k].z);
	const glm::vec3 uk1(V[k1].x, V[k1].y, V[k1].z);

	glm::vec3 nk = -glm::cross(uk - vi, uk1 - vi);
	const float len_k = glm::length(nk);
	nk.x /= len_k;
	nk.y /= len_k;
	nk.z /= len_k;

	const double energy = (nk_c.x - nk.x)*(nk_c.x - nk.x) + (nk_c.y - nk.y)*(nk_c.y - nk.y) + (nk_c.z - nk.z)*(nk_c.z - nk.z);
	return energy;
}

void Mesh::GetLocalMinimalConeVertex(const std::multimap<DirectedEdge, DirectedEdge>::const_iterator it, glm::vec3& opt_v)
{
	double min_sum = 100000000.0;

	const std::pair<std::multimap<DirectedEdge, DirectedEdge>::const_iterator, std::multimap<DirectedEdge, DirectedEdge>::const_iterator> bound =
			E_E.equal_range(it->first);
	for (std::multimap<DirectedEdge, DirectedEdge>::const_iterator iter_out = bound.first; iter_out != bound.second; ++iter_out)
	{
		const DirectedEdge& e_ij = iter_out->first;
		const DirectedEdge& u_k_kplus = iter_out->second;

		const int i = e_ij.p0;
		const int j = e_ij.p1;
		const int k = u_k_kplus.p0;
		const int k1 = u_k_kplus.p1;
		const glm::vec3 vi(V[i].x, V[i].y, V[i].z);
		const glm::vec3 vj(V[j].x, V[j].y, V[j].z);
		const glm::vec3 uk(V[k].x, V[k].y, V[k].z);
		const glm::vec3 uk1(V[k1].x, V[k1].y, V[k1].z);

		glm::vec3 nk = -glm::cross(uk - vi, uk1 - vi);
		const float len_k = glm::length(nk);
		nk.x /= len_k;
		nk.y /= len_k;
		nk.z /= len_k;

		double sum = 0.0;
		for (std::multimap<DirectedEdge, DirectedEdge>::const_iterator iter_inner = bound.first; iter_inner != bound.second; ++iter_inner)
		{
			double a = GetSingleEcone(iter_inner->first, iter_inner->second, nk);
			sum += a;
		}
		if (sum < min_sum)
		{
			glm::vec3 eij = vj - vi;
			const float len = glm::length(eij);
			opt_v.x = vi.x + len*nk.x;
			opt_v.y = vi.y + len*nk.y;
			opt_v.z = vi.z + len*nk.z;
			//////////////////////////////////////
			min_sum = sum;
		}
	}
}

void Mesh::WriteConeShape(const char* filename, const std::multimap<DirectedEdge, DirectedEdge>::const_iterator it)
{
	std::vector<Cell> coneCell;
	const std::pair<std::multimap<DirectedEdge, DirectedEdge>::const_iterator, std::multimap<DirectedEdge, DirectedEdge>::const_iterator> bound =
		E_E.equal_range(it->first);
	for (std::multimap<DirectedEdge, DirectedEdge>::const_iterator iter_out = bound.first; iter_out != bound.second; ++iter_out)
	{
		const Cell& c = C.at(iter_out->second.cellIndex);
		coneCell.push_back(c);
	}

	Mesh mesh(V, coneCell, HEXAHEDRA);
	MeshFileWriter meshWriter(mesh, filename);
	meshWriter.WriteFile();
}
void Mesh::OptimizeConeShape(int iterMaxNum)
{
	const std::multimap<DirectedEdge, DirectedEdge>::const_iterator iterBegin = E_E.begin();
	const std::multimap<DirectedEdge, DirectedEdge>::const_iterator iterEnd = E_E.end();
	for (std::multimap<DirectedEdge, DirectedEdge>::const_iterator iter = iterBegin; iter != iterEnd;/* ++iter*/)
	{
		const DirectedEdge& eij = iter->first;
		if (!V.at(eij.p0).vinfo.bSurface/* && !vertexInfo.at(eij.p1).bSurface*/)
		{
			Vertex& target_v = V.at(eij.p1);
			glm::vec3 opt_v;
			GetLocalMinimalConeVertex(iter, opt_v);
			static int i = 0;
			////////////////////////
			if (i < 10)
			{
			char s[100];
			sprintf(s, "cone_org%d.vtk", ++i);
			std::string str((const char*)s);
			WriteConeShape(str.c_str(), iter);

			target_v.x = opt_v.x;
			target_v.y = opt_v.y;
			target_v.z = opt_v.z;

			sprintf(s, "cone_opt%d.vtk", i);
			std::string str1((const char*)s);
			WriteConeShape(str1.c_str(), iter);
			}
		}

		const std::pair<std::multimap<DirectedEdge, DirectedEdge>::const_iterator, std::multimap<DirectedEdge, DirectedEdge>::const_iterator> bound =
					E_E.equal_range(iter->first);

		iter = bound.second;
	}
}

bool Mesh::VerifyStructure(const int genus)
{
	//ExtractSurface();
	if (!C.empty())
	{
		const std::vector<Cell>::const_iterator iterCellBegin = C.begin();
		const std::vector<Cell>::const_iterator iterCellEnd = C.end();
		std::vector<Cell>::const_iterator iterCell = iterCellBegin;

		if (m_cellType == HEXAHEDRA || m_cellType == TETRAHEDRA)
		{
			for (iterCell = iterCellBegin; iterCell != iterCellEnd; ++iterCell)
			{
				for (int i = 0; i < 12; i++)
				{
					const Edge e(iterCell->at(HexEdge[i][0]), iterCell->at(HexEdge[i][1]));
					E.push_back(e);
					const DirectedEdge de1(iterCell->at(HexEdge[i][0]), iterCell->at(HexEdge[i][1]));
					const DirectedEdge de2(iterCell->at(HexEdge[i][1]), iterCell->at(HexEdge[i][0]));
					DE.push_back(de1);
					DE.push_back(de2);
				}
			}
		}
	}

//	std::sort(E.begin(), E.end());
//	std::vector<Edge>::iterator iterE = std::unique(E.begin(), E.end());
//	E.resize(std::distance(E.begin(), iterE));
	if (uniqueE.empty())
	for (int i = 0; i < E.size(); i++)
	{
		const Edge& e = E.at(i);
		std::vector<Edge>::iterator iterE = std::find(uniqueE.begin(), uniqueE.end(), e);
		if (iterE == uniqueE.end())
		{
			uniqueE.push_back(e);
		}
	}

//	std::sort(DE.begin(), DE.end());
//	std::vector<DirectedEdge>::iterator iterDE = std::unique(DE.begin(), DE.end());
//	DE.resize(std::distance(DE.begin(), iterDE));

//	std::sort(SF.begin(), SF.end());
//	std::vector<SFace>::iterator iterSF = std::unique(SF.begin(), SF.end());
//	SF.resize(std::distance(SF.begin(), iterSF));

	if (uniqueF.empty())
	{
		SF.resize(F.size());
		std::copy(F.begin(), F.end(), SF.begin());
		for (int i = 0; i < SF.size(); i++)
		{
			const SFace& f = SF.at(i);
			std::vector<SFace>::iterator iterF = std::find(uniqueF.begin(), uniqueF.end(), f);
			if (iterF == uniqueF.end())
			{
				uniqueF.push_back(f);
			}
		}
	}

	long l = V.size() - uniqueE.size() + uniqueF.size() - C.size();
	if (l == 1 - genus)
	{
		std::cout << "#V:" << V.size() << " - #E:" << uniqueE.size() << " + #F:" << uniqueF.size() << " - #C:" << C.size() << " = " << l << std::endl;
		return true;
	}
	else
	{
		std::cout << "#V:" << V.size() << " - #E:" << uniqueE.size() << " + #F:" << uniqueF.size() << " - #C:" << C.size() << " = " << left << std::endl;
		return false;
	}
}

int Mesh::GetGenus()
{
	if (!C.empty())
	{
		const std::vector<Cell>::const_iterator iterCellBegin = C.begin();
		const std::vector<Cell>::const_iterator iterCellEnd = C.end();
		std::vector<Cell>::const_iterator iterCell = iterCellBegin;

		if (m_meshType == MESH_TYPE_HEXAHEDRON || m_meshType == MESH_TYPE_HEXAHEDRON_VTK)
		{
			for (iterCell = iterCellBegin; iterCell != iterCellEnd; ++iterCell)
			{
				for (int i = 0; i < 12; i++)
				{
					/*const*/ Edge e(iterCell->at(HexEdge[i][0]), iterCell->at(HexEdge[i][1]));
					if (e.p0 > e.p1)
						swap(e.p0, e.p1);
					E.push_back(e);
					const DirectedEdge de1(iterCell->at(HexEdge[i][0]), iterCell->at(HexEdge[i][1]));
					const DirectedEdge de2(iterCell->at(HexEdge[i][1]), iterCell->at(HexEdge[i][0]));
					DE.push_back(de1);
					DE.push_back(de2);
				}
			}
		}
		else if (m_meshType == MESH_TYPE_TETRAHEDRON_VTK)
		{
			for (iterCell = iterCellBegin; iterCell != iterCellEnd; ++iterCell)
			{
				for (int i = 0; i < 6; i++)
				{
					/*const*/ Edge e(iterCell->at(TetEdge[i][0]), iterCell->at(TetEdge[i][1]));
					if (e.p0 > e.p1)
						swap(e.p0, e.p1);
					E.push_back(e);
//					const DirectedEdge de1(iterCell->at(TetEdge[i][0]), iterCell->at(TetEdge[i][1]));
//					const DirectedEdge de2(iterCell->at(TetEdge[i][1]), iterCell->at(TetEdge[i][0]));
//					DE.push_back(de1);
//					DE.push_back(de2);
				}
			}
		}
	}

//	if (uniqueE.empty())
//	for (int i = 0; i < E.size(); i++)
//	{
//		const Edge& e = E.at(i);
//		std::vector<Edge>::iterator iterE = std::find(uniqueE.begin(), uniqueE.end(), e);
//		if (iterE == uniqueE.end())
//		{
//			uniqueE.push_back(e);
//		}
//	}
	if (uniqueE.empty())
	{
		uniqueE = E;
		std::sort(uniqueE.begin(), uniqueE.end());
		std::vector<Edge>::iterator iterE = std::unique(uniqueE.begin(), uniqueE.end());
		uniqueE.resize(std::distance(uniqueE.begin(), iterE));
	}

	if (uniqueF.empty())
	{
		for (int i = 0; i < F.size(); i++)
		{
			Face f = F.at(i);
			std::sort(f.begin(), f.end());
			uniqueF.push_back(f);
		}
		uniqueF.resize(uniqueF.size());
//		SF.resize(F.size());
//		std::copy(F.begin(), F.end(), SF.begin());
//		for (int i = 0; i < SF.size(); i++)
//		{
//			const SFace& f = SF.at(i);
//			std::vector<SFace>::iterator iterF = std::find(uniqueF.begin(), uniqueF.end(), f);
//			if (iterF == uniqueF.end())
//			{
//				uniqueF.push_back(f);
//			}
//		}
		std::sort(uniqueF.begin(), uniqueF.end());
		std::vector<SFace>::iterator iterSF = std::unique(uniqueF.begin(), uniqueF.end());
		uniqueF.resize(std::distance(uniqueF.begin(), iterSF));
	}

	// V - E + F - C = 1 - g
	// g = 1 - (V - E + F - C);
	m_genus = 1 - (V.size() - uniqueE.size() + uniqueF.size() - C.size());
	std::cout << "g = 1 - (V - E + F - C) = " << m_genus << std::endl;
	std::cout << "#V:" << V.size() << " - #E:" << uniqueE.size() << " + #F:" << uniqueF.size() << " - #C:" << C.size() << " = " << std::endl;
	if (m_genus >= 0)
		return m_genus;
	else
	{
		std::cout << "g = 1 - (V - E + F - C)" << std::endl;
		std::cout << "#V:" << V.size() << " - #E:" << uniqueE.size() << " + #F:" << uniqueF.size() << " - #C:" << C.size() << " = " << left << std::endl;
		return false;
	}
}

bool Mesh::VerifyElements() const
{
	bool bRet = true;
	const std::multimap<unsigned long, unsigned long>::const_iterator iterBegin = V_V.begin();
	const std::multimap<unsigned long, unsigned long>::const_iterator iterEnd = V_V.end();
	for (std::multimap<unsigned long, unsigned long>::const_iterator iter = iterBegin; iter != iterEnd; )
	{
		const std::pair<std::multimap<unsigned long, unsigned long>::const_iterator, std::multimap<unsigned long, unsigned long>::const_iterator>
			iterRet = V_V.equal_range(iter->first);
		size_t valence = std::distance(iterRet.first, iterRet.second);
		if (valence > 6)
		{
			std::cout << "Warning! " << valence << " vertices connect to the vertex " << iter->first << std::endl;
			bRet = false;
		}
		if (valence < 3)
		{
			std::cout << "Warning! Only " << valence << " vertices connect to the vertex " << iter->first << std::endl;
			bRet = false;
		}
		iter = iterRet.second;
	}

	return bRet;
}

void Mesh::GetCellAndNeighborCells()
{
	if (!C_C.empty())
		return;
	for (size_t i = 0; i < C.size(); i++)
	{
		const unsigned long cellId = i;
		std::vector<unsigned long> neighborCellIds;
		GetNeighborCells(cellId, neighborCellIds);
		for (size_t j = 0; j < neighborCellIds.size(); j++)
		{
			const std::pair<unsigned long, unsigned long> p(cellId, neighborCellIds.at(j));
			C_C.insert(p);
		}
	}
}

void Mesh::GetFaceAndNeighborFaces()
{
	if (!F_F.empty())
		return;
	for (size_t i = 0; i < surface.size(); i++)
	{
		const unsigned long faceId = i;
		std::vector<unsigned long> neighborFaceIds;
		GetNeighborFaces(faceId, neighborFaceIds);
		for (size_t j = 0; j < neighborFaceIds.size(); j++)
		{
			const std::pair<unsigned long, unsigned long> p(faceId, neighborFaceIds.at(j));
			F_F.insert(p);
		}
	}
}

static void Union(const Cell& c1, const Cell& c2, std::vector<unsigned long>& u)
{
	for (size_t j = 0; j < c1.size(); j++)
	{
		for (size_t k = 0; k < c2.size(); k++)
		{
			if (c1.at(j) == c2.at(k))
			{
				u.push_back(c1.at(j));
				break;
			}
		}
	}
}

void Mesh::GetNeighborCells(const unsigned long cellId, std::vector<unsigned long>& neighborCellIds)
{
	const Cell& cell = C.at(cellId);
	for (size_t i = 0; i < C.size(); i++)
	{
		if (i == cellId)
			continue;
		const unsigned long neighborCellId = i;
		const Cell& neighborCell = C.at(neighborCellId);
		std::vector<unsigned long> u;
		Union(cell, neighborCell, u);
		if (u.size() == 2)
		{
			neighborCellIds.push_back(neighborCellId);
		}
	}
}

void Mesh::GetNeighborFaces(const unsigned long faceId, std::vector<unsigned long>& neighborFaceIds)
{
	const Face& face = surface.at(faceId);
	for (size_t i = 0; i < surface.size(); i++)
	{
		if (i == faceId)
			continue;
		const unsigned long neighborFaceId = i;
		const Face& neighborFace = surface.at(neighborFaceId);
		std::vector<unsigned long> u;
		Union(face, neighborFace, u);
		if (u.size() == 2)
		{
			neighborFaceIds.push_back(neighborFaceId);
		}
	}
}

const glm::vec3 Mesh::GetNormal(unsigned long faceIndex)
{
	const Face& f = surface.at(faceIndex);
	const Vertex& v_0 = V.at(f.at(0));
	const Vertex& v_1 = V.at(f.at(1));
	const Vertex& v_2 = V.at(f.at(2));

	const glm::vec3 v0(v_0.x, v_0.y, v_0.z);
	const glm::vec3 v1(v_1.x, v_1.y, v_1.z);
	const glm::vec3 v2(v_2.x, v_2.y, v_2.z);

	const glm::vec3 v10 = v0 - v1;
	const glm::vec3 v12 = v2 - v1;
	const glm::vec3 normal = glm::cross(v12, v10);

	return normal;
}

void Mesh::GetFaceType()
{
	faceType.resize(surface.size());
	for (unsigned long i = 0; i < surface.size(); i++)
	{
		const glm::vec3 normal = GetNormal(i);
		faceType.at(i) = GeoUtil::GetFaceType(normal);
	}
}

void Mesh::SmoothPatches()
{
	for (int i = 0; i < m_patches.size(); i++)
	{
		Patch& patch = m_patches.at(i);
		patch.Smooth();
	}
}

void Mesh::SetFixedVertices(std::vector<unsigned long>& fixed_points)
{
	for (size_t i = 0; i < fixed_points.size(); i++)
		V.at(fixed_points.at(i)).vinfo.bNeedSmoothing = false;
}

bool Mesh::IsPointInPoints(unsigned long pointId, const std::vector<unsigned long>& Ids)
{
    return std::find(Ids.begin(), Ids.end(), pointId) != Ids.end();
}

bool Mesh::IsCutPoint(unsigned long pointId, const std::vector<unsigned long>& Ids)
{
	std::pair<std::multimap<unsigned long, unsigned long>::iterator,
	          std::multimap<unsigned long, unsigned long>::iterator> ret = V_V.equal_range(pointId);
	for (std::multimap<unsigned long, unsigned long>::iterator iter = ret.first; iter != ret.second; ++iter)
	{
		unsigned long connectedVertexIndex = iter->second;
		if (!IsPointInPoints(connectedVertexIndex, Ids))
		{
			return true;
		}
	}
	return false;
}

int Mesh::GroupCutPoints(const std::vector<unsigned long>& Ids, std::vector<std::vector<unsigned long> >& groupIds)
{
	std::vector<unsigned long> ids = Ids;
	for (std::vector<unsigned long>::iterator iterId = ids.begin(); iterId != ids.end(); /*++iterId*/)
	{
		unsigned long pointId = *iterId;
		std::vector<unsigned long> group;
		group.push_back(pointId);
		std::vector<unsigned long>::iterator iterid1 = std::find(ids.begin(), ids.end(), pointId);
		ids.erase(iterid1);
		while (true)
		{
			std::pair<std::multimap<unsigned long, unsigned long>::iterator,
					  std::multimap<unsigned long, unsigned long>::iterator> ret = V_V.equal_range(pointId);
			bool found = false;
			for (std::multimap<unsigned long, unsigned long>::iterator iter = ret.first; iter != ret.second; ++iter)
			{
				unsigned long connectedVertexIndex = iter->second;
				std::vector<unsigned long>::iterator iterid = std::find(ids.begin(), ids.end(), connectedVertexIndex);
				if (iterid != ids.end())
				{
					group.push_back(connectedVertexIndex);
					ids.erase(iterid);
					pointId = connectedVertexIndex;
					found = true;
					break;
				}
			}
			if (!found)
				break;
		}
		groupIds.push_back(group);
		iterId = ids.begin();
	}

	return groupIds.size();
}
