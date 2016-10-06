/*
 * Patch.cpp
 *
 *  Created on: Jun 19, 2015
 *      Author: cotrik
 */

#include "Patch.h"
#include "PolycubeMesh.h"
Patch::Patch()
: m_pMesh(NULL)
{
	// TODO Auto-generated constructor stub

}

Patch::Patch(const Patch& patch)
: m_v(patch.m_v)
, m_e(patch.m_e)
, m_f(patch.m_f)
, m_corners(patch.m_corners)
, m_pMesh(patch.m_pMesh)
, dir(patch.dir)
//, pos(patch.pos)
, center(patch.center)
{
	// TODO Auto-generated constructor stub

}

Patch::Patch(const PolycubeMesh& mesh)
: m_pMesh((PolycubeMesh*)&mesh)
{
	// TODO Auto-generated constructor stub

}

Patch::Patch(const PolycubeMesh& mesh, const std::vector<unsigned long>& patchVertex)
: m_pMesh((PolycubeMesh*)&mesh)
, m_v(patchVertex)
{

}
Patch::Patch(const PolycubeMesh& mesh, const std::vector<unsigned long>& f, std::vector<Edge>& e, std::vector<unsigned long>& v)
: m_pMesh((PolycubeMesh*)&mesh)
, m_f(f)
, m_e(e)
, m_v(v)
{
	GetCorners();
	GetPatchPosition();
	dir = (DIRECTION)(mesh.faceType.at(f.at(0)) - 1);
}

Patch::~Patch()
{
	// TODO Auto-generated destructor stub
}

//const bool Patch::operator < (const Patch& right) const
//{
//	for (int i = 0; i < m_corners.size(); i++)
//	{
//
//	}
//}

void Patch::GetCorners()
{
	for (int i = 0; i < m_v.size(); i++)
	{
		const unsigned long index = m_v.at(i);
		if (m_pMesh->V.at(index).vinfo.bCorner)
		{
			const Vertex& v = m_pMesh->V.at(index);
			m_corners.push_back(v);
		}
	}
	std::sort(m_corners.begin(), m_corners.end());
}

//const bool operator < (const Patch& patch) const
//{
//	if (dir < patch.dir)
//	{
//		return true;
//	}
//	else if (dir == patch.dir)
//	{
//		if (dir == X_AXIS)
//		{
//			if (fabs(center.x - patch.center.x) < 5e-4)
//				if (glm::length(center - patch.center) < 3e-3)
//					return false;
//			return false;
//		}
//		else if (dir == Y_AXIS)
//		{
//			if (fabs(center.y - patch.center.y) < 5e-4)
//				if (glm::length(center - patch.center) < 3e-3)
//					return center.z < patch.center.z;
//			return false;
//		}
//		else if (dir == Z_AXIS)
//		{
//			if (fabs(center.z - patch.center.z) < 5e-4)
//				if (glm::length(center - patch.center) < 3e-3)
//					return center.z < patch.center.z;
//			return false;
//		}
//		else
//		{
//			std::cout << "Direction ERROR!" << std::endl;
//		}
//	}
//
//	return false;
//}

const bool Patch::operator < (const Patch& patch) const
{
	if (dir != patch.dir)
	{
		return dir < patch.dir;
	}
	else //if (dir == patch.dir)
	{
		if (dir == X_AXIS)
		{
			if (fabs(center.x - patch.center.x) < 5e-4)
				if (glm::length(center - patch.center) < 3e-3)	return false;
				else return center.y < patch.center.y;
			return center.x < patch.center.x;
		}
		else if (dir == Y_AXIS)
		{
			if (fabs(center.y - patch.center.y) < 5e-4)
				if (glm::length(center - patch.center) < 3e-3)	return false;
				else return center.z < patch.center.z;
			return center.y < patch.center.y;
		}
		else if (dir == Z_AXIS)
		{
			if (fabs(center.z - patch.center.z) < 5e-4)
				if (glm::length(center - patch.center) < 3e-3)	return false;
				else return center.x < patch.center.x;
			return center.z < patch.center.z;
		}
		else
		{
			std::cout << "Direction ERROR!" << std::endl;
		}
	}

	return false;
}

const bool Patch::operator == (const Patch& patch) const
{
	return ((dir == patch.dir) && glm::length(center - patch.center) < 3e-3);
}

void Patch::GetPatchPosition()
{
	glm::vec3 c(0.0, 0.0, 0.0);
	for (int i = 0; i < m_v.size(); i++)
	{
		const unsigned long index = m_v.at(i);
		c.x += m_pMesh->V.at(index).x;
		c.y += m_pMesh->V.at(index).y;
		c.z += m_pMesh->V.at(index).z;
	}

	center.x = c.x/m_v.size();
	center.y = c.y/m_v.size();
	center.z = c.z/m_v.size();
}

void Patch::ModifyPatchesPosition()
{
    for (int i = 0; i < m_v.size(); i++)
	{
		const unsigned long index = m_v.at(i);
		if (dir == X_AXIS) m_pMesh->V.at(index).x = center.x;
		else if (dir == Y_AXIS) m_pMesh->V.at(index).y = center.y;
		else if (dir == Z_AXIS) m_pMesh->V.at(index).z = center.z;
	}
}

void Patch::ModifyPatchesPosition_ROUND()
{
    for (int i = 0; i < m_v.size(); i++)
	{
		const unsigned long index = m_v.at(i);
		if (dir == X_AXIS) m_pMesh->V.at(index).x = int(center.x/0.005)*0.005;
		else if (dir == Y_AXIS) m_pMesh->V.at(index).y = int(center.y/0.005)*0.005;
		else if (dir == Z_AXIS) m_pMesh->V.at(index).z = int(center.z/0.005)*0.005;
	}
}

void Patch::ModifyPatchesPosition_ROUND(float delta)
{
    for (int i = 0; i < m_v.size(); i++)
	{
		const unsigned long index = m_v.at(i);
		if (dir == X_AXIS) m_pMesh->V.at(index).x = int(center.x/delta)*delta;
		else if (dir == Y_AXIS) m_pMesh->V.at(index).y = int(center.y/delta)*delta;
		else if (dir == Z_AXIS) m_pMesh->V.at(index).z = int(center.z/delta)*delta;
	}
}

void Patch::ModifyPatchesPosition_ROUND2(float delta)
{
    for (int i = 0; i < m_v.size(); i++)
    {
        const unsigned long index = m_v.at(i);
        if (dir == X_AXIS) m_pMesh->V.at(index).x = int(center.x/delta)*delta;
        else if (dir == Y_AXIS) m_pMesh->V.at(index).y = int(center.y/delta)*delta;
        else if (dir == Z_AXIS) m_pMesh->V.at(index).z = int(center.z/delta)*delta;
    }
}

void Patch::GetParametrization(const Patch& triPatch, TriQuadParas& p)
{
	p.paras.resize(m_v.size());
	p.triIndex.resize(m_v.size());
	p.flag.resize(m_v.size());
	for (int i = 0; i < m_v.size(); i++)
	{
		bool flag = false;
		unsigned long surfaceVertexIndex = m_v.at(i);
		const Vertex& v = m_pMesh->V.at(surfaceVertexIndex);
		const std::vector<unsigned long>& V = triPatch.m_v;
		const std::vector<unsigned long>& F = triPatch.m_f;
		for (int j = 0; j < F.size(); j++)
		{
			unsigned long fi = F.at(j);
			const Cell& tri = triPatch.m_pMesh->surface.at(fi);
			const Vertex& p0 = triPatch.m_pMesh->V.at(tri.at(0));
			const Vertex& p1 = triPatch.m_pMesh->V.at(tri.at(1));
			const Vertex& p2 = triPatch.m_pMesh->V.at(tri.at(2));
			glm::vec2 v02; //((p0.x - p2.x), (p0.y - p2.y), (p0.z - p2.z));
			glm::vec2 v12; //((p1.x - p2.x), (p1.y - p2.y), (p1.z - p2.z));
			glm::vec2 r(v.y, v.z);
			glm::vec2 r2(p2.y, p2.z);
			if (triPatch.m_pMesh->faceType.at(fi) == FACE_X)
			{
				v02 = glm::vec2((p0.y - p2.y), (p0.z - p2.z));
				v12 = glm::vec2((p1.y - p2.y), (p1.z - p2.z));
				r = glm::vec2(v.y, v.z);
				r2 = glm::vec2(p2.y, p2.z);
			}
			else if (triPatch.m_pMesh->faceType.at(fi) == FACE_Y)
			{
				v02 = glm::vec2((p0.x - p2.x), (p0.z - p2.z));
				v12 = glm::vec2((p1.x - p2.x), (p1.z - p2.z));
				r = glm::vec2(v.x, v.z);
				r2 = glm::vec2(p2.x, p2.z);
			}
			else if (triPatch.m_pMesh->faceType.at(fi) == FACE_Z)
			{
				v02 = glm::vec2((p0.x - p2.x), (p0.y - p2.y));
				v12 = glm::vec2((p1.x - p2.x), (p1.y - p2.y));
				r = glm::vec2(v.x, v.y);
				r2 = glm::vec2(p2.x, p2.y);
			}
			glm::mat2x2 T(v02, v12);
			glm::mat2x2 T_inverse = glm::inverse(T);
			glm::vec2 rambda = T_inverse * (r - r2);
            //std::cout << "rambda.x + rambda.y = " << rambda.x + rambda.y << std::endl;
			if (rambda.x > -1e-3 && rambda.y > -1e-3
					&& (rambda.x + rambda.y) < 1.001)
			{
				p.paras.at(i) = rambda;
				p.triIndex.at(i) = j;
				p.flag.at(i) = true;
				flag = true;
				break;
			}
		}
		if (!flag)
		{
			///////////////////////
			const glm::vec3 vv(v.x, v.y, v.z);
			double closeCellIndex = 0;
			double dis = 100000000;
			for (int j = 0; j < F.size(); j++)
			{
				const Cell& tet = triPatch.m_pMesh->surface.at(F.at(j));
				const Vertex& p0 = triPatch.m_pMesh->V.at(tet.at(0));
				const Vertex& p1 = triPatch.m_pMesh->V.at(tet.at(1));
				const Vertex& p2 = triPatch.m_pMesh->V.at(tet.at(2));

				const glm::vec3 v0(p0.x, p0.y, p0.z);
				const glm::vec3 v1(p1.x, p1.y, p1.z);
				const glm::vec3 v2(p2.x, p2.y, p2.z);

				const glm::vec3 d0(v0 - vv);
				const glm::vec3 d1(v1 - vv);
				const glm::vec3 d2(v2 - vv);

				double d = glm::length(d0) + glm::length(d1) + glm::length(d2);
				if (d < dis)
				{
					dis = d;
					closeCellIndex = j;
				}
			}
			const Cell& tet = triPatch.m_pMesh->surface.at(closeCellIndex);
			const Vertex& p0 = triPatch.m_pMesh->V.at(tet.at(0));
			const Vertex& p1 = triPatch.m_pMesh->V.at(tet.at(1));
			const Vertex& p2 = triPatch.m_pMesh->V.at(tet.at(2));

			const glm::vec3 v0(p0.x, p0.y, p0.z);
			const glm::vec3 v1(p1.x, p1.y, p1.z);
			const glm::vec3 v2(p2.x, p2.y, p2.z);

			const glm::vec3 d0(v0 - vv);
			const glm::vec3 d1(v1 - vv);
			const glm::vec3 d2(v2 - vv);
			const double d = glm::length(d0) + glm::length(d1)
					+ glm::length(d2);
			//rambda.x = d0/d; rambda.y = d1/d; rambda.z = d2/d;
			///////////////////////
			//glm::vec3 rambda(0.25, 0.25, 0.5);
			glm::vec2 rambda(glm::length(d0) / d, glm::length(d1) / d);
			p.paras.at(i) = rambda;
			p.triIndex.at(i) = closeCellIndex;
			p.flag.at(i) = false;

			//p.errorPointIndices.push_back(i);
			std::cout << "vertex = " << surfaceVertexIndex << " Fail to parametrization!" << std::endl;
			//p.paras_flag.push_back(false);
		}
	}
}

void Patch::MapbackToOrigTri(const std::vector<Vertex>& orgV, const Patch& triPatch, const TriQuadParas& p, std::vector<Vertex>& quadV)
{
	//parameterization of new OrigTet to new Cube Mesh
	std::cout << "MapbackToOrigTet" << std::endl;
	//quadV.resize(m_v.size());
	for (int i = 0; i < m_v.size(); i++)
	{
//		if (p.flag.at(i)){
		//unsigned long surfaceVertexIndex = quadMesh.surfaceVertexIndices.at(i);
		const unsigned long triIndex = p.triIndex.at(i);
		const Cell& tri = triPatch.m_pMesh->surface.at(triPatch.m_f.at(triIndex));

		const Vertex& origT0 = orgV[tri.at(0)];
		const Vertex& origT1 = orgV[tri.at(1)];
		const Vertex& origT2 = orgV[tri.at(2)];

		glm::vec2 p02; // ((origT0.x - origT3.x), (origT0.y - origT3.y), (origT0.z - origT3.z));
		glm::vec2 p12; // ((origT1.x - origT3.x), (origT1.y - origT3.y), (origT1.z - origT3.z));

		const glm::vec2& rambda_ = p.paras.at(i);
		const glm::vec3 rambda(rambda_.x, rambda_.y, 1 - rambda_.x - rambda_.y);
		glm::vec3 new_r;
		new_r.x = rambda.x * origT0.x + rambda.y * origT1.x + rambda.z * origT2.x;
		new_r.y = rambda.x * origT0.y + rambda.y * origT1.y + rambda.z * origT2.y;
		new_r.z = rambda.x * origT0.z + rambda.y * origT1.z + rambda.z * origT2.z;
		Vertex baryCenter(new_r.x, new_r.y, new_r.z);
		Vertex& v = quadV.at(m_v.at(i));
		v.x = baryCenter.x; v.y = baryCenter.y; v.z = baryCenter.z;
	}
}

void Patch::GetInnerVertices(const std::vector<unsigned long>& boundaryVertices, std::vector<unsigned long>& innerVertices)
{
	for (size_t i = 0; i < m_v.size(); i++)
	{
		const unsigned long v = m_v.at(i);
		if (std::find(boundaryVertices.begin(), boundaryVertices.end(), v) == boundaryVertices.end())
		{
			innerVertices.push_back(v);
		}
	}
	innerVertices.resize(innerVertices.size());
}
void Patch::GetBoundaryVertices(std::vector<unsigned long>& boundaryVertices)
{

	for (size_t i = 0; i < m_e.size(); i++)
	{
		const Edge& e = m_e.at(i);
		boundaryVertices.push_back(e.p0);
		boundaryVertices.push_back(e.p1);
	}
	std::sort(boundaryVertices.begin(), boundaryVertices.end());
	std::vector<unsigned long>::iterator iter = std::unique(boundaryVertices.begin(), boundaryVertices.end());
	boundaryVertices.resize(std::distance(boundaryVertices.begin(), iter));
}
void Patch::Smooth()
{
	std::vector<unsigned long> innerVertices;
	std::vector<unsigned long> boundaryVertices;
	GetBoundaryVertices(boundaryVertices);
	GetInnerVertices(boundaryVertices, innerVertices);
	//std::cout << "bound: " << boundaryVertices.size() << " inner: " << innerVertices.size() << std::endl;
	for (size_t i = 0; i < innerVertices.size(); i++)
	{
		const unsigned long vi = innerVertices.at(i);
		Vertex& v = m_pMesh->V.at(vi);
		const std::pair<std::multimap<unsigned long, unsigned long>::iterator,
		std::multimap<unsigned long, unsigned long>::iterator> ret = m_pMesh->VI_CI.equal_range(vi);

		glm::vec3 sum(0.0f, 0.0f, 0.0f);
		int count = 0;
		for (std::multimap<unsigned long, unsigned long>::iterator iter = ret.first; iter != ret.second; ++iter)
		{
			const unsigned long ci = iter->second;
			const Cell& cell = m_pMesh->C.at(ci);
			glm::vec3 center;
			m_pMesh->GetCenterPointOfCell(cell, center);
			sum = sum + center;
			count++;
		}
		v.x = sum.x/count;
		v.y = sum.y/count;
		v.z = sum.z/count;
	}
}

