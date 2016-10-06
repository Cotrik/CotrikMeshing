/*
 * PolycubeMesh.cpp
 *
 *  Created on: Jan 3, 2015
 *      Author: cotrik
 */

#include <iomanip>
#include <iostream>
#include <fstream>
#include <boost/smart_ptr.hpp>
#include "PolycubeMesh.h"
#include "MeshFileWriter.h"
#include "glm/glm.hpp"
#include "glm/gtx/intersect.hpp"
#include "include.h"

PolycubeMesh::PolycubeMesh(const std::vector<Vertex>& v, const std::vector<Cell>& c,  const ElementType cellType)
: Mesh(v, c, cellType)
{

}

PolycubeMesh::PolycubeMesh(const Mesh& mesh)
: Mesh(mesh)
{
}

PolycubeMesh::PolycubeMesh()
{
	// TODO Auto-generated constructor stub

}

PolycubeMesh::~PolycubeMesh()
{
	// TODO Auto-generated destructor stub
}

void PolycubeMesh::GetPlane_Face_MultiMap()
{
	for (std::vector<Face>::iterator iterSurface = surface.begin(); iterSurface != surface.end(); ++iterSurface)
	{
		const Vertex& v0 = V[iterSurface->at(0)];
		const Vertex& v1 = V[iterSurface->at(1)];
		const Vertex& v2 = V[iterSurface->at(2)];
		const Vector p0p1(v1.x - v0.x, v1.y - v0.y, v1.z - v0.z);
		const Vector p0p2(v2.x - v0.x, v2.y - v0.y, v2.z - v0.z);
		const Vector normal(p0p1.j * p0p2.k - p0p1.k * p0p2.j,
				p0p1.k * p0p2.i - p0p1.i * p0p2.k,
				p0p1.i * p0p2.j - p0p1.j * p0p2.i);

		const Point point(v0.y, v0.y, v0.z);
		const Plane plane(point, normal);

		m_plane_Face_Multimap.insert(std::pair<Plane, Face*>(plane, &*iterSurface));
	}
}

void PolycubeMesh::GetXYZPlane_Face_MultiMap()
{
	const pf::iterator iterBegin = m_plane_Face_Multimap.begin();
	const pf::iterator iterEnd = m_plane_Face_Multimap.end();
	std::multimap<Plane, Face*>::iterator iter = iterBegin;
	// for each plane
	while (iter != iterEnd)
	{
		glm::vec3 center(0.0f, 0.0f, 0.0f);
		GetCenterPointOfCell(*iter->second, center);
		if (fabs(iter->first.a) > 0.9)
		{
			const Plane p(1.0f, 0.0f, 0.0f, center.x);
			std::pair<Plane, Face*> item(p, iter->second);
			m_xPlane_Face_Multimap.insert(item);
			const Face& face = *iter->second;
			for (unsigned long i = 0; i < face.size(); i++)
			{
				V[face.at(i)].vinfo.dir = X_AXIS;
			}
		}
		if (fabs(iter->first.b) > 0.9)
		{
			const Plane p(0.0f, 1.0f, 0.0f, center.y);
			std::pair<Plane, Face*> item(p, iter->second);
			m_yPlane_Face_Multimap.insert(item);
			const Face& face = *iter->second;
			for (unsigned long i = 0; i < face.size(); i++)
			{
				V[face.at(i)].vinfo.dir = Y_AXIS;
			}
		}
		if (fabs(iter->first.c) > 0.9)
		{
			const Plane p(0.0f, 0.0f, 1.0f, center.z);
			std::pair<Plane, Face*> item(p, iter->second);
			m_zPlane_Face_Multimap.insert(item);
			const Face& face = *iter->second;
			for (unsigned long i = 0; i < face.size(); i++)
			{
				V[face.at(i)].vinfo.dir = Z_AXIS;
			}
		}

		++iter;
	}
}

void PolycubeMesh::GetSinglePlaneFaceMultimap(const pf& plane_Face_Multimap, df& singlePlane_Face_Multimap, const float delta/* = 0.005*/)
{
	const pf::const_iterator iterBegin = plane_Face_Multimap.begin();
	const pf::const_iterator iterEnd = plane_Face_Multimap.end();
	pf::const_iterator iter = iterBegin;
	// for each plane
	while (iter != iterEnd)
	{
		std::pair<df::iterator, df::iterator> ret;

		const df::const_iterator iterSinglePlaneBegin = singlePlane_Face_Multimap.begin();
		const df::const_iterator iterSinglePlaneEnd = singlePlane_Face_Multimap.end();
		df::const_iterator iterSinglePlane = iterSinglePlaneBegin;

		while (iterSinglePlane != iterSinglePlaneEnd)
		{
			if (fabs(iter->first.d - iterSinglePlane->first) < 2*1e-5)
			{
				singlePlane_Face_Multimap.insert(std::pair<double, Face*>(iterSinglePlane->first, iter->second));
				break;
			}
			++iterSinglePlane;
		}
		if (iterSinglePlane == iterSinglePlaneEnd)
		{
			singlePlane_Face_Multimap.insert(std::pair<double, Face*>(iter->first.d, iter->second));
		}
		++iter;
	}
}

void PolycubeMesh::GetPlane_Face_Multimap(df& m_plane_Face_Multimap, const DIRECTION dir)
{
	df theSameXPlane_Face_Multimap_Temp;
	for (std::vector<Face>::iterator iterSurface = surface.begin(); iterSurface != surface.end(); ++iterSurface)
	{
		Face& face = *iterSurface;
		double mean = 0;
		for (Face::iterator iterFace = face.begin(); iterFace != face.end(); ++iterFace)
		{
			if (dir == X_AXIS)
				mean += V[*iterFace].x;
			else if (dir == Y_AXIS)
				mean += V[*iterFace].y;
			else if (dir == Z_AXIS)
				mean += V[*iterFace].z;
		}
		theSameXPlane_Face_Multimap_Temp.insert(std::pair<double, Face*>(mean / face.size(), &face));
	}

	const df::const_iterator iterBegin = theSameXPlane_Face_Multimap_Temp.begin();
	const df::const_iterator iterEnd = theSameXPlane_Face_Multimap_Temp.end();
	df::const_iterator iterLastime = iterBegin;
	for (df::const_iterator iter = iterBegin; iter != iterEnd; ++iter)
	{
		if (fabs(iter->first - iterLastime->first) < 1e-5)
		{
			m_plane_Face_Multimap.insert(std::pair<double, Face*>(iterLastime->first, iter->second));
		}
		else
		{
			m_plane_Face_Multimap.insert(std::pair<double, Face*>(iter->first, iter->second));
			iterLastime = iter;
		}
	}
}

void PolycubeMesh::GetBoundaryLinesOfMeshPlane(const std::vector<Face*>& vecSinglePlaneFace, std::vector<Line>& vecBoundaryLine)
{
	const std::vector<Face*>::const_iterator iterFaceBegin = vecSinglePlaneFace.begin();
	const std::vector<Face*>::const_iterator iterFaceEnd = vecSinglePlaneFace.end();
	for (std::vector<Face*>::const_iterator iterFace = iterFaceBegin;	iterFace != iterFaceEnd; ++iterFace)
	{
		const Face& face = **iterFace;
		const size_t faceSize = face.size();
		for (unsigned int i = 0; i < faceSize; i++)
		{
			const Line line(face.at(i), face.at((i+1)%faceSize));
			std::vector<Line>::iterator iterLine = std::find(vecBoundaryLine.begin(), vecBoundaryLine.end(), line);
			if (iterLine != vecBoundaryLine.end())
			{
				iterLine = vecBoundaryLine.erase(iterLine);
			}
			else
			{
				vecBoundaryLine.push_back(line);
			}
		}
	}
}

void PolycubeMesh::GetPlaneBoundaryLine(const df& singlePlane_Face_Multimap, std::vector<Line>& vecBoundaryLine)
{
	const df::const_iterator iterSinglePlaneBegin = singlePlane_Face_Multimap.begin();
	const df::const_iterator iterSinglePlaneEnd = singlePlane_Face_Multimap.end();
	df::const_iterator iterSinglePlane = iterSinglePlaneBegin;
	while (iterSinglePlane != iterSinglePlaneEnd)
	{
		std::vector<Face*> vecSinglePlaneFace;

		std::pair<df::const_iterator, df::const_iterator> ret =
				singlePlane_Face_Multimap.equal_range(iterSinglePlane->first);
		for (iterSinglePlane = ret.first; iterSinglePlane != ret.second; ++iterSinglePlane)
		{
			vecSinglePlaneFace.push_back(iterSinglePlane->second);
		}

		GetBoundaryLinesOfMeshPlane(vecSinglePlaneFace, vecBoundaryLine);
	}
}

void PolycubeMesh::GetPlaneMultimap(const std::vector<Line>& vecBoundaryLine, vl& planeLineMultimap)
{
	const std::vector<Line>::const_iterator iterBegin = vecBoundaryLine.begin();
	const std::vector<Line>::const_iterator iterEnd = vecBoundaryLine.end();
	for (std::vector<Line>::const_iterator iter = iterBegin; iter != iterEnd; ++iter)
	{
		const Vertex& p0 = V[iter->p0];
		const Vertex& p1 = V[iter->p1];
		Vector v(p1.x - p0.x, p1.y - p0.y, p1.z - p0.z);
		if (fabs(v.i) < 1e-4) v.i = 0;
		if (fabs(v.j) < 1e-4) v.j = 0;
		if (fabs(v.k) < 1e-4) v.k = 0;
		unsigned int index = std::distance(iterBegin, iter);
		const Line& line = vecBoundaryLine.at(index);
		planeLineMultimap.insert(std::pair<Vector, Line*>(v, (Line*)&line));
	}
}

void PolycubeMesh::Getxy_xzPlaneLineMultimap(const vl& xPlaneLineMultimap, dl& xyPlaneLineMultimap, dl& xzPlaneLineMultimap)
{
	const vl::const_iterator iterBegin = xPlaneLineMultimap.begin();
	const vl::const_iterator iterEnd = xPlaneLineMultimap.end();
	for (std::multimap<Vector, Line*>::const_iterator iter = iterBegin; iter != iterEnd; ++iter)
	{
		if (iter->first.j == 0)
		{
			xyPlaneLineMultimap.insert(std::pair<double, Line*>(iter->first.k, iter->second));
			const Line& line = *iter->second;
			V[line.p0].vinfo.dir = Z_LINE;
			V[line.p1].vinfo.dir = Z_LINE;
		}
		if (iter->first.k == 0)
		{
			xzPlaneLineMultimap.insert(std::pair<double, Line*>(iter->first.j, iter->second));
			const Line& line = *iter->second;
			V[line.p0].vinfo.dir = Y_LINE;
			V[line.p1].vinfo.dir = Y_LINE;
		}
	}
}

void PolycubeMesh::Getyx_yzPlaneLineMultimap(const vl& yPlaneLineMultimap, dl& yxPlaneLineMultimap, dl& yzPlaneLineMultimap)
{
	const vl::const_iterator iterBegin = yPlaneLineMultimap.begin();
	const vl::const_iterator iterEnd = yPlaneLineMultimap.end();
	for (std::multimap<Vector, Line*>::const_iterator iter = iterBegin; iter != iterEnd; ++iter)
	{
		if (iter->first.i == 0)
		{
			yxPlaneLineMultimap.insert(std::pair<double, Line*>(iter->first.k, iter->second));
			const Line& line = *iter->second;
			V[line.p0].vinfo.dir = Z_LINE;
			V[line.p1].vinfo.dir = Z_LINE;
		}
		if (iter->first.k == 0)
		{
			yzPlaneLineMultimap.insert(std::pair<double, Line*>(iter->first.i, iter->second));
			const Line& line = *iter->second;
			V[line.p0].vinfo.dir = X_LINE;
			V[line.p1].vinfo.dir = X_LINE;
		}
	}
}

void PolycubeMesh::Getzx_zyPlaneLineMultimap(const vl& zPlaneLineMultimap, dl& zxPlaneLineMultimap, dl& zyPlaneLineMultimap)
{
	const vl::const_iterator iterBegin = zPlaneLineMultimap.begin();
	const vl::const_iterator iterEnd = zPlaneLineMultimap.end();
	for (std::multimap<Vector, Line*>::const_iterator iter = iterBegin; iter != iterEnd; ++iter)
	{
		if (iter->first.i == 0)
		{
			zxPlaneLineMultimap.insert(std::pair<double, Line*>(iter->first.j, iter->second));
		}
		if (iter->first.j == 0)
		{
			zyPlaneLineMultimap.insert(std::pair<double, Line*>(iter->first.i, iter->second));
		}
	}
}

void PolycubeMesh::GetXYPlaneCorner(const dl& xyPlaneLineMultimap, std::vector<unsigned long>& vec_xyPlaneCorner)
{
	const dl::const_iterator iterBegin = xyPlaneLineMultimap.begin();
	const dl::const_iterator iterEnd = xyPlaneLineMultimap.end();
	for (dl::const_iterator iter = iterBegin; iter != iterEnd; ++iter)
	{
		std::vector<unsigned long>::iterator iterPoint = std::find(vec_xyPlaneCorner.begin(), vec_xyPlaneCorner.end(), (*iter->second).p0);
		if (iterPoint != vec_xyPlaneCorner.end())
		{
			iterPoint = vec_xyPlaneCorner.erase(iterPoint);
		}
		else
		{
			vec_xyPlaneCorner.push_back((*iter->second).p0);
		}

		iterPoint = std::find(vec_xyPlaneCorner.begin(), vec_xyPlaneCorner.end(), (*iter->second).p1);
		if (iterPoint != vec_xyPlaneCorner.end())
		{
			iterPoint = vec_xyPlaneCorner.erase(iterPoint);
		}
		else
		{
			vec_xyPlaneCorner.push_back((*iter->second).p1);
		}
	}
}

void PolycubeMesh::GetCorners()
{
	if (m_plane_Face_Multimap.empty())
	{
		//AjustSurfaceToIntegerPlane(0.005);
		GetPlane_Face_MultiMap();
		GetXYZPlane_Face_MultiMap();

		GetSinglePlaneFaceMultimap(m_xPlane_Face_Multimap, m_xSinglePlane_Face_Multimap);
		GetSinglePlaneFaceMultimap(m_yPlane_Face_Multimap, m_ySinglePlane_Face_Multimap);
		GetSinglePlaneFaceMultimap(m_zPlane_Face_Multimap, m_zSinglePlane_Face_Multimap);

		GetPlaneBoundaryLine(m_xSinglePlane_Face_Multimap, m_vecXBoundaryLine);
		GetPlaneBoundaryLine(m_ySinglePlane_Face_Multimap, m_vecYBoundaryLine);
		GetPlaneBoundaryLine(m_zSinglePlane_Face_Multimap, m_vecZBoundaryLine);

		GetPlaneMultimap(m_vecXBoundaryLine, m_xPlaneLineMultimap);
		GetPlaneMultimap(m_vecYBoundaryLine, m_yPlaneLineMultimap);
		GetPlaneMultimap(m_vecZBoundaryLine, m_zPlaneLineMultimap);

		std::multimap<double, Line*> yxPlaneLineMultimap;
		Getxy_xzPlaneLineMultimap(m_xPlaneLineMultimap, m_xyPlaneLineMultimap, m_xzPlaneLineMultimap);
		Getyx_yzPlaneLineMultimap(m_yPlaneLineMultimap,   yxPlaneLineMultimap, m_yzPlaneLineMultimap);

		GetXYPlaneCorner(m_xyPlaneLineMultimap, m_vec_xyPlaneCorner);

		std::cout << "m_vec_xyPlaneCorner: " << m_vec_xyPlaneCorner.size() << std::endl;
	}
}

//std::ofstream nodeInfo("node_info.txt");
void AdjustCoordinateValue(const char* dir, const float& delta, const glm::vec3& normal, const Vertex& minVertex,
		Vertex* pVertex0, Vertex* pVertex1, Vertex* pVertex2)
{
	float position = fabs(pVertex0->x - minVertex.x)/delta;
	float unfittedPortion = position - int(position);
	if (unfittedPortion < 0.1f)
		return;
	float compensation = (1.0f - unfittedPortion + 0.10f)*delta;

	if (normal.x > 0)
	{
		pVertex0->vinfo.compensation.x = compensation;
		pVertex1->vinfo.compensation.x = compensation;
		pVertex2->vinfo.compensation.x = compensation;
	}
	else
	{
		pVertex0->vinfo.compensation.x = -compensation;
		pVertex1->vinfo.compensation.x = -compensation;
		pVertex2->vinfo.compensation.x = -compensation;
	}
}

const float BOUNDARY_ERROR = 0.1f;
float DELTA = 0;
void PolycubeMesh::AdjustPolycubeSurfaceToHoldAllSmallCubesVolume(const float& delta, const Vertex& minVertex, const char* outputFilename)
{
	DELTA = delta;
	const std::vector<Face>::const_iterator iterSurfaceBegin = surface.begin();
	const std::vector<Face>::const_iterator iterSurfaceEnd = surface.end();

	const std::vector<Vertex>::const_iterator iterVertexBegin = V.begin();
	const std::vector<Vertex>::const_iterator iterVertexEnd = V.end();

	int i = 0;
	for (std::vector<Face>::const_iterator iterSurface = iterSurfaceBegin; iterSurface != iterSurfaceEnd; ++iterSurface)
	{
		Vertex* pVertex0 = &V[iterSurface->at(0)];
		Vertex* pVertex1 = &V[iterSurface->at(1)];
		Vertex* pVertex2 = &V[iterSurface->at(2)];

		const glm::vec3 v0(pVertex0->x, pVertex0->y, pVertex0->z);
		const glm::vec3 v1(pVertex1->x, pVertex1->y, pVertex1->z);
		const glm::vec3 v2(pVertex2->x, pVertex2->y, pVertex2->z);

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

		if (dir == X_AXIS)
		{
			float position0 = fabs(pVertex0->x - minVertex.x)/delta;
			float position1 = fabs(pVertex1->x - minVertex.x)/delta;
			float position2 = fabs(pVertex2->x - minVertex.x)/delta;
			float unfittedPortion0 = position0 - int(position0);
			float unfittedPortion1 = position1 - int(position1);
			float unfittedPortion2 = position2 - int(position2);
			if (normal.x > 0)
			{
				float compensation0 = (1.0f - unfittedPortion0 + BOUNDARY_ERROR)*delta;
				float compensation1 = (1.0f - unfittedPortion1 + BOUNDARY_ERROR)*delta;
				float compensation2 = (1.0f - unfittedPortion2 + BOUNDARY_ERROR)*delta;
				pVertex0->vinfo.compensation.x = compensation0 > delta ? compensation0 - delta : compensation0;
				pVertex1->vinfo.compensation.x = compensation1 > delta ? compensation1 - delta : compensation1;
				pVertex2->vinfo.compensation.x = compensation2 > delta ? compensation2 - delta : compensation2;
			}
			else
			{
				float compensation0 = (-unfittedPortion0 - BOUNDARY_ERROR)*delta;
				float compensation1 = (-unfittedPortion1 - BOUNDARY_ERROR)*delta;
				float compensation2 = (-unfittedPortion2 - BOUNDARY_ERROR)*delta;
				pVertex0->vinfo.compensation.x = compensation0 < -delta ? compensation0 + delta : compensation0;
				pVertex1->vinfo.compensation.x = compensation1 < -delta ? compensation1 + delta : compensation1;
				pVertex2->vinfo.compensation.x = compensation2 < -delta ? compensation2 + delta : compensation2;
			}
		}
		else if (dir == Y_AXIS)
		{
			float position0 = fabs(pVertex0->y - minVertex.y)/delta;
			float position1 = fabs(pVertex1->y - minVertex.y)/delta;
			float position2 = fabs(pVertex2->y - minVertex.y)/delta;
			float unfittedPortion0 = position0 - int(position0);
			float unfittedPortion1 = position1 - int(position1);
			float unfittedPortion2 = position2 - int(position2);

			if (normal.y > 0)
			{
				float compensation0 = (1.0f - unfittedPortion0 + BOUNDARY_ERROR)*delta;
				float compensation1 = (1.0f - unfittedPortion1 + BOUNDARY_ERROR)*delta;
				float compensation2 = (1.0f - unfittedPortion2 + BOUNDARY_ERROR)*delta;
				pVertex0->vinfo.compensation.y = compensation0 > delta ? compensation0 - delta : compensation0;
				pVertex1->vinfo.compensation.y = compensation1 > delta ? compensation1 - delta : compensation1;
				pVertex2->vinfo.compensation.y = compensation2 > delta ? compensation2 - delta : compensation2;
			}
			else
			{
				float compensation0 = (-unfittedPortion0 - BOUNDARY_ERROR)*delta;
				float compensation1 = (-unfittedPortion1 - BOUNDARY_ERROR)*delta;
				float compensation2 = (-unfittedPortion2 - BOUNDARY_ERROR)*delta;
				pVertex0->vinfo.compensation.y = compensation0 < -delta ? compensation0 + delta : compensation0;
				pVertex1->vinfo.compensation.y = compensation1 < -delta ? compensation1 + delta : compensation1;
				pVertex2->vinfo.compensation.y = compensation2 < -delta ? compensation2 + delta : compensation2;
			}
		}
		else if (dir == Z_AXIS)
		{
			float position0 = fabs(pVertex0->z - minVertex.z)/delta;
			float position1 = fabs(pVertex1->z - minVertex.z)/delta;
			float position2 = fabs(pVertex2->z - minVertex.z)/delta;
			float unfittedPortion0 = position0 - int(position0);
			float unfittedPortion1 = position1 - int(position1);
			float unfittedPortion2 = position2 - int(position2);
			if (normal.z > 0)
			{
				float compensation0 = (1.0f - unfittedPortion0 + BOUNDARY_ERROR)*delta;
				float compensation1 = (1.0f - unfittedPortion1 + BOUNDARY_ERROR)*delta;
				float compensation2 = (1.0f - unfittedPortion2 + BOUNDARY_ERROR)*delta;
				pVertex0->vinfo.compensation.z = compensation0 > delta ? compensation0 - delta : compensation0;
				pVertex1->vinfo.compensation.z = compensation1 > delta ? compensation1 - delta : compensation1;
				pVertex2->vinfo.compensation.z = compensation2 > delta ? compensation2 - delta : compensation2;
			}
			else
			{
				float compensation0 = (-unfittedPortion0 - BOUNDARY_ERROR)*delta;
				float compensation1 = (-unfittedPortion1 - BOUNDARY_ERROR)*delta;
				float compensation2 = (-unfittedPortion2 - BOUNDARY_ERROR)*delta;
				pVertex0->vinfo.compensation.z = compensation0 < -delta ? compensation0 + delta : compensation0;
				pVertex1->vinfo.compensation.z = compensation1 < -delta ? compensation1 + delta : compensation1;
				pVertex2->vinfo.compensation.z = compensation2 < -delta ? compensation2 + delta : compensation2;
			}
		}
//		nodeInfo << "Triangle" << i++ << " " << std::fixed << std::setprecision(5) << " normal(" << normal.x << " "
//				<< normal.y << " " << normal.z << ") " << std::endl;
//		nodeInfo << std::fixed << std::setprecision(5) << "node_" << iterSurface->at(0) << "("
//				<< pVertex0->x - minVertex.x << " " << pVertex0->y - minVertex.y << " " << pVertex0->z - minVertex.z
//				<< ") compensation(" << pVertex0->vinfo.compensation.x << " " << pVertex0->vinfo.compensation.y << " "
//				<< pVertex0->vinfo.compensation.z << ")" << std::endl;
//		nodeInfo << std::fixed << std::setprecision(5) << "node_" << iterSurface->at(1) << "("
//				<< pVertex1->x - minVertex.x << " " << pVertex1->y - minVertex.y << " " << pVertex1->z - minVertex.z
//				<< ") compensation(" << pVertex1->vinfo.compensation.x << " " << pVertex1->vinfo.compensation.y << " "
//				<< pVertex1->vinfo.compensation.z << ")" << std::endl;
//		nodeInfo << std::fixed << std::setprecision(5) << "node_" << iterSurface->at(2) << "("
//				<< pVertex2->x - minVertex.x << " " << pVertex2->y - minVertex.y << " " << pVertex2->z - minVertex.z
//				<< ") compensation(" << pVertex2->vinfo.compensation.x << " " << pVertex2->vinfo.compensation.y << " "
//				<< pVertex2->vinfo.compensation.z << ")" << std::endl;
	}
//	for (std::vector<Face>::const_iterator iterSurface = iterSurfaceBegin; iterSurface != iterSurfaceEnd; ++iterSurface)
//	{
//		Vertex* pVertex0 = &V[iterSurface->at(0)];
//		Vertex* pVertex1 = &V[iterSurface->at(1)];
//		Vertex* pVertex2 = &V[iterSurface->at(2)];
//
//		if (!pVertex0->vinfo.bAjusted)
//		{
//			pVertex0->x += pVertex0->vinfo.compensation.x;
//			pVertex0->y += pVertex0->vinfo.compensation.y;
//			pVertex0->z += pVertex0->vinfo.compensation.z;
//			pVertex0->vinfo.bAjusted = true;
//			nodeInfo << std::fixed << std::setprecision(5) << "--node_" << iterSurface->at(0) << "("
//					<< pVertex0->x - minVertex.x << " " << pVertex0->y - minVertex.y << " " << pVertex0->z - minVertex.z
//					<< ") compensation(" << pVertex0->vinfo.compensation.x << " " << pVertex0->vinfo.compensation.y << " "
//					<< pVertex0->vinfo.compensation.z << ")" << std::endl;
//		}
//		if (!pVertex1->vinfo.bAjusted)
//		{
//			pVertex1->x += pVertex1->vinfo.compensation.x;
//			pVertex1->y += pVertex1->vinfo.compensation.y;
//			pVertex1->z += pVertex1->vinfo.compensation.z;
//			pVertex1->vinfo.bAjusted = true;
//			nodeInfo << std::fixed << std::setprecision(5) << "--node_" << iterSurface->at(1) << "("
//					<< pVertex1->x - minVertex.x << " " << pVertex1->y - minVertex.y << " " << pVertex1->z - minVertex.z
//					<< ") compensation(" << pVertex1->vinfo.compensation.x << " " << pVertex1->vinfo.compensation.y << " "
//					<< pVertex1->vinfo.compensation.z << ")" << std::endl;
//		}
//		if (!pVertex2->vinfo.bAjusted)
//		{
//			pVertex2->x += pVertex2->vinfo.compensation.x;
//			pVertex2->y += pVertex2->vinfo.compensation.y;
//			pVertex2->z += pVertex2->vinfo.compensation.z;
//			pVertex2->vinfo.bAjusted = true;
//			nodeInfo << std::fixed << std::setprecision(5) << "--node_" << iterSurface->at(2) << "("
//					<< pVertex2->x - minVertex.x << " " << pVertex2->y - minVertex.y << " " << pVertex2->z - minVertex.z
//					<< ") compensation(" << pVertex2->vinfo.compensation.x << " " << pVertex2->vinfo.compensation.y << " "
//					<< pVertex2->vinfo.compensation.z << ")" << std::endl;
//		}
//	}
}

void PolycubeMesh::AdjustPolycubeSurfaceToHoldAllSmallCubesVolume1(const float& delta, const Vertex& minVertex, const char* outputFilename)
{
	DELTA = delta;
	const std::vector<Face>::const_iterator iterSurfaceBegin = surface.begin();
	const std::vector<Face>::const_iterator iterSurfaceEnd = surface.end();

	const std::vector<Vertex>::const_iterator iterVertexBegin = V.begin();
	const std::vector<Vertex>::const_iterator iterVertexEnd = V.end();

	int i = 0;
	for (std::vector<Face>::const_iterator iterSurface = iterSurfaceBegin; iterSurface != iterSurfaceEnd; ++iterSurface)
	{
		Vertex* pVertex0 = &V[iterSurface->at(0)];
		Vertex* pVertex1 = &V[iterSurface->at(1)];
		Vertex* pVertex2 = &V[iterSurface->at(2)];

		const glm::vec3 v0(pVertex0->x, pVertex0->y, pVertex0->z);
		const glm::vec3 v1(pVertex1->x, pVertex1->y, pVertex1->z);
		const glm::vec3 v2(pVertex2->x, pVertex2->y, pVertex2->z);

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

		if (dir == X_AXIS)
		{
			if (normal.x > 0)
			{
				float compensation = BOUNDARY_ERROR * delta;
				pVertex0->vinfo.compensation.x = compensation;
				pVertex1->vinfo.compensation.x = compensation;
				pVertex2->vinfo.compensation.x = compensation;
			}
			else
			{
				float compensation = 0.0 - BOUNDARY_ERROR * delta;
				pVertex0->vinfo.compensation.x = compensation;
				pVertex1->vinfo.compensation.x = compensation;
				pVertex2->vinfo.compensation.x = compensation;
			}
		}
		else if (dir == Y_AXIS)
		{
			if (normal.y > 0)
			{
				float compensation = BOUNDARY_ERROR * delta;
				pVertex0->vinfo.compensation.y = compensation;
				pVertex1->vinfo.compensation.y = compensation;
				pVertex2->vinfo.compensation.y = compensation;
			}
			else
			{
				float compensation = 0.0 - BOUNDARY_ERROR * delta;
				pVertex0->vinfo.compensation.y = compensation;
				pVertex1->vinfo.compensation.y = compensation;
				pVertex2->vinfo.compensation.y = compensation;
			}
		}
		else if (dir == Z_AXIS)
		{
			if (normal.z > 0)
			{
				float compensation = BOUNDARY_ERROR * delta;
				pVertex0->vinfo.compensation.z = compensation;
				pVertex1->vinfo.compensation.z = compensation;
				pVertex2->vinfo.compensation.z = compensation;
			}
			else
			{
				float compensation = 0.0 - BOUNDARY_ERROR * delta;
				pVertex0->vinfo.compensation.z = compensation;
				pVertex1->vinfo.compensation.z = compensation;
				pVertex2->vinfo.compensation.z = compensation;
			}
		}
	}
}

void PolycubeMesh::AdjustPolycubeSurfaceToHoldAllSmallCubesVolume1_Reverse(const float& delta, const Vertex& minVertex, const char* outputFilename)
{
	DELTA = delta;
	const std::vector<Face>::const_iterator iterSurfaceBegin = surface.begin();
	const std::vector<Face>::const_iterator iterSurfaceEnd = surface.end();

	const std::vector<Vertex>::const_iterator iterVertexBegin = V.begin();
	const std::vector<Vertex>::const_iterator iterVertexEnd = V.end();

	int i = 0;
	for (std::vector<Face>::const_iterator iterSurface = iterSurfaceBegin; iterSurface != iterSurfaceEnd; ++iterSurface)
	{
		Vertex* pVertex0 = &V[iterSurface->at(0)];
		Vertex* pVertex1 = &V[iterSurface->at(1)];
		Vertex* pVertex2 = &V[iterSurface->at(2)];

		const glm::vec3 v0(pVertex0->x, pVertex0->y, pVertex0->z);
		const glm::vec3 v1(pVertex1->x, pVertex1->y, pVertex1->z);
		const glm::vec3 v2(pVertex2->x, pVertex2->y, pVertex2->z);

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

		if (dir == X_AXIS)
		{
			if (normal.x > 0)
			{
				float compensation = 0.0 - BOUNDARY_ERROR * delta;
				pVertex0->vinfo.compensation.x = compensation;
				pVertex1->vinfo.compensation.x = compensation;
				pVertex2->vinfo.compensation.x = compensation;
			}
			else
			{
				float compensation = BOUNDARY_ERROR * delta;
				pVertex0->vinfo.compensation.x = compensation;
				pVertex1->vinfo.compensation.x = compensation;
				pVertex2->vinfo.compensation.x = compensation;
			}
		}
		else if (dir == Y_AXIS)
		{
			if (normal.y > 0)
			{
				float compensation = 0.0 - BOUNDARY_ERROR * delta;
				pVertex0->vinfo.compensation.y = compensation;
				pVertex1->vinfo.compensation.y = compensation;
				pVertex2->vinfo.compensation.y = compensation;
			}
			else
			{
				float compensation = BOUNDARY_ERROR * delta;
				pVertex0->vinfo.compensation.y = compensation;
				pVertex1->vinfo.compensation.y = compensation;
				pVertex2->vinfo.compensation.y = compensation;
			}
		}
		else if (dir == Z_AXIS)
		{
			if (normal.z > 0)
			{
				float compensation = 0.0 - BOUNDARY_ERROR * delta;
				pVertex0->vinfo.compensation.z = compensation;
				pVertex1->vinfo.compensation.z = compensation;
				pVertex2->vinfo.compensation.z = compensation;
			}
			else
			{
				float compensation = BOUNDARY_ERROR * delta;
				pVertex0->vinfo.compensation.z = compensation;
				pVertex1->vinfo.compensation.z = compensation;
				pVertex2->vinfo.compensation.z = compensation;
			}
		}
	}
}

void PolycubeMesh::AdjustPolycubeSurfaceToHoldAllSmallCubesVolume1_ROUND(const float& delta, const Vertex& minVertex, const char* outputFilename)
{
	DELTA = delta;
	const std::vector<Face>::const_iterator iterSurfaceBegin = surface.begin();
	const std::vector<Face>::const_iterator iterSurfaceEnd = surface.end();

	const std::vector<Vertex>::const_iterator iterVertexBegin = V.begin();
	const std::vector<Vertex>::const_iterator iterVertexEnd = V.end();

	int i = 0;
	for (std::vector<Face>::const_iterator iterSurface = iterSurfaceBegin; iterSurface != iterSurfaceEnd; ++iterSurface)
	{
		Vertex* pVertex0 = &V[iterSurface->at(0)];
		Vertex* pVertex1 = &V[iterSurface->at(1)];
		Vertex* pVertex2 = &V[iterSurface->at(2)];

		const glm::vec3 v0(pVertex0->x, pVertex0->y, pVertex0->z);
		const glm::vec3 v1(pVertex1->x, pVertex1->y, pVertex1->z);
		const glm::vec3 v2(pVertex2->x, pVertex2->y, pVertex2->z);

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

		if (dir == X_AXIS)
		{
			if (normal.x > 0)
			{
				float compensation = delta;
				pVertex0->vinfo.compensation.x = compensation;
				pVertex1->vinfo.compensation.x = compensation;
				pVertex2->vinfo.compensation.x = compensation;
			}
			else
			{
				float compensation = -delta;
				pVertex0->vinfo.compensation.x = compensation;
				pVertex1->vinfo.compensation.x = compensation;
				pVertex2->vinfo.compensation.x = compensation;
			}
		}
		else if (dir == Y_AXIS)
		{
			if (normal.y > 0)
			{
				float compensation = delta;
				pVertex0->vinfo.compensation.y = compensation;
				pVertex1->vinfo.compensation.y = compensation;
				pVertex2->vinfo.compensation.y = compensation;
			}
			else
			{
				float compensation = -delta;
				pVertex0->vinfo.compensation.y = compensation;
				pVertex1->vinfo.compensation.y = compensation;
				pVertex2->vinfo.compensation.y = compensation;
			}
		}
		else if (dir == Z_AXIS)
		{
			if (normal.z > 0)
			{
				float compensation = delta;
				pVertex0->vinfo.compensation.z = compensation;
				pVertex1->vinfo.compensation.z = compensation;
				pVertex2->vinfo.compensation.z = compensation;
			}
			else
			{
				float compensation = -delta;
				pVertex0->vinfo.compensation.z = compensation;
				pVertex1->vinfo.compensation.z = compensation;
				pVertex2->vinfo.compensation.z = compensation;
			}
		}
	}
}

void PolycubeMesh::AdjustPolycubeSurfaceToHoldAllSmallCubesVolume2_ROUND(const float& delta, const Vertex& minVertex, const char* outputFilename)
{
    DELTA = delta;
    const std::vector<Face>::const_iterator iterSurfaceBegin = surface.begin();
    const std::vector<Face>::const_iterator iterSurfaceEnd = surface.end();

    const std::vector<Vertex>::const_iterator iterVertexBegin = V.begin();
    const std::vector<Vertex>::const_iterator iterVertexEnd = V.end();

    int i = 0;
    for (std::vector<Face>::const_iterator iterSurface = iterSurfaceBegin; iterSurface != iterSurfaceEnd; ++iterSurface)
    {
        Vertex* pVertex0 = &V[iterSurface->at(0)];
        Vertex* pVertex1 = &V[iterSurface->at(1)];
        Vertex* pVertex2 = &V[iterSurface->at(2)];

        const glm::vec3 v0(pVertex0->x, pVertex0->y, pVertex0->z);
        const glm::vec3 v1(pVertex1->x, pVertex1->y, pVertex1->z);
        const glm::vec3 v2(pVertex2->x, pVertex2->y, pVertex2->z);

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

        if (dir == X_AXIS)
        {
            if (normal.x > 0)
            {
                float compensation = 0;
                if (pVertex0->x < 0)
                    compensation = -delta;
                pVertex0->vinfo.compensation.x = compensation;
                pVertex1->vinfo.compensation.x = compensation;
                pVertex2->vinfo.compensation.x = compensation;
            }
            else
            {
                float compensation = delta;
                if (pVertex0->x < 0)
                    compensation = 0;
                pVertex0->vinfo.compensation.x = compensation;
                pVertex1->vinfo.compensation.x = compensation;
                pVertex2->vinfo.compensation.x = compensation;
            }
        }
        else if (dir == Y_AXIS)
        {
            if (normal.y > 0)
            {
                float compensation = 0;
                if (pVertex0->y < 0)
                    compensation = -delta;
                pVertex0->vinfo.compensation.y = compensation;
                pVertex1->vinfo.compensation.y = compensation;
                pVertex2->vinfo.compensation.y = compensation;
            }
            else
            {
                float compensation = delta;
                if (pVertex0->y < 0)
                    compensation = 0;
                pVertex0->vinfo.compensation.y = compensation;
                pVertex1->vinfo.compensation.y = compensation;
                pVertex2->vinfo.compensation.y = compensation;
            }
        }
        else if (dir == Z_AXIS)
        {
            if (normal.z > 0)
            {
                float compensation = 0;
                if (pVertex0->z < 0)
                    compensation = -delta;
                pVertex0->vinfo.compensation.z = compensation;
                pVertex1->vinfo.compensation.z = compensation;
                pVertex2->vinfo.compensation.z = compensation;
            }
            else
            {
                float compensation = delta;
                if (pVertex0->z < 0)
                    compensation = 0;
                pVertex0->vinfo.compensation.z = compensation;
                pVertex1->vinfo.compensation.z = compensation;
                pVertex2->vinfo.compensation.z = compensation;
            }
        }
    }
}

void PolycubeMesh::AdjustPolycubeSurfaceToHoldAllSmallCubesVolume_round(const float& delta, const Vertex& minVertex, const char* outputFilename)
{
	DELTA = delta;
	const std::vector<Face>::const_iterator iterSurfaceBegin = surface.begin();
	const std::vector<Face>::const_iterator iterSurfaceEnd = surface.end();

	const std::vector<Vertex>::const_iterator iterVertexBegin = V.begin();
	const std::vector<Vertex>::const_iterator iterVertexEnd = V.end();

	int i = 0;
	for (std::vector<Face>::const_iterator iterSurface = iterSurfaceBegin; iterSurface != iterSurfaceEnd; ++iterSurface)
	{
		Vertex* pVertex0 = &V[iterSurface->at(0)];
		Vertex* pVertex1 = &V[iterSurface->at(1)];
		Vertex* pVertex2 = &V[iterSurface->at(2)];

		const glm::vec3 v0(pVertex0->x, pVertex0->y, pVertex0->z);
		const glm::vec3 v1(pVertex1->x, pVertex1->y, pVertex1->z);
		const glm::vec3 v2(pVertex2->x, pVertex2->y, pVertex2->z);

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

		if (dir == X_AXIS)
		{
			if (normal.x > 0)
			{
				float compensation = BOUNDARY_ERROR * delta;
				pVertex0->vinfo.compensation.x = compensation;
				pVertex1->vinfo.compensation.x = compensation;
				pVertex2->vinfo.compensation.x = compensation;
			}
			else
			{
				float compensation = 0.0 - BOUNDARY_ERROR * delta;
				pVertex0->vinfo.compensation.x = compensation;
				pVertex1->vinfo.compensation.x = compensation;
				pVertex2->vinfo.compensation.x = compensation;
			}
		}
		else if (dir == Y_AXIS)
		{
			if (normal.y > 0)
			{
				float compensation = BOUNDARY_ERROR * delta;
				pVertex0->vinfo.compensation.y = compensation;
				pVertex1->vinfo.compensation.y = compensation;
				pVertex2->vinfo.compensation.y = compensation;
			}
			else
			{
				float compensation = 0.0 - BOUNDARY_ERROR * delta;
				pVertex0->vinfo.compensation.y = compensation;
				pVertex1->vinfo.compensation.y = compensation;
				pVertex2->vinfo.compensation.y = compensation;
			}
		}
		else if (dir == Z_AXIS)
		{
			if (normal.z > 0)
			{
				float compensation = BOUNDARY_ERROR * delta;
				pVertex0->vinfo.compensation.z = compensation;
				pVertex1->vinfo.compensation.z = compensation;
				pVertex2->vinfo.compensation.z = compensation;
			}
			else
			{
				float compensation = 0.0 - BOUNDARY_ERROR * delta;
				pVertex0->vinfo.compensation.z = compensation;
				pVertex1->vinfo.compensation.z = compensation;
				pVertex2->vinfo.compensation.z = compensation;
			}
		}
	}
}

void PolycubeMesh::AdjustPolycubeSurfaceToHoldAllSmallCubesVolume_quad(const float& delta, const Vertex& minVertex, const char* outputFilename)
{
	DELTA = delta;
	const std::vector<Face>::const_iterator iterSurfaceBegin = surface.begin();
	const std::vector<Face>::const_iterator iterSurfaceEnd = surface.end();

	const std::vector<Vertex>::const_iterator iterVertexBegin = V.begin();
	const std::vector<Vertex>::const_iterator iterVertexEnd = V.end();

	int i = 0;
	for (std::vector<Face>::const_iterator iterSurface = iterSurfaceBegin; iterSurface != iterSurfaceEnd; ++iterSurface)
	{
		Vertex* pVertex0 = &V[iterSurface->at(0)];
		Vertex* pVertex1 = &V[iterSurface->at(1)];
		Vertex* pVertex2 = &V[iterSurface->at(2)];
		Vertex* pVertex3 = &V[iterSurface->at(3)];

		const glm::vec3 v0(pVertex0->x, pVertex0->y, pVertex0->z);
		const glm::vec3 v1(pVertex1->x, pVertex1->y, pVertex1->z);
		const glm::vec3 v2(pVertex2->x, pVertex2->y, pVertex2->z);
		const glm::vec3 v3(pVertex3->x, pVertex3->y, pVertex3->z);

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

		if (dir == X_AXIS)
		{
			if (normal.x < 0)
			{
				float compensation = BOUNDARY_ERROR * delta;
				pVertex0->vinfo.compensation.x = compensation;
				pVertex1->vinfo.compensation.x = compensation;
				pVertex2->vinfo.compensation.x = compensation;
				pVertex3->vinfo.compensation.x = compensation;
			}
			else
			{
				float compensation = 0.0 - BOUNDARY_ERROR * delta;
				pVertex0->vinfo.compensation.x = compensation;
				pVertex1->vinfo.compensation.x = compensation;
				pVertex2->vinfo.compensation.x = compensation;
				pVertex3->vinfo.compensation.x = compensation;
			}
		}
		else if (dir == Y_AXIS)
		{
			if (normal.y < 0)
			{
				float compensation = BOUNDARY_ERROR * delta;
				pVertex0->vinfo.compensation.y = compensation;
				pVertex1->vinfo.compensation.y = compensation;
				pVertex2->vinfo.compensation.y = compensation;
				pVertex3->vinfo.compensation.y = compensation;
			}
			else
			{
				float compensation = 0.0 - BOUNDARY_ERROR * delta;
				pVertex0->vinfo.compensation.y = compensation;
				pVertex1->vinfo.compensation.y = compensation;
				pVertex2->vinfo.compensation.y = compensation;
				pVertex3->vinfo.compensation.y = compensation;
			}
		}
		else if (dir == Z_AXIS)
		{
			if (normal.z < 0)
			{
				float compensation = BOUNDARY_ERROR * delta;
				pVertex0->vinfo.compensation.z = compensation;
				pVertex1->vinfo.compensation.z = compensation;
				pVertex2->vinfo.compensation.z = compensation;
				pVertex3->vinfo.compensation.z = compensation;
			}
			else
			{
				float compensation = 0.0 - BOUNDARY_ERROR * delta;
				pVertex0->vinfo.compensation.z = compensation;
				pVertex1->vinfo.compensation.z = compensation;
				pVertex2->vinfo.compensation.z = compensation;
				pVertex3->vinfo.compensation.z = compensation;
			}
		}
	}
}

void PolycubeMesh::RoundSurface(const double eps/* = 1e-3*/)
{
    Magnify(eps);
    // Round Surface
    for (size_t id = 0; id < V.size(); id++) {
        Vertex& v = V.at(id);
        if (v.vinfo.bSurface) {
            v.x = int(v.x);
            v.y = int(v.y);
            v.z = int(v.z);
        }
    }
}

void PolycubeMesh::Magnify(const double eps/* = 1e-3*/)
{
    if (eps <= 0)
        std::cout << "eps Error!\n";

    ExtractSurface();
    int magnification = 1.0/eps;
    for (size_t id = 0; id < V.size(); id++) {
        Vertex& v = V.at(id);
        v.x *= magnification;
        v.y *= magnification;
        v.z *= magnification;
    }
}

void PolycubeMesh::RoundSurfaceAxis(const double eps)
{
    Magnify(eps);
    const std::vector<Face>::const_iterator iterSurfaceBegin = surface.begin();
    const std::vector<Face>::const_iterator iterSurfaceEnd = surface.end();

    const std::vector<Vertex>::const_iterator iterVertexBegin = V.begin();
    const std::vector<Vertex>::const_iterator iterVertexEnd = V.end();

    int i = 0;
    for (std::vector<Face>::const_iterator iterSurface = iterSurfaceBegin; iterSurface != iterSurfaceEnd; ++iterSurface)
    {
        Vertex* pVertex0 = &V[iterSurface->at(0)];
        Vertex* pVertex1 = &V[iterSurface->at(1)];
        Vertex* pVertex2 = &V[iterSurface->at(2)];

        const glm::vec3 v0(pVertex0->x, pVertex0->y, pVertex0->z);
        const glm::vec3 v1(pVertex1->x, pVertex1->y, pVertex1->z);
        const glm::vec3 v2(pVertex2->x, pVertex2->y, pVertex2->z);

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

        if (dir == X_AXIS)
        {
            pVertex0->x = int(pVertex0->x);
            pVertex1->x = int(pVertex1->x);
            pVertex2->x = int(pVertex2->x);
        }
        else if (dir == Y_AXIS)
        {
            pVertex0->y = int(pVertex0->y);
            pVertex1->y = int(pVertex1->y);
            pVertex2->y = int(pVertex2->y);
        }
        else if (dir == Z_AXIS)
        {
            pVertex0->z = int(pVertex0->z);
            pVertex1->z = int(pVertex1->z);
            pVertex2->z = int(pVertex2->z);
        }
    }
}
void PolycubeMesh::GenerateHexVertices(std::vector<Vertex>& hexV)
{
    //Hex vertices are generated first for the PolyCube corners, then edges, then facets and finally volume
    std::vector<int> cornerIds;
    std::vector<Edge> edges;
    std::vector<int> facetIds;
    GenerateHexVerticesForCorners(hexV, cornerIds);
    GenerateHexVerticesForEdges(hexV, edges);
    GenerateHexVerticesForFacets(hexV, facetIds);
}

void PolycubeMesh::GenerateHexVerticesForCorners(std::vector<Vertex>& hexV, std::vector<int>& cornerIds)
{
    for (size_t i = 0; i < m_vec_xyPlaneCorner.size(); i++) {
        unsigned long id = m_vec_xyPlaneCorner.at(i);
        cornerIds.push_back(id);
        const Vertex& v = V.at(id);
        hexV.push_back(v);
    }
}
void PolycubeMesh::GenerateHexVerticesForEdges(std::vector<Vertex>& hexV, std::vector<Edge>& edges)
{

}
void PolycubeMesh::GenerateHexVerticesForFacets(std::vector<Vertex>& hexV, std::vector<int>& facetIds)
{

}
void PolycubeMesh::AddCompensation()
{
	const std::vector<Face>::const_iterator iterSurfaceBegin = surface.begin();
	const std::vector<Face>::const_iterator iterSurfaceEnd = surface.end();
	for (std::vector<Face>::const_iterator iterSurface = iterSurfaceBegin; iterSurface != iterSurfaceEnd; ++iterSurface)
	{
		Vertex* pVertex0 = &V[iterSurface->at(0)];
		Vertex* pVertex1 = &V[iterSurface->at(1)];
		Vertex* pVertex2 = &V[iterSurface->at(2)];

		if (!pVertex0->vinfo.bAjusted)
		{
			pVertex0->x += pVertex0->vinfo.compensation.x;
			pVertex0->y += pVertex0->vinfo.compensation.y;
			pVertex0->z += pVertex0->vinfo.compensation.z;
			pVertex0->vinfo.bAjusted = true;
//			nodeInfo << std::fixed << std::setprecision(5) << "--node_" << iterSurface->at(0) << "("
//					<< (pVertex0->x - m_minVertex.x)/DELTA << " " << (pVertex0->y - m_minVertex.y)/DELTA << " " << (pVertex0->z - m_minVertex.z)/DELTA
//					<< ") compensation(" << pVertex0->vinfo.compensation.x << " " << pVertex0->vinfo.compensation.y << " "
//					<< pVertex0->vinfo.compensation.z << ")" << std::endl;
		}
		if (!pVertex1->vinfo.bAjusted)
		{
			pVertex1->x += pVertex1->vinfo.compensation.x;
			pVertex1->y += pVertex1->vinfo.compensation.y;
			pVertex1->z += pVertex1->vinfo.compensation.z;
			pVertex1->vinfo.bAjusted = true;
//			nodeInfo << std::fixed << std::setprecision(5) << "--node_" << iterSurface->at(1) << "("
//					<< (pVertex1->x - m_minVertex.x)/DELTA << " " << (pVertex1->y - m_minVertex.y)/DELTA << " " << (pVertex1->z - m_minVertex.z)/DELTA
//					<< ") compensation(" << pVertex1->vinfo.compensation.x << " " << pVertex1->vinfo.compensation.y << " "
//					<< pVertex1->vinfo.compensation.z << ")" << std::endl;
		}
		if (!pVertex2->vinfo.bAjusted)
		{
			pVertex2->x += pVertex2->vinfo.compensation.x;
			pVertex2->y += pVertex2->vinfo.compensation.y;
			pVertex2->z += pVertex2->vinfo.compensation.z;
			pVertex2->vinfo.bAjusted = true;
//			nodeInfo << std::fixed << std::setprecision(5) << "--node_" << iterSurface->at(2) << "("
//					<< (pVertex2->x - m_minVertex.x)/DELTA << " " << (pVertex2->y - m_minVertex.y)/DELTA << " " << (pVertex2->z - m_minVertex.z)/DELTA
//					<< ") compensation(" << pVertex2->vinfo.compensation.x << " " << pVertex2->vinfo.compensation.y << " "
//					<< pVertex2->vinfo.compensation.z << ")" << std::endl;
		}
	}
}

void PolycubeMesh::AddCompensation_round()
{
	const std::vector<Face>::const_iterator iterSurfaceBegin = surface.begin();
	const std::vector<Face>::const_iterator iterSurfaceEnd = surface.end();
	for (std::vector<Face>::const_iterator iterSurface = iterSurfaceBegin; iterSurface != iterSurfaceEnd; ++iterSurface)
	{
		Vertex* pVertex0 = &V[iterSurface->at(0)];
		Vertex* pVertex1 = &V[iterSurface->at(1)];
		Vertex* pVertex2 = &V[iterSurface->at(2)];

		if (!pVertex0->vinfo.bAjusted)
		{
			pVertex0->x += pVertex0->vinfo.compensation.x;
			pVertex0->y += pVertex0->vinfo.compensation.y;
			pVertex0->z += pVertex0->vinfo.compensation.z;
			pVertex0->vinfo.bAjusted = true;
#if _DEBUG
			nodeInfo << std::fixed << std::setprecision(5) << "--node_" << iterSurface->at(0) << "("
					<< (pVertex0->x - m_minVertex.x)/DELTA << " " << (pVertex0->y - m_minVertex.y)/DELTA << " " << (pVertex0->z - m_minVertex.z)/DELTA
					<< ") compensation(" << pVertex0->vinfo.compensation.x << " " << pVertex0->vinfo.compensation.y << " "
					<< pVertex0->vinfo.compensation.z << ")" << std::endl;
#endif
		}
		if (!pVertex1->vinfo.bAjusted)
		{
			pVertex1->x += pVertex1->vinfo.compensation.x;
			pVertex1->y += pVertex1->vinfo.compensation.y;
			pVertex1->z += pVertex1->vinfo.compensation.z;
			pVertex1->vinfo.bAjusted = true;
#if _DEBUG
			nodeInfo << std::fixed << std::setprecision(5) << "--node_" << iterSurface->at(1) << "("
					<< (pVertex1->x - m_minVertex.x)/DELTA << " " << (pVertex1->y - m_minVertex.y)/DELTA << " " << (pVertex1->z - m_minVertex.z)/DELTA
					<< ") compensation(" << pVertex1->vinfo.compensation.x << " " << pVertex1->vinfo.compensation.y << " "
					<< pVertex1->vinfo.compensation.z << ")" << std::endl;
#endif
		}
		if (!pVertex2->vinfo.bAjusted)
		{
			pVertex2->x += pVertex2->vinfo.compensation.x;
			pVertex2->y += pVertex2->vinfo.compensation.y;
			pVertex2->z += pVertex2->vinfo.compensation.z;
			pVertex2->vinfo.bAjusted = true;
#if _DEBUG
			nodeInfo << std::fixed << std::setprecision(5) << "--node_" << iterSurface->at(2) << "("
					<< (pVertex2->x - m_minVertex.x)/DELTA << " " << (pVertex2->y - m_minVertex.y)/DELTA << " " << (pVertex2->z - m_minVertex.z)/DELTA
					<< ") compensation(" << pVertex2->vinfo.compensation.x << " " << pVertex2->vinfo.compensation.y << " "
					<< pVertex2->vinfo.compensation.z << ")" << std::endl;
#endif
		}
	}
}

void PolycubeMesh::AddCompensation_quad()
{
	const std::vector<Face>::const_iterator iterSurfaceBegin = surface.begin();
	const std::vector<Face>::const_iterator iterSurfaceEnd = surface.end();
	for (std::vector<Face>::const_iterator iterSurface = iterSurfaceBegin; iterSurface != iterSurfaceEnd; ++iterSurface)
	{
		Vertex* pVertex0 = &V[iterSurface->at(0)];
		Vertex* pVertex1 = &V[iterSurface->at(1)];
		Vertex* pVertex2 = &V[iterSurface->at(2)];
		Vertex* pVertex3 = &V[iterSurface->at(3)];

		if (!pVertex0->vinfo.bAjusted)
		{
			pVertex0->x += pVertex0->vinfo.compensation.x;
			pVertex0->y += pVertex0->vinfo.compensation.y;
			pVertex0->z += pVertex0->vinfo.compensation.z;
			pVertex0->vinfo.bAjusted = true;
		}
		if (!pVertex1->vinfo.bAjusted)
		{
			pVertex1->x += pVertex1->vinfo.compensation.x;
			pVertex1->y += pVertex1->vinfo.compensation.y;
			pVertex1->z += pVertex1->vinfo.compensation.z;
			pVertex1->vinfo.bAjusted = true;
		}
		if (!pVertex2->vinfo.bAjusted)
		{
			pVertex2->x += pVertex2->vinfo.compensation.x;
			pVertex2->y += pVertex2->vinfo.compensation.y;
			pVertex2->z += pVertex2->vinfo.compensation.z;
			pVertex2->vinfo.bAjusted = true;
		}
		if (!pVertex3->vinfo.bAjusted)
		{
			pVertex3->x += pVertex3->vinfo.compensation.x;
			pVertex3->y += pVertex3->vinfo.compensation.y;
			pVertex3->z += pVertex3->vinfo.compensation.z;
			pVertex3->vinfo.bAjusted = true;
		}
	}
}

void PolycubeMesh::AddCompensation_quad_round()
{
	const std::vector<Face>::const_iterator iterSurfaceBegin = surface.begin();
	const std::vector<Face>::const_iterator iterSurfaceEnd = surface.end();
	for (std::vector<Face>::const_iterator iterSurface = iterSurfaceBegin; iterSurface != iterSurfaceEnd; ++iterSurface)
	{
		Vertex* pVertex0 = &V[iterSurface->at(0)];
		Vertex* pVertex1 = &V[iterSurface->at(1)];
		Vertex* pVertex2 = &V[iterSurface->at(2)];
		Vertex* pVertex3 = &V[iterSurface->at(3)];

		if (!pVertex0->vinfo.bAjusted)
		{
			pVertex0->x += pVertex0->vinfo.compensation.x > 0 ? 1.0f : -1.0f;
			pVertex0->y += pVertex0->vinfo.compensation.y > 0 ? 1.0f : -1.0f;
			pVertex0->z += pVertex0->vinfo.compensation.z > 0 ? 1.0f : -1.0f;
			pVertex0->vinfo.bAjusted = true;
		}
		if (!pVertex1->vinfo.bAjusted)
		{
			pVertex1->x += pVertex1->vinfo.compensation.x > 0 ? 1.0f : -1.0f;
			pVertex1->y += pVertex1->vinfo.compensation.y > 0 ? 1.0f : -1.0f;
			pVertex1->z += pVertex1->vinfo.compensation.z > 0 ? 1.0f : -1.0f;
			pVertex1->vinfo.bAjusted = true;
		}
		if (!pVertex2->vinfo.bAjusted)
		{
			pVertex2->x += pVertex2->vinfo.compensation.x > 0 ? 1.0f : -1.0f;
			pVertex2->y += pVertex2->vinfo.compensation.y > 0 ? 1.0f : -1.0f;
			pVertex2->z += pVertex2->vinfo.compensation.z > 0 ? 1.0f : -1.0f;
			pVertex2->vinfo.bAjusted = true;
		}
		if (!pVertex3->vinfo.bAjusted)
		{
			pVertex3->x += pVertex3->vinfo.compensation.x > 0 ? 1.0f : -1.0f;
			pVertex3->y += pVertex3->vinfo.compensation.y > 0 ? 1.0f : -1.0f;
			pVertex3->z += pVertex3->vinfo.compensation.z > 0 ? 1.0f : -1.0f;
			pVertex3->vinfo.bAjusted = true;
		}
	}
}

//void PolycubeMesh::GetMaxMinCoordinates()
//{
////	for (std::vector<unsigned long>::iterator iterPoint = m_vec_xyPlaneCorner.begin(); iterPoint != m_vec_xyPlaneCorner.end(); ++iterPoint)
////	{
////		GeoUtil::GetMaxMinCoordinateValueOfVertex(V[*iterPoint], m_maxVertex, m_minVertex);
////	}
//	for (std::vector<Vertex>::iterator iter = V.begin(); iter != V.end(); ++iter)
//	{
//		if (vertexInfo.at(std::distance(V.begin(), iter)).bSurface)
//		GeoUtil::GetMaxMinCoordinateValueOfVertex(*iter, m_maxVertex, m_minVertex);
//	}
//}

void PolycubeMesh::GetVertexInfo()
{
	cout << "m_vecXBoundaryLine size = " << m_vecXBoundaryLine.size() << std::endl;
	for (unsigned int i = 0; i < m_vecXBoundaryLine.size(); i++)
	{
		const Line& line = m_vecXBoundaryLine[i];
		V[line.p0].vinfo.bBoundaryLine = true;
		V[line.p1].vinfo.bBoundaryLine = true;
	}
	cout << "m_vecYBoundaryLine size = " << m_vecYBoundaryLine.size() << std::endl;
	for (unsigned int i = 0; i < m_vecYBoundaryLine.size(); i++)
	{
		const Line& line = m_vecYBoundaryLine[i];
		V[line.p0].vinfo.bBoundaryLine = true;
		V[line.p1].vinfo.bBoundaryLine = true;
	}
	cout << "m_vecZBoundaryLine size = " << m_vecZBoundaryLine.size() << std::endl;
	for (unsigned int i = 0; i < m_vecZBoundaryLine.size(); i++)
	{
		const Line& line = m_vecZBoundaryLine[i];
		V[line.p0].vinfo.bBoundaryLine = true;
		V[line.p1].vinfo.bBoundaryLine = true;
	}

	for (unsigned int i = 0; i < m_vec_xyPlaneCorner.size(); i++)
	{
		const unsigned long& p = m_vec_xyPlaneCorner[i];
		V[p].vinfo.bCorner = true;
	}
}

void PolycubeMesh::SmoothBoundaryLine(const unsigned int time)
{
	GetVertexInfo();
	int count = time;
	while (count-- != 0)
	{
		for (int i = 0; i < edgePatches.size(); i++)
		{
			std::vector<unsigned long> patchEdgeVertices;
			std::vector<Edge>& edge_patch = edgePatches.at(i);
			for (int k = 0; k < edge_patch.size(); k++)
			{
				const Edge& edge = edge_patch.at(k);
				patchEdgeVertices.push_back(edge.p1);
				patchEdgeVertices.push_back(edge.p0);
			}
			std::sort(patchEdgeVertices.begin(), patchEdgeVertices.end());
			std::vector<unsigned long>::iterator iter = std::unique(patchEdgeVertices.begin(), patchEdgeVertices.end());
			patchEdgeVertices.resize(std::distance(patchEdgeVertices.begin(), iter));
			//////////////////////////////////////////////////
			if ( i == 43)
			{
				std::cout << "edgepatch " << 43 << " vertex :";
				for (int n = 0; n < patchEdgeVertices.size(); n++)
				{
					std::cout << patchEdgeVertices.at(n) << ", ";
				}
				std::cout << std::endl;
			}
			/////////////////////////////////////////////////////
			for (int k = 0; k < patchEdgeVertices.size(); k++)
			{
				unsigned long index = patchEdgeVertices.at(k);
//				bool b = false;
//				if (index == 7590)
//				{
//					b = true;
//				}
				std::vector<unsigned long> jointConnectedVertices;
				if (V.at(index).vinfo.bCorner)
				{
					//std::cout << "Corner --------" << index << std::endl;
					continue;
				}
				for (int j = 0; j < edge_patch.size(); j++)
				{
					const Edge& edge = edge_patch.at(j);
					if (edge.p0 == index)
					{
						jointConnectedVertices.push_back(edge.p1);
					}
					else if (edge.p1 == index)
					{
						jointConnectedVertices.push_back(edge.p0);
					}
				}
				if (index == 7590)
				{
					std::cout << "7590(" << V.at(index).x << ", " << V.at(index).y << ", " << V.at(index).z << ")" << std::endl;
					std::cout << jointConnectedVertices.at(0) << "(" << V.at(jointConnectedVertices.at(0)).x << ", " << V.at(jointConnectedVertices.at(0)).y << ", " << V.at(jointConnectedVertices.at(0)).z << ")" << std::endl;
					std::cout << jointConnectedVertices.at(1) << "(" << V.at(jointConnectedVertices.at(1)).x << ", " << V.at(jointConnectedVertices.at(1)).y << ", " << V.at(jointConnectedVertices.at(1)).z << ")" << std::endl;
				}
				else if (index == 8163)
				{
					std::cout << "8163(" << V.at(index).x << ", " << V.at(index).y << ", " << V.at(index).z << ")" << std::endl;
					std::cout << jointConnectedVertices.at(0) << "(" << V.at(jointConnectedVertices.at(0)).x << ", " << V.at(jointConnectedVertices.at(0)).y << ", " << V.at(jointConnectedVertices.at(0)).z << ")" << std::endl;
					std::cout << jointConnectedVertices.at(1) << "(" << V.at(jointConnectedVertices.at(1)).x << ", " << V.at(jointConnectedVertices.at(1)).y << ", " << V.at(jointConnectedVertices.at(1)).z << ")" << std::endl;
				}
				V.at(index).x = (V[jointConnectedVertices.at(0)].x + V[jointConnectedVertices.at(1)].x) * 0.5;
				V.at(index).y = (V[jointConnectedVertices.at(0)].y + V[jointConnectedVertices.at(1)].y) * 0.5;
				V.at(index).z = (V[jointConnectedVertices.at(0)].z + V[jointConnectedVertices.at(1)].z) * 0.5;
				if (index == 7590)
				{
					std::cout << "7590(" << V.at(index).x << ", " << V.at(index).y << ", " << V.at(index).z << ")" << std::endl;
				}
				else if (index == 8163)
				{
					std::cout << "8163(" << V.at(index).x << ", " << V.at(index).y << ", " << V.at(index).z << ")" << std::endl;
				}
			}
		}

	}
}

void PolycubeMesh::SmoothBoundaryLine(const int patch, const unsigned int time)
{
	GetVertexInfo();
	int count = time;
	while (count-- != 0)
	{
		std::vector<unsigned long> patchEdgeVertices;
		std::vector<Edge>& edge_patch = edgePatches.at(patch);
		for (int i = 0; i < edge_patch.size(); i++)
		{
			std::vector<unsigned long> patchEdgeVertices;
			//std::vector<Edge>& edge_patch = edgePatches.at(i);
			const Edge& edge = edge_patch.at(i);
			V.at(edge.p1).vinfo.bBoundaryLine = true;
			V.at(edge.p0).vinfo.bBoundaryLine = true;
		}
		for (int k = 0; k < edge_patch.size(); k++)
		{
			const Edge& edge = edge_patch.at(k);
			patchEdgeVertices.push_back(edge.p1);
			patchEdgeVertices.push_back(edge.p0);
		}
		std::sort(patchEdgeVertices.begin(), patchEdgeVertices.end());
		std::vector<unsigned long>::iterator iter = std::unique(patchEdgeVertices.begin(), patchEdgeVertices.end());
		patchEdgeVertices.resize(std::distance(patchEdgeVertices.begin(), iter));

		for (int k = 0; k < patchEdgeVertices.size(); k++)
		{
			unsigned long index = patchEdgeVertices.at(k);
			std::vector<unsigned long> jointConnectedVertices;
			if (V.at(index).vinfo.bCorner)
				continue;
			for (int j = 0; j < edge_patch.size(); j++)
			{
				const Edge& edge = edge_patch.at(j);
				if (edge.p0 == index)
				{
					jointConnectedVertices.push_back(edge.p1);
				}
				else if (edge.p1 == index)
				{
					jointConnectedVertices.push_back(edge.p0);
				}
			}
			V.at(index).x = (V[jointConnectedVertices.at(0)].x + V[jointConnectedVertices.at(1)].x) * 0.5;
			V.at(index).y = (V[jointConnectedVertices.at(0)].y + V[jointConnectedVertices.at(1)].y) * 0.5;
			V.at(index).z = (V[jointConnectedVertices.at(0)].z + V[jointConnectedVertices.at(1)].z) * 0.5;
		}
	}
}
//void PolycubeMesh::SmoothBoundaryLine(const unsigned int time)
//{
//	GetVertexInfo();
//	int count = time;
//	while(count-- != 0)
//	{
//		std::map<unsigned long, glm::vec3> map_point_meanCoordinates;
//		std::multimap<unsigned long, Edge>::iterator iterBegin = V_E.begin();
//		std::multimap<unsigned long, Edge>::iterator iterEnd = V_E.end();
//		for (std::multimap<unsigned long, Edge>::iterator iter = iterBegin; iter != iterEnd;/* ++iter*/)
//		{
//			if (vertexInfo[iter->first].bBoundaryLine && !vertexInfo[iter->first].bCorner)
//			{
//				int pointCount = 0;
//				const unsigned long vertexIndex = iter->first;
//				glm::vec3 mean(0, 0, 0);
//				std::pair<std::multimap<unsigned long, Edge>::iterator, std::multimap<unsigned long, Edge>::iterator> ret = V_E.equal_range(iter->first);
//				for (iter = ret.first; iter != ret.second; ++iter)
//				{
//					const Edge& e = iter->second;
//					if (IsEdgeOnBondaryLine(e))
//					{
//						glm::vec3 m(0.0f, 0.0f, 0.0f);
//						GetCenterPointOfEdge(e, m);
//						mean += m;
//						pointCount++;
//					}
//				}
//				mean.x = mean.x/pointCount;
//				mean.y = mean.y/pointCount;
//				mean.z = mean.z/pointCount;
//				map_point_meanCoordinates.insert(std::pair<unsigned long, glm::vec3>(vertexIndex, mean));
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
//			Vertex& v = V[iter->first];
//			v.x = iter->second.x;
//			v.y = iter->second.y;
//			v.z = iter->second.z;
//		}
//	}
//}
void PolycubeMesh::CurvatureWeightingSmoothBoundaryLine(const unsigned int time)
{
	GetVertexInfo();
	int count = time;
	while(count-- != 0)
	{
		std::map<unsigned long, glm::vec3> map_point_meanCoordinates;
		std::multimap<unsigned long, Edge>::iterator iterBegin = V_E.begin();
		std::multimap<unsigned long, Edge>::iterator iterEnd = V_E.end();
		for (std::multimap<unsigned long, Edge>::iterator iter = iterBegin; iter != iterEnd;/* ++iter*/)
		{
			if (V[iter->first].vinfo.bBoundaryLine && !V[iter->first].vinfo.bCorner)
			{
				int pointCount = 0;
				const unsigned long vertexIndex = iter->first;
				int patch_label = vertexLabel.at(vertexIndex).at(0);
				glm::vec3 mean(0, 0, 0);
				std::pair<std::multimap<unsigned long, Edge>::iterator, std::multimap<unsigned long, Edge>::iterator> ret = V_E.equal_range(iter->first);
				for (iter = ret.first; iter != ret.second; ++iter)
				{
					const Edge& e = iter->second;
					//if (IsEdgeOnBondaryLine(e))
					if (IsEdgeOnPatch(e, patch_label))
					{
						glm::vec3 m(0.0f, 0.0f, 0.0f);
						GetCurvatureCenterPointOfEdge(e, m);
						mean += m;
						pointCount++;
					}
				}
				mean.x = mean.x/pointCount;
				mean.y = mean.y/pointCount;
				mean.z = mean.z/pointCount;
				map_point_meanCoordinates.insert(std::pair<unsigned long, glm::vec3>(vertexIndex, mean));
			}
			else
			{
				++iter;
			}
		}

		std::map<unsigned long, glm::vec3>::iterator iterB = map_point_meanCoordinates.begin();
		std::map<unsigned long, glm::vec3>::iterator iterE = map_point_meanCoordinates.end();
		for (std::map<unsigned long, glm::vec3>::iterator iter = iterB; iter != iterE; ++iter)
		{
			Vertex& v = V[iter->first];
			v.x = iter->second.x;
			v.y = iter->second.y;
			v.z = iter->second.z;
		}
	}
}

void PolycubeMesh::LaplacianWeightingSmoothBoundaryLine(const unsigned int time)
{
	GetVertexInfo();
	int count = time;
	while(count-- != 0)
	{
		std::map<unsigned long, glm::vec3> map_point_meanCoordinates;
		std::multimap<unsigned long, Edge>::iterator iterBegin = V_E.begin();
		std::multimap<unsigned long, Edge>::iterator iterEnd = V_E.end();
		for (std::multimap<unsigned long, Edge>::iterator iter = iterBegin; iter != iterEnd;/* ++iter*/)
		{
			if (V[iter->first].vinfo.bBoundaryLine && !V[iter->first].vinfo.bCorner)
			{
				int pointCount = 0;
				const unsigned long vertexIndex = iter->first;
				int patch_label = vertexLabel.at(vertexIndex).at(0);
				glm::vec3 mean(0, 0, 0);
				std::pair<std::multimap<unsigned long, Edge>::iterator, std::multimap<unsigned long, Edge>::iterator> ret = V_E.equal_range(iter->first);
				for (iter = ret.first; iter != ret.second; ++iter)
				{
					const Edge& e = iter->second;
					if (IsEdgeOnBondaryLine(e))
					//if (IsEdgeOnPatch(e, patch_label))
					{
						glm::vec3 m(0.0f, 0.0f, 0.0f);
						GetCenterPointOfEdge(e, m);
						mean += m;
						pointCount++;
					}
				}
				mean.x = mean.x/pointCount;
				mean.y = mean.y/pointCount;
				mean.z = mean.z/pointCount;
				map_point_meanCoordinates.insert(std::pair<unsigned long, glm::vec3>(vertexIndex, mean));
			}
			else
			{
				++iter;
			}
		}

		std::map<unsigned long, glm::vec3>::iterator iterB = map_point_meanCoordinates.begin();
		std::map<unsigned long, glm::vec3>::iterator iterE = map_point_meanCoordinates.end();
		for (std::map<unsigned long, glm::vec3>::iterator iter = iterB; iter != iterE; ++iter)
		{
			Vertex& v = V[iter->first];
			v.x = iter->second.x;
			v.y = iter->second.y;
			v.z = iter->second.z;
		}
	}
}

void PolycubeMesh::SmoothSurface(const unsigned int time)
{
	// 根据patches smooth
	GetVertexInfo();
	int count = time;
	while(count-- != 0)
	{
		std::map<unsigned long, glm::vec3> map_point_meanCoordinates;
		std::multimap<unsigned long, SFace>::iterator iterBegin = V_F.begin();
		std::multimap<unsigned long, SFace>::iterator iterEnd = V_F.end();
		for (std::multimap<unsigned long, SFace>::iterator iter = iterBegin; iter != iterEnd;/* ++iter*/)
		{
			if (V[iter->first].vinfo.bSurface && !V[iter->first].vinfo.bBoundaryLine)
			{
				int pointCount = 0;
				const unsigned long vertexIndex = iter->first;
				int patch_label = vertexLabel.at(vertexIndex).at(0);
				glm::vec3 mean(0, 0, 0);
				std::pair<std::multimap<unsigned long, SFace>::iterator, std::multimap<unsigned long, SFace>::iterator> ret = V_F.equal_range(iter->first);
				for (iter = ret.first; iter != ret.second; ++iter)
				{
					const Face& f = iter->second;
					if (IsFaceOnPatch(f, patch_label))
					{
						glm::vec3 m(0.0f, 0.0f, 0.0f);
						GetCenterPointOfCell(f, m);
						mean += m;
						pointCount++;
//						if (vertexIndex == 2143)
//						{
//							const unsigned long vertexIndex0 = f.at(0);
//							const unsigned long vertexIndex1 = f.at(1);
//							const unsigned long vertexIndex2 = f.at(2);
//							cout << "V[" << vertexIndex0 << "] = (" << V[vertexIndex0].x << ", " << V[vertexIndex0].y << ", " << V[vertexIndex0].z << ")" << std::endl;
//							cout << "V[" << vertexIndex1 << "] = (" << V[vertexIndex1].x << ", " << V[vertexIndex1].y << ", " << V[vertexIndex1].z << ")" << std::endl;
//							cout << "V[" << vertexIndex2 << "] = (" << V[vertexIndex2].x << ", " << V[vertexIndex2].y << ", " << V[vertexIndex2].z << ")" << std::endl;
//
//							cout << "V[m] = (" << m.x << ", " << m.y << ", " << m.z << ")" << std::endl;
//						}
					}
				}
				mean.x = mean.x/pointCount;
				mean.y = mean.y/pointCount;
				mean.z = mean.z/pointCount;
				map_point_meanCoordinates.insert(std::pair<unsigned long, glm::vec3>(vertexIndex, mean));
//				if (vertexIndex == 2143)
//				{
//					cout << "V[2143] = (" << V[vertexIndex].x << ", " << V[vertexIndex].y << ", " << V[vertexIndex].z << ")" << std::endl;
//					cout << "V[mean] = (" << mean.x << ", " << mean.y << ", " << mean.z << ")" << std::endl;
//				}
			}
			else
			{
				++iter;
			}
		}

		std::map<unsigned long, glm::vec3>::iterator iterB = map_point_meanCoordinates.begin();
		std::map<unsigned long, glm::vec3>::iterator iterE = map_point_meanCoordinates.end();
		for (std::map<unsigned long, glm::vec3>::iterator iter = iterB; iter != iterE; ++iter)
		{
			Vertex& v = V[iter->first];
			v.x = iter->second.x;
			v.y = iter->second.y;
			v.z = iter->second.z;
		}
	}
}

void PolycubeMesh::SmoothCurvature(const unsigned int time)
{
	GetVertexInfo();
	int count = time;
	while(count-- != 0)
	{
		std::map<unsigned long, double> map_point_meanCurvature;
		std::multimap<unsigned long, SFace>::iterator iterBegin = V_F.begin();
		std::multimap<unsigned long, SFace>::iterator iterEnd = V_F.end();
		for (std::multimap<unsigned long, SFace>::iterator iter = iterBegin; iter != iterEnd;/* ++iter*/)
		{
			if (V[iter->first].vinfo.bSurface && !V[iter->first].vinfo.bBoundaryLine)
			{
				int pointCount = 0;
				const unsigned long vertexIndex = iter->first;
				double mean = 0;
				std::pair<std::multimap<unsigned long, SFace>::iterator, std::multimap<unsigned long, SFace>::iterator> ret = V_F.equal_range(iter->first);
				for (iter = ret.first; iter != ret.second; ++iter)
				{
					const Face& f = iter->second;
					if (IsFaceOnSurface(f))
					{
						mean += GetCenterCurvatureOfCell(f);
						pointCount++;
					}
				}
				mean = mean/pointCount;
				map_point_meanCurvature.insert(std::pair<unsigned long, double>(vertexIndex, mean));
			}
			else
			{
				++iter;
			}
		}

		std::map<unsigned long, double>::iterator iterB = map_point_meanCurvature.begin();
		std::map<unsigned long, double>::iterator iterE = map_point_meanCurvature.end();
		for (std::map<unsigned long, double>::iterator iter = iterB; iter != iterE; ++iter)
		{
			V[iter->first].vinfo.curvature = iter->second;
		}
	}
	char filename[64] = {0};
	sprintf(filename, "smoothcurvature%d.txt", time);
	std::ofstream sc(filename);
	for (unsigned long i = 0; i < V.size(); ++i)
	{
		if (V[i].vinfo.bSurface)
		{
			sc << V[i].vinfo.curvature << std::endl;
		}
		else
		{
			sc << 0 << std::endl;
		}
	}
	sc.close();
}

void PolycubeMesh::CurvatureWeightingSmoothSurface(const unsigned int time)
{
	GetVertexInfo();
	int count = time;
	while(count-- != 0)
	{
		std::map<unsigned long, glm::vec3> map_point_meanCoordinates;
		std::multimap<unsigned long, SFace>::iterator iterBegin = V_F.begin();
		std::multimap<unsigned long, SFace>::iterator iterEnd = V_F.end();
		for (std::multimap<unsigned long, SFace>::iterator iter = iterBegin; iter != iterEnd;/* ++iter*/)
		{
			if (V[iter->first].vinfo.bSurface && !V[iter->first].vinfo.bBoundaryLine)
			{
				int pointCount = 0;
				const unsigned long vertexIndex = iter->first;
				glm::vec3 mean(0, 0, 0);
				std::pair<std::multimap<unsigned long, SFace>::iterator, std::multimap<unsigned long, SFace>::iterator> ret = V_F.equal_range(iter->first);
				for (iter = ret.first; iter != ret.second; ++iter)
				{
					const Face& f = iter->second;
					if (IsFaceOnSurface(f))
					{
						glm::vec3 m(0.0f, 0.0f, 0.0f);
						GetCurvatureCenterPointOfCell(f, m);
						mean += m;
						pointCount++;
					}
				}
				mean.x = mean.x/pointCount;
				mean.y = mean.y/pointCount;
				mean.z = mean.z/pointCount;
				map_point_meanCoordinates.insert(std::pair<unsigned long, glm::vec3>(vertexIndex, mean));
			}
			else
			{
				++iter;
			}
		}

		std::map<unsigned long, glm::vec3>::iterator iterB = map_point_meanCoordinates.begin();
		std::map<unsigned long, glm::vec3>::iterator iterE = map_point_meanCoordinates.end();
		for (std::map<unsigned long, glm::vec3>::iterator iter = iterB; iter != iterE; ++iter)
		{
			Vertex& v = V[iter->first];
			v.x = iter->second.x;
			v.y = iter->second.y;
			v.z = iter->second.z;
		}
	}
}

bool PolycubeMesh::IsVerticesOnTheSamePlane(const std::vector<unsigned long>& surfaceVertices, DIRECTION& d)
{
	DIRECTION dir = UNKNOWN_AXIS;
	glm::vec3 v(0.0f, 0.0f, 0.0f);
	unsigned long index = 0;
	for (unsigned long i = 0; i < surfaceVertices.size(); i++)
	{
		if (V[i].vinfo.bSurface && !V[i].vinfo.bBoundaryLine)
		{
			dir = V[i].vinfo.dir;
			v.x = V[i].x;
			v.y = V[i].y;
			v.z = V[i].z;
			index = i;
			break;
		}
	}

	if (dir != UNKNOWN_AXIS)
	{
		d = dir;
		std::pair<std::multimap<unsigned long, SFace>::iterator, std::multimap<unsigned long, SFace>::iterator > ret = V_F.equal_range(index);
		for (std::multimap<unsigned long, SFace>::iterator iter = ret.first; iter != ret.second; ++iter)
		{
			const Face& face = iter->second;
			if (IsFaceOnSurface(face))
			{
				const Vertex* pVertex0 = &V[face.at(0)];
				const Vertex* pVertex1 = &V[face.at(1)];
				const Vertex* pVertex2 = &V[face.at(2)];

				const glm::vec3 v0(pVertex0->x, pVertex0->y, pVertex0->z);
				const glm::vec3 v1(pVertex1->x, pVertex1->y, pVertex1->z);
				const glm::vec3 v2(pVertex2->x, pVertex2->y, pVertex2->z);

				const glm::vec3 v10 = v0 - v1;
				const glm::vec3 v12 = v2 - v1;
				const glm::vec3 normal = glm::cross(v12, v10);
				if (dir == X_AXIS){
					if (normal.x > 0) { dir = X_AXIS_PLUS;}
					else{ dir = X_AXIS_MINUS;}
				}
				else if (dir == Y_AXIS)	{
					if (normal.y > 0) { dir = Y_AXIS_PLUS;}
					else{dir = Y_AXIS_MINUS;}
				}
				else if (dir == Z_AXIS)	{
					if (normal.z > 0) { dir = Z_AXIS_PLUS;}
					else { dir = Z_AXIS_MINUS;}
				}
				break;
			}
		}
		for (unsigned long i = 0; i < surfaceVertices.size(); i++)
		{
			if (V[i].vinfo.bSurface && !V[i].vinfo.bBoundaryLine)
			{
				if (dir == X_AXIS)
				{
					if (fabs(v.x - V[i].x) > 1e-2)
					{
						return false;
					}
				}
				else if (dir == Y_AXIS)
				{
					if (fabs(v.y - V[i].y) > 1e-2)
					{
						return false;
					}
				}
				else if (dir == Z_AXIS)
				{
					if (fabs(v.z - V[i].z) > 1e-2)
					{
						return false;
					}
				}
			}
		}
	}

	return true;
}

struct DirVec
{
	DIRECTION dir;
	glm::vec3 vec;
};

const DirVec dirVec[6] =
{
	{X_AXIS_PLUS, glm::vec3(1.0, 0.0, 0.0)},
	{Y_AXIS_PLUS, glm::vec3(0.0, 1.0, 0.0)},
	{Z_AXIS_PLUS, glm::vec3(0.0, 0.0, 1.0)},
	{X_AXIS_MINUS, glm::vec3(-1.0, 0.0, 0.0)},
	{Y_AXIS_MINUS, glm::vec3(0.0, -1.0, 0.0)},
	{Z_AXIS_MINUS, glm::vec3(0.0, 0.0, -1.0)}
};

void PolycubeMesh::SmoothBoundaryCell(const unsigned int time)
{
	GetVertexInfo();
	int count = time;
	while(count-- != 0)
	{
		//std::map<unsigned long, glm::vec3> map_point_meanCoordinates;
		std::multimap<unsigned long, Edge>::iterator iterBegin = V_E.begin();
		std::multimap<unsigned long, Edge>::iterator iterEnd = V_E.end();
		for (std::multimap<unsigned long, Edge>::iterator iter = iterBegin; iter != iterEnd;/* ++iter*/)
		{
			if (!V[iter->first].vinfo.bSurface && V[iter->first].vinfo.bBoundaryCell)
			{
				int pointCount = 0;
				const unsigned long vertexIndex = iter->first;
				glm::vec3 mean(0, 0, 0);
				std::vector<unsigned long> surfaceVertices;
				std::pair<std::multimap<unsigned long, Edge>::iterator, std::multimap<unsigned long, Edge>::iterator> ret = V_E.equal_range(iter->first);
				for (iter = ret.first; iter != ret.second; ++iter)
				{
					const Edge& e = iter->second;
					if (V[e.p0].vinfo.bSurface)
					{
						surfaceVertices.push_back(e.p0);
					}
					else if (V[e.p1].vinfo.bSurface)
					{
						surfaceVertices.push_back(e.p1);
					}
				}
				DIRECTION dir = UNKNOWN_AXIS;
				if (IsVerticesOnTheSamePlane(surfaceVertices, dir))
				{
					GetCenterPointOfCell(surfaceVertices, mean);
					if (dir != UNKNOWN_AXIS)
					{
						Vertex& v = V[iter->first];
						glm::vec3 orig_pos(v.x, v.y, v.z);
						glm::vec3 unit_normal = glm::vec3(0.0, 0.0, 0.0) - dirVec[dir - X_AXIS_PLUS].vec;
						const double s = glm::dot(orig_pos - mean, unit_normal);
						glm::vec3 new_pos = glm::vec3(s*unit_normal.x, s*unit_normal.y, s*unit_normal.z) + mean;
						v.x = new_pos.x;
						v.y = new_pos.y;
						v.z = new_pos.z;
					}
				}
//				m_map_point_meanCoordinates.insert(std::pair<unsigned long, glm::vec3>(vertexIndex, mean));
			}
			else
			{
				++iter;
			}
		}

//		std::map<unsigned long, glm::vec3>::iterator iterB = m_map_point_meanCoordinates.begin();
//		std::map<unsigned long, glm::vec3>::iterator iterE = m_map_point_meanCoordinates.end();
//		for (std::map<unsigned long, glm::vec3>::iterator iter = iterB; iter != iterE; ++iter)
//		{
//			Vertex& v = V[iter->first];
//			v.x = iter->second.x;
//			v.y = iter->second.y;
//			v.z = iter->second.z;
//		}
	}
}

void PolycubeMesh::SmoothVolume(const unsigned int time)
{
	//GetVertexInfo();
	int count = time;
	while(count-- != 0)
	{
		std::map<unsigned long, glm::vec3> map_point_meanCoordinates;
		std::multimap<unsigned long, unsigned long>::iterator iterBegin = VI_CI.begin();
		std::multimap<unsigned long, unsigned long>::iterator iterEnd = VI_CI.end();
		for (std::multimap<unsigned long, unsigned long>::iterator iter = iterBegin; iter != iterEnd;/* ++iter*/)
		{
			unsigned long index = iter->first;
			if (!V[iter->first].vinfo.bSurface && V[iter->first].vinfo.bBoundaryCell)
			{
				int pointCount = 0;
				glm::vec3 mean(0, 0, 0);
				std::pair<std::multimap<unsigned long, unsigned long>::iterator, std::multimap<unsigned long, unsigned long>::iterator> ret = VI_CI.equal_range(iter->first);
				for (iter = ret.first; iter != ret.second; ++iter)
				{
					const Cell& c = C.at(iter->second);
					glm::vec3 m(0.0f, 0.0f, 0.0f);
					GetCenterPointOfCell(c, m);
					mean += m;
					pointCount++;
//					if (index == 25905)
//					{
//						std::cout << " point 25905 connected cell centroid :(" << m.x << ", "<< m.y << ", "<< m.z << ")"<< std::endl;
//					}
				}
				mean.x = mean.x/pointCount;
				mean.y = mean.y/pointCount;
				mean.z = mean.z/pointCount;
//				if (index == 25905)
//				{
//					const Vertex& v = V.at(index);
//					std::cout << " point 25905 orig :(" << v.x  << ", "<< v.y << ", " << v.z << ")" << std::endl;
//					std::cout << " point 25905 mean :(" << mean.x  << ", "<< mean.y << ", " << mean.z << ")" << std::endl;
//				}
				map_point_meanCoordinates.insert(std::pair<unsigned long, glm::vec3>(index, mean));
			}
			else
			{
				++iter;
			}
		}

		std::map<unsigned long, glm::vec3>::iterator iterB = map_point_meanCoordinates.begin();
		std::map<unsigned long, glm::vec3>::iterator iterE = map_point_meanCoordinates.end();
		for (std::map<unsigned long, glm::vec3>::iterator iter = iterB; iter != iterE; ++iter)
		{
			if (!V[iter->first].vinfo.bSurface){
			Vertex& v = V[iter->first];
			v.x = iter->second.x;
			v.y = iter->second.y;
			v.z = iter->second.z;}
		}
	}
}

//void PolycubeMesh::SmoothVolume(const unsigned int time)
//{
//	GetVertexInfo();
//	int count = time;
//	while(count-- != 0)
//	{
//		std::map<unsigned long, glm::vec3> map_point_meanCoordinates;
//		std::multimap<unsigned long, Cell>::iterator iterBegin = V_C.begin();
//		std::multimap<unsigned long, Cell>::iterator iterEnd = V_C.end();
//		for (std::multimap<unsigned long, Cell>::iterator iter = iterBegin; iter != iterEnd;/* ++iter*/)
//		{
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
//				}
//				mean.x = mean.x/pointCount;
//				mean.y = mean.y/pointCount;
//				mean.z = mean.z/pointCount;
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

void PolycubeMesh::Smooth(const unsigned int time)
{
	SmoothBoundaryLine(time);
	SmoothSurface(time);
	SmoothVolume(time);
}

void PolycubeMesh::CurvatureWeightingSmooth(const unsigned int time)
{
	CurvatureWeightingSmoothBoundaryLine(time);
	CurvatureWeightingSmoothSurface(time);
	//CurvatureWeightingSmoothVolume(time);
}

void PolycubeMesh::AjustSurfaceToIntegerPlane(const float delta)
{
	for (std::vector<Face>::iterator iterSurface = surface.begin(); iterSurface != surface.end(); ++iterSurface)
	{
		const Face& face = *iterSurface;
		const size_t faceSize = face.size();
		glm::vec3 center(0.0f, 0.0f, 0.0f);
		GetCenterPointOfCell(face, center);

		const Vertex& v0 = V[iterSurface->at(0)];
		const Vertex& v1 = V[iterSurface->at(1)];
		const Vertex& v2 = V[iterSurface->at(2)];
		const Vector p0p1(v1.x - v0.x, v1.y - v0.y, v1.z - v0.z);
		const Vector p0p2(v2.x - v0.x, v2.y - v0.y, v2.z - v0.z);
		const Vector normal(p0p1.j * p0p2.k - p0p1.k * p0p2.j,
				p0p1.k * p0p2.i - p0p1.i * p0p2.k,
				p0p1.i * p0p2.j - p0p1.j * p0p2.i);

		if (fabs(normal.i) > 0.9)
		{
			for (unsigned int i = 0; i < faceSize; i++)
			{
				V[face.at(i)].x = int(center.x/delta)*delta;
			}
		}
		else if (fabs(normal.j) > 0.9)
		{
			for (unsigned int i = 0; i < faceSize; i++)
			{
				V[face.at(i)].y = int(center.y/delta)*delta;
			}
		}
		else if (fabs(normal.k) > 0.9)
		{
			for (unsigned int i = 0; i < faceSize; i++)
			{
				V[face.at(i)].z = int(center.z/delta)*delta;
			}
		}
	}
}

const double alpha = 1e-6;
void PolycubeMesh::GetE_ij_and_F_ij(const unsigned long i, const unsigned long j, double& E_ij, glm::vec3& F_ij)
{
	const glm::vec3 Vi(V[i].x, V[i].y, V[i].z);
	const glm::vec3 Vj(V[j].x, V[j].y, V[j].z);
	const glm::vec3 Vij(Vi - Vj);
	const double r = glm::length(Vij);
	const double sigma = V[j].vinfo.curvature;
//	double E_ij = sigma/tan(0.5*PI*r/sigma) + 0.5*PI*r/sigma - 0.5*PI;
	E_ij = alpha*sigma/(r*r);
	F_ij.x = E_ij*Vij.x;
	F_ij.y = E_ij*Vij.y;
	F_ij.z = E_ij*Vij.z;
	return;
}

void PolycubeMesh::GetSurfaceNeignbor(const unsigned long vertexIndex, std::vector<unsigned long>& vecNeighbor)
{
	std::pair<std::multimap<unsigned long, Edge>::iterator, std::multimap<unsigned long, Edge>::iterator> ret =
			V_E.equal_range(vertexIndex);
	for (std::multimap<unsigned long, Edge>::iterator iter = ret.first; iter != ret.second; ++iter)
	{
		if (V.at(vertexIndex).vinfo.bSurface)
		{
			const Edge& e = iter->second;
			if (e.p0 == vertexIndex)
			{
				vecNeighbor.push_back(e.p1);
			}
			else
			{
				vecNeighbor.push_back(e.p0);
			}
		}
	}
	std::sort(vecNeighbor.begin(), vecNeighbor.end());
	std::vector<unsigned long>::iterator iter = std::unique(vecNeighbor.begin(), vecNeighbor.end());
	vecNeighbor.resize(std::distance(vecNeighbor.begin(), iter));
}

void PolycubeMesh::GetBoundaryLineNeignbor(const unsigned long vertexIndex, std::vector<unsigned long>& vecNeighbor)
{
	std::pair<std::multimap<unsigned long, Edge>::iterator, std::multimap<unsigned long, Edge>::iterator> ret =
			V_E.equal_range(vertexIndex);
	for (std::multimap<unsigned long, Edge>::iterator iter = ret.first; iter != ret.second; ++iter)
	{
		const Edge& e = iter->second;
		if (IsEdgeOnBondaryLine(e))
		{
			unsigned long i = vertexIndex;
			unsigned long j = e.p0;
			if (e.p0 == vertexIndex)
			{
				vecNeighbor.push_back(e.p1);
			}
			else
			{
				vecNeighbor.push_back(e.p0);
			}
		}
	}
	std::sort(vecNeighbor.begin(), vecNeighbor.end());
	std::vector<unsigned long>::iterator iter = std::unique(vecNeighbor.begin(), vecNeighbor.end());
	vecNeighbor.resize(std::distance(vecNeighbor.begin(), iter));
}

void PolycubeMesh::ParticleSmooth(const unsigned int time)
{
	unsigned int count = time;
	while (count-- != 0)
	{
		const size_t n = V.size();
		std::vector<glm::vec3> V_V;
		for (unsigned long i = 0; i < n; i++)
		{
			if (V[i].vinfo.bSurface && !V[i].vinfo.bBoundaryLine)
			{
				std::vector<unsigned long> vecNeighbor;
				GetSurfaceNeignbor(i, vecNeighbor);
				double E_i = 0;
				glm::vec3 F_i(0.0, 0.0, 0.0);
				for (unsigned long j = 0; j < vecNeighbor.size(); j++)
				{
					double E_ij = 0.0f;
					glm::vec3 F_ij(0.0, 0.0, 0.0);
					GetE_ij_and_F_ij(i, vecNeighbor[j], E_ij, F_ij);
					E_i += E_ij;
					F_i += F_ij;
				}
				double minDistance = 1000000;
				glm::vec3 vi(V[i].x, V[i].y, V[i].z);
				for (unsigned long j = 0; j < vecNeighbor.size(); j++)
				{
					unsigned long _j = vecNeighbor[j];
					glm::vec3 vj(V[_j].x, V[_j].y, V[_j].z);
				    double d_ij = glm::length(vi - vj);
				    if (d_ij < minDistance)
				    	minDistance = d_ij;
				}
				double displacement = glm::length(F_i);
				if (displacement > minDistance)
				{
					F_i.x = F_i.x*minDistance/displacement;
					F_i.y = F_i.y*minDistance/displacement;
					F_i.z = F_i.z*minDistance/displacement;
				}

				V[i].x += 0.166*F_i.x;
				V[i].y += 0.166*F_i.y;
				V[i].z += 0.166*F_i.z;
				//V_V.push_back(F_i);
			}
			else if (V[i].vinfo.bBoundaryLine && !V[i].vinfo.bCorner)
			{
				std::vector<unsigned long> vecNeighbor;
				GetBoundaryLineNeignbor(i, vecNeighbor);
				double E_i = 0;
				glm::vec3 F_i(0.0, 0.0, 0.0);
				for (unsigned long j = 0; j < vecNeighbor.size(); j++)
				{
					double E_ij = 0.0f;
					glm::vec3 F_ij(0.0, 0.0, 0.0);
					GetE_ij_and_F_ij(i, vecNeighbor[j], E_ij, F_ij);
					E_i += E_ij;
					F_i += F_ij;
				}
				double minDistance = 1000000;
				glm::vec3 vi(V[i].x, V[i].y, V[i].z);
				for (unsigned long j = 0; j < vecNeighbor.size(); j++)
				{
					unsigned long _j = vecNeighbor[j];
					glm::vec3 vj(V[_j].x, V[_j].y, V[_j].z);
				    double d_ij = glm::length(vi - vj);
				    if (d_ij < minDistance)
				    	minDistance = d_ij;
				}
				double displacement = glm::length(F_i);
				if (displacement > minDistance)
				{
					F_i.x = F_i.x*minDistance/displacement;
					F_i.y = F_i.y*minDistance/displacement;
					F_i.z = F_i.z*minDistance/displacement;
				}

				V[i].x += 0.166*F_i.x;
				V[i].y += 0.166*F_i.y;
				V[i].z += 0.166*F_i.z;
				//V_V.push_back(F_i);
			}
			else
			{
				;//V_V.push_back(glm::vec3(0.0, 0.0, 0.0));
			}
		}
	}
}

void PolycubeMesh::LabelSurfaceFace()
{
	faceLabel.clear();
	faceLabel.resize(surface.size(), 0);
	int currentLabel = 1;
	LabelSinglePlaneFace(m_xSinglePlane_Face_Multimap, currentLabel); m_xLabel = currentLabel - 1;
	LabelSinglePlaneFace(m_ySinglePlane_Face_Multimap, currentLabel); m_yLabel = currentLabel - 1;
	LabelSinglePlaneFace(m_zSinglePlane_Face_Multimap, currentLabel); m_zLabel = currentLabel - 1;
	m_patchNumber = currentLabel - 1;

	std::cout << " x patches: " << m_xLabel  << ": ( 1 - " << m_xLabel << " )" << std::endl;
	std::cout << " y patches: " << m_yLabel  << ": ( " << m_xLabel + 1 << " - " << m_yLabel << " )" << std::endl;
	std::cout << " z patches: " << m_zLabel  << ": ( " << m_yLabel + 1 << " - " << m_zLabel << " )" << std::endl;
}

void PolycubeMesh::LabelSinglePlaneFace(const df& singlePlane_Face_Multimap, int& label)
{
	const df::const_iterator iterSinglePlaneBegin = singlePlane_Face_Multimap.begin();
	const df::const_iterator iterSinglePlaneEnd = singlePlane_Face_Multimap.end();
	df::const_iterator iterSinglePlane = iterSinglePlaneBegin;
	while (iterSinglePlane != iterSinglePlaneEnd)
	{
		std::vector<Face*> vecSinglePlaneFace;
		std::pair<df::const_iterator, df::const_iterator> ret = singlePlane_Face_Multimap.equal_range(iterSinglePlane->first);
		int count = 0;
		for (iterSinglePlane = ret.first; iterSinglePlane != ret.second; ++iterSinglePlane)
		{
			const size_t PF_SIZE = surface.size();
			for (int i = 0; i < PF_SIZE; i++)
			{
				if (iterSinglePlane->second == &surface.at(i))
				{
					faceLabel.at(i) = label;
					count++;
					//std::cout << "i = " << i << " label = " << label << std::endl;
					break;
				}
			}
		}
#ifndef _TEST
		label++;
#else
		if (count > 2)
			label++;
		else
		{
			for (iterSinglePlane = ret.first; iterSinglePlane != ret.second; ++iterSinglePlane)
			{
				const size_t PF_SIZE = surface.size();
				for (int i = 0; i < PF_SIZE; i++)
				{
					if (iterSinglePlane->second == &surface.at(i))
					{
						faceLabel.at(i) = label - 1;
						break;
					}
				}
			}
		}
#endif
	}
}
void PolycubeMesh::GetAreaPatches()
{
	if (faceLabel.empty())
	{
		std::cout << "faceLabel is empty" << std::endl;
		GetPatchesLabel();
	}
	facePatches.clear();
	facePatches.resize(m_patchNumber);
	const size_t PF_SIZE = surface.size();
	for (int i = 0; i < PF_SIZE; i++)
	{
		if (faceLabel.at(i) != 0)
		{
			unsigned long facePatchIndex = faceLabel.at(i) - 1;
			facePatches.at(facePatchIndex).push_back(i);
		}
	}
}

void PolycubeMesh::GetFacePatches()
{
	facePatches.clear();
	facePatches.resize(m_patchNumber);
	const size_t PF_SIZE = surface.size();
	for (int i = 0; i < PF_SIZE; i++)
	{
		if (faceLabel.at(i) != 0)
		{
			unsigned long facePatchIndex = faceLabel.at(i) - 1;
			facePatches.at(facePatchIndex).push_back(i);
		}
	}
}

void PolycubeMesh::GetEdgePatches()
{
	edgePatches.clear();
//	edgePatches.resize(faceLable.size());
	const size_t facePatches_SIZE = facePatches.size();
	for (int i = 0; i < facePatches_SIZE; i++)
	{
		std::vector<Edge> edgePatch;
		GetBoundaryEdge(facePatches.at(i), edgePatch);
		edgePatches.push_back(edgePatch);
	}
}

void PolycubeMesh::GetBoundaryEdge(const std::vector<unsigned long>& faceIndex, std::vector<Edge>& edgePatch)
{
	for (int i = 0; i < faceIndex.size(); i++)
	{
		const Cell& face = surface.at(faceIndex.at(i));
		if (face.size() == 3)
		{
			Edge e1(face.at(0), face.at(1));
			Edge e2(face.at(1), face.at(2));
			Edge e3(face.at(2), face.at(0));
			std::vector<Edge>::iterator iter = std::find(edgePatch.begin(), edgePatch.end(), e1);
			if (iter != edgePatch.end()) iter = edgePatch.erase(iter);
			else edgePatch.push_back(e1);

			iter = std::find(edgePatch.begin(), edgePatch.end(), e2);
			if (iter != edgePatch.end()) iter = edgePatch.erase(iter);
			else edgePatch.push_back(e2);

			iter = std::find(edgePatch.begin(), edgePatch.end(), e3);
			if (iter != edgePatch.end()) iter = edgePatch.erase(iter);
			else edgePatch.push_back(e3);
		}
		if (face.size() == 4)
		{
			Edge e1(face.at(0), face.at(1));
			Edge e2(face.at(1), face.at(2));
			Edge e3(face.at(2), face.at(3));
			Edge e4(face.at(3), face.at(0));
			std::vector<Edge>::iterator iter = std::find(edgePatch.begin(), edgePatch.end(), e1);
			if (iter != edgePatch.end()) iter = edgePatch.erase(iter);
			else edgePatch.push_back(e1);

			iter = std::find(edgePatch.begin(), edgePatch.end(), e2);
			if (iter != edgePatch.end()) iter = edgePatch.erase(iter);
			else edgePatch.push_back(e2);

			iter = std::find(edgePatch.begin(), edgePatch.end(), e3);
			if (iter != edgePatch.end()) iter = edgePatch.erase(iter);
			else edgePatch.push_back(e3);

			iter = std::find(edgePatch.begin(), edgePatch.end(), e4);
			if (iter != edgePatch.end()) iter = edgePatch.erase(iter);
			else edgePatch.push_back(e4);
		}
	}
}

void PolycubeMesh::GetVertexPatches()
{
//	GetSurfaceVertexIndices();
//	vertexPatches.clear();
//	vertexPatches.resize(surfaceVertexIndices.size());

	vertexPatches.clear();
	vertexPatches.resize(edgePatches.size());

	vertexLabel.clear();
	vertexLabel.resize(V.size());

	const size_t Patches_SIZE = facePatches.size();
	// for each patch
	for (int i = 0; i < Patches_SIZE; i++)
	{
		const std::vector<unsigned long>& facePatch = facePatches.at(i);
		for (int j = 0; j < facePatch.size(); j++)
		{
			const Face& face = surface.at(facePatch.at(j));
			vertexPatches.at(i).push_back(face.at(0));
			vertexPatches.at(i).push_back(face.at(1));
			vertexPatches.at(i).push_back(face.at(2));
			if (face.size() == 4)
				vertexPatches.at(i).push_back(face.at(3));

			vertexLabel.at(face.at(0)).push_back(i + 1);
			vertexLabel.at(face.at(1)).push_back(i + 1);
			vertexLabel.at(face.at(2)).push_back(i + 1);
			if (face.size() == 4)
				vertexLabel.at(face.at(3)).push_back(i + 1);
		}
	}
	for (int i = 0; i < Patches_SIZE; i++)
	{
		std::vector<unsigned long>& vertexPatch = vertexPatches.at(i);
		std::sort(vertexPatch.begin(), vertexPatch.end());
		std::vector<unsigned long>::iterator iter = std::unique(vertexPatch.begin(), vertexPatch.end());
		size_t size = std::distance(vertexPatch.begin(), iter);
		vertexPatch.resize(size);
	}

	for (int i = 0; i < V.size(); i++)
	{
		if (V.at(i).vinfo.bSurface)
		{
			std::vector<int>& vertex_label = vertexLabel.at(i);
			std::sort(vertex_label.begin(), vertex_label.end());
			std::vector<int>::iterator iter = std::unique(vertex_label.begin(), vertex_label.end());
			size_t size = std::distance(vertex_label.begin(), iter);
			vertex_label.resize(size);
		}
	}
}

int PolycubeMesh::DividePatch(const std::vector<unsigned long>& vertexPatch,
		std::vector<std::vector<unsigned long> >& vertex_patches, const DIRECTION dir)
{
	if (vertexPatch.empty() || vertexPatch.size() < 3)
		return 0;

	std::vector<unsigned long> connectedPatch;
	connectedPatch.push_back(vertexPatch.at(0));
	size_t currentSize = 0;
	while(currentSize != connectedPatch.size())
	{
		currentSize = connectedPatch.size();
		std::vector<unsigned long> connectedRegion = connectedPatch;
		for (int j = 0; j < connectedPatch.size(); j++)
		{
			unsigned long vertexIndex = connectedPatch.at(j);
			const Vertex& vertex = V.at(vertexIndex);
			std::pair<std::multimap<unsigned long, unsigned long>::iterator, std::multimap<unsigned long, unsigned long>::iterator>
			ret = V_V.equal_range(vertexIndex);
			for (std::multimap<unsigned long, unsigned long>::iterator iter = ret.first; iter != ret.second; ++iter)
			{
				unsigned long connectedVertexIndex = iter->second;
				if (!V.at(connectedVertexIndex).vinfo.bSurface)
					continue;
				std::vector<unsigned long>::const_iterator iterPos = std::find(vertexPatch.begin(), vertexPatch.end(), connectedVertexIndex);
				if (iterPos == vertexPatch.end())
					continue;

				const Vertex& connectedVertex = V.at(connectedVertexIndex);
				if (dir == X_AXIS)
				{
					if (fabs(connectedVertex.x - vertex.x) < 1e-6)
					{
						connectedRegion.push_back(connectedVertexIndex);
					}
				}
				else if (dir == Y_AXIS)
				{
					if (fabs(connectedVertex.y - vertex.y) < 1e-6)
					{
						connectedRegion.push_back(connectedVertexIndex);
					}
				}
				else if (dir == Z_AXIS)
				{
					if (fabs(connectedVertex.z - vertex.z) < 1e-6)
					{
						connectedRegion.push_back(connectedVertexIndex);
					}
				}
			}
		}
		std::sort(connectedRegion.begin(), connectedRegion.end());
		std::vector<unsigned long>::iterator iter = std::unique(connectedRegion.begin(), connectedRegion.end());
		connectedRegion.resize(std::distance(connectedRegion.begin(), iter));
		connectedPatch = connectedRegion;
	}

	if (connectedPatch.size() > 3)
	vertex_patches.push_back(connectedPatch);

	if (connectedPatch.size() == vertexPatch.size())
		return 1;
	else
	{
		std::vector<unsigned long> annother_connectedPatch = vertexPatch;
		for (std::vector<unsigned long>::iterator iter = annother_connectedPatch.begin(); iter != annother_connectedPatch.end(); /*++iter*/)
		{
			bool bFound = false;
			for (int j = 0; j < connectedPatch.size(); j++)
			{
				if (*iter == connectedPatch.at(j))
				{
					iter = annother_connectedPatch.erase(iter);
					bFound = true;
					break;
				}
			}
			if (!bFound)
				++iter;
		}
		return DividePatch(annother_connectedPatch, vertex_patches, dir) + 1;
	}
}
void PrintPatch(const std::vector<unsigned long>& vertexPatch, const std::vector<std::vector<unsigned long> >& vertex_patches, int i)
{
	if (vertex_patches.size() > 1)
	{
		std::cout << "patch " << i << " has " << vertex_patches.size() << " patches, size " << vertexPatch.size() << ", ";// << std::endl;
		for (int j = 0; j < vertex_patches.size(); j++)
		{
			//std::cout << "---subpatch " << j << "'s size: " <<  vertex_patches.at(j).size() << std::endl;
			std::cout  << j << ":" << vertex_patches.at(j).size() << " ";
		}
		std::cout << std::endl;
	}
}
void PolycubeMesh::CheckPatchesConnection()
{
	m_vertex_patches.clear();
	m_vertex_patches_num = 0;
	m_vertex_patches_num_x = 0;
	m_vertex_patches_num_y = 0;
	m_vertex_patches_num_z = 0;

	for (int i = 0; i < m_xLabel; i++)
	{
		const std::vector<unsigned long>& vertexPatch = vertexPatches.at(i);
		std::vector<std::vector<unsigned long> > vertex_patches;
		DividePatch(vertexPatch, vertex_patches, X_AXIS);
		PrintPatch(vertexPatch, vertex_patches, i);

		m_vertex_patches.push_back(vertex_patches);
		m_vertex_patches_num += vertex_patches.size();

		for (int j = 0; j < vertex_patches.size(); j++)
			m_total_vertex_patches.push_back(vertex_patches.at(j));
	}
	m_vertex_patches_num_x = m_vertex_patches_num;

	for (int i = m_xLabel; i < m_yLabel; i++)
	{
		const std::vector<unsigned long>& vertexPatch = vertexPatches.at(i);
		std::vector<std::vector<unsigned long> > vertex_patches;
		DividePatch(vertexPatch, vertex_patches, Y_AXIS);
		PrintPatch(vertexPatch, vertex_patches, i);

		m_vertex_patches.push_back(vertex_patches);
		m_vertex_patches_num += vertex_patches.size();

		for (int j = 0; j < vertex_patches.size(); j++)
			m_total_vertex_patches.push_back(vertex_patches.at(j));
	}
	m_vertex_patches_num_y = m_vertex_patches_num;

	for (int i = m_yLabel; i < m_zLabel; i++)
	{
		const std::vector<unsigned long>& vertexPatch = vertexPatches.at(i);
		std::vector<std::vector<unsigned long> > vertex_patches;
		DividePatch(vertexPatch, vertex_patches, Z_AXIS);
		PrintPatch(vertexPatch, vertex_patches, i);

		m_vertex_patches.push_back(vertex_patches);
		m_vertex_patches_num += vertex_patches.size();

		for (int j = 0; j < vertex_patches.size(); j++)
			m_total_vertex_patches.push_back(vertex_patches.at(j));
	}
	m_vertex_patches_num_z = m_vertex_patches_num;

	std::cout << "##### m_vertex_patches_num is " << m_vertex_patches_num << std::endl;
	std::cout << "##### m_vertex_patches_num_x is " << m_vertex_patches_num_x << std::endl;
	std::cout << "##### m_vertex_patches_num_y is " << m_vertex_patches_num_y << std::endl;
	std::cout << "##### m_vertex_patches_num_z is " << m_vertex_patches_num_z << std::endl;

//	const std::vector<unsigned long>& vertex_patch_ = m_total_vertex_patches.at(0);
//	Patch patch_(*this, vertex_patch_);
//	m_patches.resize(m_total_vertex_patches.size(), patch_);
	for (int i = 0; i < m_total_vertex_patches.size(); i++)
	{
		const std::vector<unsigned long>& vertex_patch = m_total_vertex_patches.at(i);
		Patch patch(*this, vertex_patch);
		if (i < m_vertex_patches_num_x) patch.dir = X_AXIS;
		else if (i >= m_vertex_patches_num_x && i < m_vertex_patches_num_y) patch.dir = Y_AXIS;
		else if (i >= m_vertex_patches_num_y && i < m_vertex_patches_num_z) patch.dir = Z_AXIS;
		patch.GetCorners();
#ifdef _DEBUG
		std::cout << "corners number: " << i << " --- " << patch.m_corners.size() << std::endl;
#endif
		m_patches.push_back(patch);
		//m_patches.at(i) = patch;
	}
}

void PolycubeMesh::SortPatches(const PolycubeMesh& polycubeTetMesh)
{
	// the patch which is the closest to the corners of the patch of polycube Tet mesh is the corresponding patch
	std::vector<Patch> patch;
	for (int i = 0; i < polycubeTetMesh.m_vertex_patches_num_x; i++)
	{
		double min_distance = 100.0;
		unsigned long index = 0;
		for (int j = 0; j < m_vertex_patches_num_x; j++)
		{
			if (m_patches.at(j).m_corners.size() == polycubeTetMesh.m_patches.at(i).m_corners.size())
			{
				double distance = 0;
				for (int k = 0; k < m_patches.at(j).m_corners.size(); k++)
				{
					const Vertex& v = m_patches.at(j).m_corners.at(k) - polycubeTetMesh.m_patches.at(i).m_corners.at(k);
					distance += fabs(v.Length());
				}
				if (distance < min_distance)
				{
					min_distance = distance;
					index = j;
				}
			}
		}
		patch.push_back(m_patches.at(index));
	}

	for (int i = polycubeTetMesh.m_vertex_patches_num_x; i < polycubeTetMesh.m_vertex_patches_num_y; i++)
	{
		double min_distance = 100.0;
		unsigned long index = 0;
		for (int j = m_vertex_patches_num_x; j < m_vertex_patches_num_y; j++)
		{
			if (m_patches.at(j).m_corners.size() == polycubeTetMesh.m_patches.at(i).m_corners.size())
			{
				double distance = 0;
				for (int k = 0; k < m_patches.at(j).m_corners.size(); k++)
				{
					const Vertex& v = m_patches.at(j).m_corners.at(k) - polycubeTetMesh.m_patches.at(i).m_corners.at(k);
					distance += fabs(v.Length());
				}
				if (distance < min_distance)
				{
					min_distance = distance;
					index = j;
				}
			}
		}
		patch.push_back(m_patches.at(index));
	}

	for (int i = polycubeTetMesh.m_vertex_patches_num_y; i < polycubeTetMesh.m_vertex_patches_num_z; i++)
	{
		double min_distance = 100.0;
		unsigned long index = 0;
		for (int j = m_vertex_patches_num_y; j < m_vertex_patches_num_z; j++)
		{
			if (m_patches.at(j).m_corners.size() == polycubeTetMesh.m_patches.at(i).m_corners.size())
			{
				double distance = 0;
				for (int k = 0; k < m_patches.at(j).m_corners.size(); k++)
				{
					const Vertex& v = m_patches.at(j).m_corners.at(k) - polycubeTetMesh.m_patches.at(i).m_corners.at(k);
					distance += fabs(v.Length());
				}
				if (distance < min_distance)
				{
					min_distance = distance;
					index = j;
				}
			}
		}
		patch.push_back(m_patches.at(index));
	}

	patch.resize(patch.size());
	m_patches = patch;
}

void PolycubeMesh::SortPatches()
{
	std::sort(m_patches.begin(), m_patches.end());
}
//////////////////////////////////////
// Get patch labels, including the label whose patches are adjacent to the patch
void PolycubeMesh::GetPatchesLabel()
{
	patchLabel.clear();
	const size_t edgePatches_SIZE = edgePatches.size();
	// for each patch
	for (int i = 0; i < edgePatches_SIZE; i++)
	{
		const std::vector<Edge>& edgePatch = edgePatches.at(i);
		std::vector<int> label;
		for (int j = 0; j < edgePatch.size(); j++)
		{
			const Edge& edge = edgePatch.at(j);
			std::vector<int>& vertexPatch0 = vertexLabel.at(edge.p0);
			std::vector<int>& vertexPatch1 = vertexLabel.at(edge.p1);
			for (int k = 0; k < vertexPatch0.size(); k++)
			{
				label.push_back(vertexPatch0.at(k));
			}
			for (int k = 0; k < vertexPatch1.size(); k++)
			{
				label.push_back(vertexPatch1.at(k));
			}
		}
		std::sort(label.begin(), label.end());
		std::vector<int>::iterator iter = std::unique(label.begin(), label.end());
		size_t size = std::distance(label.begin(), iter);
		label.resize(size);
		patchLabel.push_back(label);
	}
	patchLabel.resize(patchLabel.size());
}
void PolycubeMesh::RelabelIsolatedPatch()
{
	for (int i = 0; i < patchLabel.size(); i++)
	{
		const std::vector<Edge>& edgePatch = edgePatches.at(i);
		std::vector<unsigned long> vec_index;
		for (int j = 0; j < edgePatch.size(); j++)
		{
			const Edge& edge = edgePatch.at(j);
			vec_index.push_back(edge.p0);
			vec_index.push_back(edge.p1);
		}
		std::sort(vec_index.begin(), vec_index.end());
		for (int j = 0; j < vec_index.size() - 2; j++)
		{
			if(vec_index.at(j) == vec_index.at(j + 1) && vec_index.at(j + 1) == vec_index.at(j + 2))
			{
				unsigned long jointIndex = vec_index.at(j);
				std::cout << "IsolatedPatch joint : " << vec_index.at(j) << std::endl;
/*
				// randomly select an edge which is connected to the join. Then walk along the edge path whose the direction is opposite to the join
				// and judge whether the size of vertexlabel is increased. If yes, stop, select another edge which is connected to the join and walk again.
				// If no, then the path walked by is the isolated edge patch.
				// step 1 : find vertices that is connected to the joint. store them in vector.
				const size_t jointLabelSize = vertexLabel.at(jointIndex);
				std::vector<unsigned long> jointConnectedVertices;
				for (int k = 0; k < edgePatch.size(); k++)
				{
					const Edge& edge = edgePatch.at(k);
					if (edge.p0 == jointIndex)
					{
						jointConnectedVertices.push_back(edge.p1);
					}
					else if (edge.p1 == jointIndex)
					{
						jointConnectedVertices.push_back(edge.p0);
					}
				}
				// step 2 : pop a vertex [v] from vector, and find vertices that are connected to [v],
				//          remove those vertices walked by, remaining vertex is the next we walk to.
				std::vector<edge> edgesOfIsolatedPatch;
				for (int l = 0; l < jointConnectedVertices.size(); l++)
				{
					unsigned long walkedByVertex = jointIndex;
					unsigned long currentVertex = jointConnectedVertices.at(l);

					edgesOfIsolatedPatch.clear();
					bool bFoundIsolatedPatch = false;
					for (int k = 0; k < edgePatch.size(); k++)
					{
						// step 3 : judge whether the size of vertexlabel is increased. if yes, break.
						//          otherwise judge whether it is joint, if yes, we get the isolated patch. break;
						if (vertexLabel.at(currentVertex) > jointLabelSize)
						{
							break;
						}

						const Edge& edge = edgePatch.at(k);
						if (edge.p0 == currentVertex || edge.p1 != walkedByVertex)
						{
							Edge isolatedEdge(walkedByVertex, currentVertex);
							edgesOfIsolatedPatch.push_back(isolatedEdge);

							walkedByVertex = currentVertex;
							currentVertex = edge.p1;
						}
						else if (edge.p1 == currentVertex || edge.p0 != walkedByVertex)
						{
							Edge isolatedEdge(walkedByVertex, currentVertex);
							edgesOfIsolatedPatch.push_back(isolatedEdge);

							walkedByVertex = currentVertex;
							currentVertex = edge.p0;
						}

						if (currentVertex == joint)
						{
							Edge isolatedEdge(walkedByVertex, currentVertex);
							edgesOfIsolatedPatch.push_back(isolatedEdge);
							bFoundIsolatedPatch = true;
						}
					}

					if (bFoundIsolatedPatch)
					{
						std::count << "Obtained IsolatedPatch" << std::endl;
						// step 4 : find the edge boundary that is wrapped this isolated patch.
						int correctLabel = 0;
						std::vector<int>& jointLabel = vertexLabel.at(jointIndex);
						const Edge& e = edgesOfIsolatedPatch.at(0);
						for (int n = 0; n < 2; n++)
						{
							bool hole = true;
							bool convex = false;
							const std::vector<unsigned long>& face_patch = facePatches.at(jointLabel.at(n) - 1);
							for (int w = 0; w < face_patch.size(); w++)
							{
								const Face& f = surface.at(face_patch.at(w));
								Edge e0(f[0], f[1]);
								Edge e1(f[1], f[2]);
								Edge e2(f[2], f[0]);
								bool flag = false;
								unsigned long indexOfVertexThatIsInsideTheBoundary = 0;
								// find a point
								if (e0 == e)
								{
									flag = true;
									indexOfVertexThatIsInsideTheBoundary = f[2];
								}
								if (e1 == e)
								{
									flag = true;
									indexOfVertexThatIsInsideTheBoundary = f[0];
								}
								if (e2 == e)
								{
									flag = true;
									indexOfVertexThatIsInsideTheBoundary = f[1];
								}

								if (flag)
								{
									const std::vector<unsigned long>& vertex_patch = vertexPatches.at(jointLabel.at(n) - 1);
									if (std::find(vertex_patch.begin(), vertex_patch.end(), indexOfVertexThatIsInsideTheBoundary) != vertex_patch.end())
									{
										hole = false;
										correctLabel = n == 0? jointLabel.at(1) :jointLabel.at(0);
										std::cout << "correct label is " << correctLabel << std::endl;
										break;
									}
									else
									{
										hole = true;
									}
								}
							}
						}

						// step
						break;
					}
				}
*/
				// already found joint;
				break;
			}
		}
	}
}
void PolycubeMesh::ComputeWedgeTwoEdgePathlengths(const int patchLabelIndex, double* pathLength, int* path_label)
{
	std::vector<unsigned long>& face_patch = facePatches.at(patchLabelIndex);
	int faceIndex = *face_patch.begin();
	int face_label = faceLabel.at(faceIndex);
	int current_path_index = 0;
	for (int j = 0; j < 3; j++)
	{
		// compute edge path lengths
		const int label = patchLabel.at(patchLabelIndex).at(j);
		if (label == face_label)
			continue;
		double path_length = 0.0;
		for (int k = 0; k < edgePatches.at(patchLabelIndex).size(); k++)
		{
			const Edge& edge = edgePatches.at(patchLabelIndex).at(k);
			const std::vector<int>& v0_label = vertexLabel.at(edge.p0);
			const std::vector<int>& v1_label = vertexLabel.at(edge.p1);
			if (std::find(v0_label.begin(), v0_label.end(), label) != v0_label.end()
					&& std::find(v1_label.begin(), v1_label.end(), label) != v1_label.end())
			{
				glm::vec3 v0(V[edge.p0].x, V[edge.p0].y, V[edge.p0].z);
				glm::vec3 v1(V[edge.p1].x, V[edge.p1].y, V[edge.p1].z);
				path_length += glm::length(v0 - v1);
			}
		}
		pathLength[current_path_index] = path_length;
		path_label[current_path_index] = label;
		current_path_index++;
	}
}
void PolycubeMesh::RelabelWedgePatch()
{
	for (int i = 0; i < patchLabel.size(); i++)
	{
		// judge whether it is a wedge
		std::vector<unsigned long>& face_patch = facePatches.at(i);
		if (patchLabel.at(i).size() == 3)
		{
			double pathLength[2] = {0.0, 0.0};
			int  path_label[2] = {0, 0};
			ComputeWedgeTwoEdgePathlengths(i, pathLength, path_label);

			int faceIndex = *face_patch.begin();
			int face_label = faceLabel.at(faceIndex);
			FACE_TYPE face_type = faceType.at(faceIndex);

			int max_path_label = pathLength[0] > pathLength[1] ? path_label[0] : path_label[1];
			int min_path_label = pathLength[0] > pathLength[1] ? path_label[1] : path_label[0];
			std::cout << "max_path_label " << max_path_label << std::endl;
			std::cout << "min_path_label " << min_path_label << std::endl;
			// modify patch face's label
			int max_path_patch_index = max_path_label - 1;
			std::vector<unsigned long>& max_path_face_patch = facePatches.at(max_path_patch_index);
			unsigned long max_path_face_index = *max_path_face_patch.begin();
			FACE_TYPE max_path_face_type = faceType.at(max_path_patch_index);
			const Face& max_path_face = surface.at(max_path_face_index);
			const Vertex& v0 = V.at(max_path_face.at(0));
			for (int j = 0; j < face_patch.size(); j++)
			{
				max_path_face_patch.push_back(face_patch.at(j));
			}
			for (int j = 0; j < faceLabel.size(); j++)
			{
				const Face& f = surface.at(j);
				Vertex& p0 = V.at(f.at(0));
				Vertex& p1 = V.at(f.at(1));
				Vertex& p2 = V.at(f.at(2));
				if (faceLabel.at(j) == face_label)
				{
					///////////////////////////////
					// change vertexInfo clear bBoundaryLine, clear bCorner
					VertexInfo& vi0 = V.at(f.at(0)).vinfo;
					VertexInfo& vi1 = V.at(f.at(1)).vinfo;
					VertexInfo& vi2 = V.at(f.at(2)).vinfo;

//					if (vi0.bCorner) std::cout << "Corner " << f.at(0) << std::endl;
//					if (vi1.bCorner) std::cout << "Corner " << f.at(1) << std::endl;
//					if (vi2.bCorner) std::cout << "Corner " << f.at(2) << std::endl;
//					if (vi0.bBoundaryLine) std::cout << "Boundaryline " << f.at(0) << std::endl;
//					if (vi1.bBoundaryLine) std::cout << "Boundaryline " << f.at(1) << std::endl;
//					if (vi2.bBoundaryLine) std::cout << "Boundaryline " << f.at(2) << std::endl;
//					vi0.bCorner = false;
//					vi1.bCorner = false;
//					vi2.bCorner = false;
					std::vector<int>& vertex_label0 = vertexLabel.at(f.at(0));
					std::vector<int>& vertex_label1 = vertexLabel.at(f.at(1));
					std::vector<int>& vertex_label2 = vertexLabel.at(f.at(2));
					std::vector<int>::iterator iter0 = std::find(vertex_label0.begin(), vertex_label0.end(), max_path_label);
					std::vector<int>::iterator iter1 = std::find(vertex_label1.begin(), vertex_label1.end(), max_path_label);
					std::vector<int>::iterator iter2 = std::find(vertex_label2.begin(), vertex_label2.end(), max_path_label);

					std::vector<int>::iterator iter0m = std::find(vertex_label0.begin(), vertex_label0.end(), min_path_label);
					std::vector<int>::iterator iter1m = std::find(vertex_label1.begin(), vertex_label1.end(), min_path_label);
					std::vector<int>::iterator iter2m = std::find(vertex_label2.begin(), vertex_label2.end(), min_path_label);
					if (iter0 != vertex_label0.end() && iter0m == vertex_label0.end())
					{
						if (vi0.bBoundaryLine)
						std::cout << "Clear Boundaryline " << f.at(0) << std::endl;
						vi0.bBoundaryLine = false;
						vertex_label0.erase(iter0);
					}
					if (iter1 != vertex_label1.end() && iter1m == vertex_label1.end())
					{
						if (vi1.bBoundaryLine)
						std::cout << "Clear Boundaryline " << f.at(1) << std::endl;
						vi1.bBoundaryLine = false;
						vertex_label1.erase(iter1);
					}
					if (iter2 != vertex_label2.end() && iter2m == vertex_label2.end())
					{
						if (vi2.bBoundaryLine)
						std::cout << "Clear Boundaryline " << f.at(2) << std::endl;
						vi2.bBoundaryLine = false;
						vertex_label2.erase(iter2);
					}
					///////////////////////////////
					if (iter0 != vertex_label0.end() && iter0m != vertex_label0.end())
					{
						if (vi0.bCorner){
							std::vector<unsigned long>::iterator iterCorner = std::find(m_vec_xyPlaneCorner.begin(), m_vec_xyPlaneCorner.end(), f.at(0));
							m_vec_xyPlaneCorner.erase(iterCorner);
						std::cout << "Clear Corner " << f.at(0) << std::endl;
						}
						vi0.bCorner = false;
						vertex_label0.erase(iter0);
					}
					if (iter1 != vertex_label1.end() && iter1m != vertex_label1.end())
					{
						if (vi1.bCorner){
							std::vector<unsigned long>::iterator iterCorner = std::find(m_vec_xyPlaneCorner.begin(), m_vec_xyPlaneCorner.end(), f.at(1));
							m_vec_xyPlaneCorner.erase(iterCorner);
						std::cout << "Clear Corner " << f.at(1) << std::endl;
						}
						vi1.bCorner = false;
						vertex_label1.erase(iter1);
					}
					if (iter2 != vertex_label2.end() && iter2m != vertex_label2.end())
					{
						if (vi2.bCorner){
							std::vector<unsigned long>::iterator iterCorner = std::find(m_vec_xyPlaneCorner.begin(), m_vec_xyPlaneCorner.end(), f.at(2));
							m_vec_xyPlaneCorner.erase(iterCorner);
						std::cout << "Clear Corner " << f.at(2) << std::endl;
						}
						vi2.bCorner = false;
						vertex_label2.erase(iter2);
					}
					///////////////////////////////
					faceLabel.at(j) = max_path_label;
					if (max_path_face_type == FACE_X)
					{
						p0.x = v0.x;
						p1.x = v0.x;
						p2.x = v0.x;
					}
					else if (max_path_face_type == FACE_Y)
					{
						p0.y = v0.y;
						p1.y = v0.y;
						p2.y = v0.y;
					}
					else if (max_path_face_type == FACE_Z)
					{
						p0.z = v0.z;
						p1.z = v0.z;
						p2.z = v0.z;
					}
				}
			}
			///////////////////////////////////////////////////////////////
			std::vector<unsigned long> edge_vertex;
			const std::vector<Edge>& edge_path = edgePatches.at(max_path_label - 1);
			for (int n = 0; n < edge_path.size(); n++)
			{
				const Edge& edge = edge_path.at(n);
				if (edge.p0 == 8163 || edge.p1 == 8163)
				{
					std::cout << "8163 is on edgepatch's boundary" << std::endl;
				}
				if (edge.p0 == 7590 || edge.p1 == 7590)
				{
					std::cout << "7590 is on edgepatch's boundary" << std::endl;
				}
				edge_vertex.push_back(edge.p0);
				edge_vertex.push_back(edge.p1);
			}
			std::sort(edge_vertex.begin(), edge_vertex.end());
			std::vector<unsigned long>::iterator iterunique = std::unique(edge_vertex.begin(), edge_vertex.end());
			edge_vertex.resize(std::distance(edge_vertex.begin(), iterunique));
			std::cout << "edgepatch " << max_path_label << " vertex :";
			for (int n = 0; n < edge_vertex.size(); n++)
			{
				std::cout << edge_vertex.at(n) << ", ";
			}
			std::cout << std::endl;
			face_patch.clear();
			///////////////////////////////////////////////////////////////
//			for (int j = 0; j < patchLabel.size(); j++)
//			{
//				std::vector<int>& patch = patchLabel.at(j);
//				const std::vector<int>::iterator iter = std::find(patch.begin(), patch.end(), face_label);
//				if (iter != patch.end())
//				{
//					patch.erase(iter);
//				}
//			}
		}
	}
	GetEdgePatches();
	GetVertexPatches();
	GetPatchesLabel();

	//SmoothBoundaryLine(43, 1);
}
void PolycubeMesh::RelabelSurfaceFacePerFace()
{
	// Get three neighbor triangle's faceType
	for (int i = 0; i < surface.size(); i++)
	{
		const Face& tri = surface.at(i);
		int face_label[3] = {0, 0, 0};
		bool bFoundEdge0 = false;
		bool bFoundEdge1 = false;
		bool bFoundEdge2 = false;
		for (int j = 0; j < surface.size(); j++)
		{
			bool bFoundPoint0 = false;
			bool bFoundPoint1 = false;
			bool bFoundPoint2 = false;
			const Face& face = surface.at(j);
			if (std::find(face.begin(), face.end(), tri.at(0)) != face.end())
			{
				bFoundPoint0 = true;
			}
			if (std::find(face.begin(), face.end(), tri.at(1)) != face.end())
			{
				bFoundPoint1 = true;
			}
			if (std::find(face.begin(), face.end(), tri.at(2)) != face.end())
			{
				bFoundPoint2 = true;
			}
			if (bFoundPoint0 && bFoundPoint1 && bFoundPoint2)
			{
				continue;
			}
			if (bFoundPoint0 && bFoundPoint1)
			{
				face_label[0] = faceLabel[j];
				bFoundEdge0 = true;
			}
			else if (bFoundPoint1 && bFoundPoint2)
			{
				face_label[1] = faceLabel[j];
				bFoundEdge0 = true;
			}
			else if (bFoundPoint0 && bFoundPoint2)
			{
				face_label[2] = faceLabel[j];
				bFoundEdge0 = true;
			}

			if (bFoundEdge0 && bFoundEdge1 && bFoundEdge2)
			{
				break;
			}
		}

		if (face_label[0] == face_label[1])
		{
			if (faceLabel.at(i) != face_label[0])
			{
				std::cout << " origin face label " << faceLabel.at(i) << " change: " << face_label[0] << "," << face_label[1] << "," << face_label[2] << std::endl;
				faceLabel.at(i) = face_label[0];
				std::vector<unsigned long>& face_patch = facePatches.at(face_label[0] - 1);
				face_patch.push_back(i);
			}
		}
		else if (face_label[1] == face_label[2])
		{
			if (faceLabel.at(i) != face_label[1])
			{
				std::cout << " origin face label " << faceLabel.at(i) << " change: " << face_label[0] << "," << face_label[1] << "," << face_label[2] << std::endl;
				faceLabel.at(i) = face_label[1];
				std::vector<unsigned long>& face_patch = facePatches.at(face_label[1] - 1);
				face_patch.push_back(i);
			}
		}
		else if (face_label[2] == face_label[0])
		{
			if (faceLabel.at(i) != face_label[2])
			{
				std::cout << " origin face label " << faceLabel.at(i) << " change: " << face_label[0] << "," << face_label[1] << "," << face_label[2] << std::endl;
				faceLabel.at(i) = face_label[2];
				std::vector<unsigned long>& face_patch = facePatches.at(face_label[2] - 1);
				face_patch.push_back(i);
			}
		}
	}

	GetEdgePatches();
	GetVertexPatches();
	GetPatchesLabel();
}
//void PolycubeMesh::RelabelSurfaceFacePerFace()
//{
//	// Get three neighbor triangle's faceType
//	for (int i = 0; i < surface.size(); i++)
//	{
//		const Cell& tri = surface.at(i);
//		Edge e1(tri.at(0), tri.at(1));
//		Edge e2(tri.at(1), tri.at(2));
//		Edge e3(tri.at(2), tri.at(0));
//		FACE_TYPE face_type[3] =
//		{ FACE_UNKNOWN, FACE_UNKNOWN, FACE_UNKNOWN };
//		std::pair<std::multimap<unsigned long, Face>::iterator, std::multimap<unsigned long, Face>::iterator> ret1 = V_F.equal_range(tri.at(0));
//		for (std::multimap<unsigned long, Face>::iterator iter = ret1.first; iter != ret1.second; ++iter)
//		{
//			Face& f = iter->second;
//			if (IsFaceOnSurface(f))
//			{
//				if (std::find(f.begin(), f.end(), tri.at(0)) != f.end() && std::find(f.begin(), f.end(), tri.at(1)) != f.end()
//						&& std::find(f.begin(), f.end(), tri.at(2)) != f.end())
//					continue;
//				if (std::find(f.begin(), f.end(), e1.p0) != f.end() && std::find(f.begin(), f.end(), e1.p1) != f.end())
//				{
//					const Vertex* pVertex0 = &V[f.at(0)];
//					const Vertex* pVertex1 = &V[f.at(1)];
//					const Vertex* pVertex2 = &V[f.at(2)];
//
//					const glm::vec3 v0(pVertex0->x, pVertex0->y, pVertex0->z);
//					const glm::vec3 v1(pVertex1->x, pVertex1->y, pVertex1->z);
//					const glm::vec3 v2(pVertex2->x, pVertex2->y, pVertex2->z);
//
//					const glm::vec3 v10 = v0 - v1;
//					const glm::vec3 v12 = v2 - v1;
//					const glm::vec3 normal = glm::cross(v12, v10);
//
//					face_type[0] = GetFaceType(normal);
//					break;
//				}
//			}
//		}
//		std::pair<std::multimap<unsigned long, Face>::iterator, std::multimap<unsigned long, Face>::iterator> ret2 = V_F.equal_range(
//				tri.at(1));
//		for (std::multimap<unsigned long, Face>::iterator iter = ret2.first; iter != ret2.second; ++iter)
//		{
//			Face& f = iter->second;
//			if (IsFaceOnSurface(f))
//			{
//				if (std::find(f.begin(), f.end(), tri.at(0)) != f.end() && std::find(f.begin(), f.end(), tri.at(1)) != f.end()
//						&& std::find(f.begin(), f.end(), tri.at(2)) != f.end())
//					continue;
//				if (std::find(f.begin(), f.end(), e2.p0) != f.end() && std::find(f.begin(), f.end(), e2.p1) != f.end())
//				{
//					const Vertex* pVertex0 = &V[f.at(0)];
//					const Vertex* pVertex1 = &V[f.at(1)];
//					const Vertex* pVertex2 = &V[f.at(2)];
//
//					const glm::vec3 v0(pVertex0->x, pVertex0->y, pVertex0->z);
//					const glm::vec3 v1(pVertex1->x, pVertex1->y, pVertex1->z);
//					const glm::vec3 v2(pVertex2->x, pVertex2->y, pVertex2->z);
//
//					const glm::vec3 v10 = v0 - v1;
//					const glm::vec3 v12 = v2 - v1;
//					const glm::vec3 normal = glm::cross(v12, v10);
//
//					face_type[1] = GetFaceType(normal);
//					break;
//				}
//			}
//		}
//
//		std::pair<std::multimap<unsigned long, Face>::iterator, std::multimap<unsigned long, Face>::iterator> ret3 = V_F.equal_range(
//				tri.at(0));
//		for (std::multimap<unsigned long, Face>::iterator iter = ret3.first; iter != ret3.second; ++iter)
//		{
//			Face& f = iter->second;
//			if (IsFaceOnSurface(f))
//			{
//				if (std::find(f.begin(), f.end(), tri.at(0)) != f.end() && std::find(f.begin(), f.end(), tri.at(1)) != f.end()
//						&& std::find(f.begin(), f.end(), tri.at(2)) != f.end())
//					continue;
//				if (std::find(f.begin(), f.end(), e3.p0) != f.end() && std::find(f.begin(), f.end(), e3.p1) != f.end())
//				{
//					const Vertex* pVertex0 = &V[f.at(0)];
//					const Vertex* pVertex1 = &V[f.at(1)];
//					const Vertex* pVertex2 = &V[f.at(2)];
//
//					const glm::vec3 v0(pVertex0->x, pVertex0->y, pVertex0->z);
//					const glm::vec3 v1(pVertex1->x, pVertex1->y, pVertex1->z);
//					const glm::vec3 v2(pVertex2->x, pVertex2->y, pVertex2->z);
//
//					const glm::vec3 v10 = v0 - v1;
//					const glm::vec3 v12 = v2 - v1;
//					const glm::vec3 normal = glm::cross(v12, v10);
//
//					face_type[2] = GetFaceType(normal);
//					break;
//				}
//			}
//		}
//		//			std::cout << " face :" << faceType[0] << "," << faceType[1] << "," << faceType[2] <<std::endl;
//		int count[7] =
//		{ 0, 0, 0, 0, 0, 0, 0 };
//		for (int j = 0; j < 3; j++)
//		{
//			count[face_type[j]]++;
//		}
//		for (int j = 0; j < FACE_OTHER; j++)
//		{
//			if (count[j] >= 2)
//			{
//				if (faceType.at(i) != j)
//				{
//					std::cout << " origin face index " << i << " face type :" << (int)faceType.at(i) << " change: " << face_type[0] << "," << face_type[1]
//							<< "," << face_type[2] << std::endl;
//				}
//				faceType.at(i) = (FACE_TYPE) j;
//				break;
//			}
//		}
//	}
//}

void PolycubeMesh::HexMeshing(const char* outputPolycubeHexFilename, const float delta, const DIRECTION dir)
{
	for (std::vector<Vertex>::iterator iterPoint = V.begin(); iterPoint != V.end(); ++iterPoint)
	{
		GeoUtil::GetMaxMinCoordinateValueOfVertex(*iterPoint, m_maxVertex, m_minVertex);
	}
	Vector axis(m_maxVertex.x - m_minVertex.x, m_maxVertex.y - m_minVertex.y, m_maxVertex.z - m_minVertex.z);
	DivideHexahedronIntoSmallerOnes(m_minVertex, axis, delta, delta, delta, m_vecPolycubeHexVertex, m_vecPolycubeHexCell);
	removeInvalidHexahedronsInBox_new(delta, dir);
	WriteHexahedralmesh(outputPolycubeHexFilename);
}

void PolycubeMesh::GetPatchMinDistance()
{
	std::vector<double> patchPos;
//	for (int i = 0; i < patchLabel.size(); i++)
//	{
//		std::vector<unsigned long>& facePatch = facePatches.at(i);
//		std::vector<unsigned long>& vetexPatch = vertexPatches.at(i);
//		double x = 0.0, y = 0.0, z = 0.0;
//		const size_t vetexPatchSize = vetexPatch.size();
//		for (int j = 0; j < vetexPatchSize; j++)
//		{
//			const Vertex& v = V.at(vetexPatch.at(j));
//			x += v.x;
//			y += v.y;
//			z += v.z;
//		}
//		if (vetexPatchSize != 0)
//		{
//			x /= vetexPatchSize;
//			y /= vetexPatchSize;
//			z /= vetexPatchSize;
//		    //----------------------------------------
//			const FACE_TYPE face_Type = faceType.at(facePatch.at(0));
//
//			if (face_Type == FACE_X)
//			{
//				patchPos.push_back(x);
//			}
//			else if (face_Type == FACE_Y)
//			{
//				patchPos.push_back(y);
//			}
//			else if (face_Type == FACE_Z)
//			{
//				patchPos.push_back(z);
//			}
//		}
//	}

	///////////
	for (std::vector<Patch>::iterator iter = m_patches.begin(); iter != m_patches.end(); ++iter)
	{
		if (iter->dir == X_AXIS)
		{
			patchPos.push_back(iter->center.x);
		}
		else if (iter->dir == Y_AXIS)
		{
			patchPos.push_back(iter->center.y);
		}
		else if (iter->dir == Z_AXIS)
		{
			patchPos.push_back(iter->center.z);
		}
	}
	///////////
	double min_distance = 1;
	double min_distance_x = 1;
	double min_distance_y = 1;
	double min_distance_z = 1;
	double sum_d = 0;
	for (int i = 1; i < m_xLabel; i++)
	{
		double distance = patchPos.at(i) - patchPos.at(i - 1);
//		std::cout << "patch distance = " << patchPos.at(i) << " - " << patchPos.at(i - 1) << " = " << distance << std::endl;
//		m_xpatchInterval.push_back(distance);
		if (distance > 1e-3)
		{
			m_xpatchInterval.push_back(distance);
			sum_d = 0;
		}
		else
			if (!m_xpatchInterval.empty())	m_xpatchInterval.at(m_xpatchInterval.size() - 1) += distance + sum_d;
			else sum_d += distance;
		if (distance < min_distance && distance > 1e-3)
			min_distance = distance;
		if (distance < min_distance_x && distance > 1e-3)
			min_distance_x = distance;
	}
	sum_d = 0;
	for (int i = m_xLabel + 1; i < m_yLabel; i++)
	{
		double distance = patchPos.at(i) - patchPos.at(i - 1);
//		std::cout << "patch distance = " << patchPos.at(i) << " - " << patchPos.at(i - 1) << " = " << distance << std::endl;
//		m_ypatchInterval.push_back(distance);
		if (distance > 1e-3)
		{
			m_ypatchInterval.push_back(distance);
			sum_d = 0;
		}
		else
			if (!m_ypatchInterval.empty())	m_ypatchInterval.at(m_ypatchInterval.size() - 1) += distance + sum_d;
			else sum_d += distance;
		if (distance < min_distance && distance > 1e-3)
			min_distance = distance;
		if (distance < min_distance_y && distance > 1e-3)
			min_distance_y = distance;
	}
	sum_d = 0;
	for (int i = m_yLabel + 1; i < m_zLabel; i++)
	{
		double distance = patchPos.at(i) - patchPos.at(i - 1);
//		std::cout << "patch distance = " << patchPos.at(i) << " - " << patchPos.at(i - 1) << " = " << distance << std::endl;
//		m_zpatchInterval.push_back(distance);
		if (distance > 1e-3)
		{
			m_zpatchInterval.push_back(distance);
			sum_d = 0;
		}
		else if (sum_d > 1e-3)
		{
			m_zpatchInterval.push_back(sum_d);
			sum_d = 0;
		}
		else
		{
			if (!m_zpatchInterval.empty())	m_zpatchInterval.at(m_zpatchInterval.size() - 1) += distance;
			else sum_d += distance;
		}
		if (distance < min_distance && distance > 1e-3)
			min_distance = distance;
		if (distance < min_distance_z && distance > 1e-3)
			min_distance_z = distance;
	}
	cout << "patch min_distance = " << min_distance << std::endl;
	cout << "patch min_distance_x = " << min_distance_x << std::endl;
	cout << "patch min_distance_y = " << min_distance_y << std::endl;
	cout << "patch min_distance_z = " << min_distance_z << std::endl;

	m_patchPos.clear();
	m_patchPos = patchPos;
	m_min_distance = min_distance;
	m_min_distance_x = min_distance_x;
	m_min_distance_y = min_distance_y;
	m_min_distance_z = min_distance_z;
}

void PolycubeMesh::GetPatchMinDistance_round()
{
    std::vector<double> patchPos;
    for (std::vector<Patch>::iterator iter = m_patches.begin(); iter != m_patches.end(); ++iter)
        if (iter->dir == X_AXIS)
            patchPos.push_back(iter->center.x);
        else if (iter->dir == Y_AXIS)
            patchPos.push_back(iter->center.y);
        else if (iter->dir == Z_AXIS)
            patchPos.push_back(iter->center.z);
    ///////////
    double min_distance = 1;
    double min_distance_x = 1;
    double min_distance_y = 1;
    double min_distance_z = 1;
    double sum_d = 0;
    for (int i = 1; i < m_xLabel; i++)
    {
        double distance = patchPos.at(i) - patchPos.at(i - 1);
//      std::cout << "patch distance = " << patchPos.at(i) << " - " << patchPos.at(i - 1) << " = " << distance << std::endl;
//      m_xpatchInterval.push_back(distance);
        if (distance > 1e-3)
        {
            m_xpatchInterval.push_back(distance);
            sum_d = 0;
        }
        else
            if (!m_xpatchInterval.empty())  m_xpatchInterval.at(m_xpatchInterval.size() - 1) += distance + sum_d;
            else sum_d += distance;
        if (distance < min_distance && distance > 1e-3)
            min_distance = distance;
        if (distance < min_distance_x && distance > 1e-3)
            min_distance_x = distance;
    }
    sum_d = 0;
    for (int i = m_xLabel + 1; i < m_yLabel; i++)
    {
        double distance = patchPos.at(i) - patchPos.at(i - 1);
//      std::cout << "patch distance = " << patchPos.at(i) << " - " << patchPos.at(i - 1) << " = " << distance << std::endl;
//      m_ypatchInterval.push_back(distance);
        if (distance > 1e-3)
        {
            m_ypatchInterval.push_back(distance);
            sum_d = 0;
        }
        else
            if (!m_ypatchInterval.empty())  m_ypatchInterval.at(m_ypatchInterval.size() - 1) += distance + sum_d;
            else sum_d += distance;
        if (distance < min_distance && distance > 1e-3)
            min_distance = distance;
        if (distance < min_distance_y && distance > 1e-3)
            min_distance_y = distance;
    }
    sum_d = 0;
    for (int i = m_yLabel + 1; i < m_zLabel; i++)
    {
        double distance = patchPos.at(i) - patchPos.at(i - 1);
//      std::cout << "patch distance = " << patchPos.at(i) << " - " << patchPos.at(i - 1) << " = " << distance << std::endl;
//      m_zpatchInterval.push_back(distance);
        if (distance > 1e-3)
        {
            m_zpatchInterval.push_back(distance);
            sum_d = 0;
        }
        else if (sum_d > 1e-3)
        {
            m_zpatchInterval.push_back(sum_d);
            sum_d = 0;
        }
        else
        {
            if (!m_zpatchInterval.empty())  m_zpatchInterval.at(m_zpatchInterval.size() - 1) += distance;
            else sum_d += distance;
        }
        if (distance < min_distance && distance > 1e-3)
            min_distance = distance;
        if (distance < min_distance_z && distance > 1e-3)
            min_distance_z = distance;
    }
    cout << "patch min_distance = " << min_distance << std::endl;
    cout << "patch min_distance_x = " << min_distance_x << std::endl;
    cout << "patch min_distance_y = " << min_distance_y << std::endl;
    cout << "patch min_distance_z = " << min_distance_z << std::endl;

    m_patchPos.clear();
    m_patchPos = patchPos;
    m_min_distance = min_distance;
    m_min_distance_x = min_distance_x;
    m_min_distance_y = min_distance_y;
    m_min_distance_z = min_distance_z;
}

double g_limit = 0.008;
double g_upper_size = 0.008;
double g_lower_size = 0.004;

static void GetInterval(const double m_min_distance, const std::vector<double>& patchInterval, std::vector<double>& hexGridInterval)
{
	for (int i = 0; i < patchInterval.size(); i++)
	{
		double choise = 0.05;
		const double d = patchInterval.at(i);
		if (g_lower_size == g_upper_size)
		{
			double max = g_limit > m_min_distance ? g_limit : m_min_distance;
//			choise = patchInterval.at(i) > max ? max : m_min_distance;
			choise = patchInterval.at(i) > max ? max : patchInterval.at(i);
		}
		else
		{
			if (d < m_min_distance) choise = m_min_distance;
			else if (d > m_min_distance && d < g_lower_size) choise = d;
			else if (d > g_lower_size && d < g_upper_size) choise = g_lower_size;
			else if (d > g_upper_size) choise = g_upper_size;
		}

		int num = round(patchInterval.at(i)/choise);
		double interval = patchInterval.at(i)/num;
		for (int j = 0; j < num; j++)
		{
			if (j == num - 1) interval = patchInterval.at(i) - (num - 1)*interval;
			hexGridInterval.push_back(interval);
		}
	}
}

static void GetInterval_aniso(const double m_min_distance, const std::vector<double>& patchInterval, std::vector<double>& hexGridInterval,
        const double scale = 1.0, const float delta = 0.002f)
{
    for (int i = 0; i < patchInterval.size(); i++)
    {
        int num = ceil(patchInterval.at(i)/(delta * scale));
        double interval = patchInterval.at(i)/num;
        for (int j = 0; j < num; j++)
        {
            if (j == num - 1) interval = patchInterval.at(i) - (num - 1)*interval;
            hexGridInterval.push_back(interval);
        }
    }
}
unsigned int PolycubeMesh::DividePolycube(const Vertex& originPoint, std::vector<Vertex>& vecVertex, std::vector<Cell>& vecHexahedron)
{
	GetPatchMinDistance();

	g_limit = g_upper_size;
	unsigned int prevPointNumber = vecVertex.size();
	unsigned int PointNumber = prevPointNumber;
	std::vector<double> xInterval;
	std::vector<double> yInterval;
	std::vector<double> zInterval;

//	if (m_min_distance < 0.002)
//	{
//		m_min_distance = 0.002;
//	}
	if (m_min_distance > g_lower_size)
	{
		swap(m_min_distance, g_lower_size);
	}
	if (m_min_distance_x > g_lower_size)
	{
		m_min_distance_x = g_lower_size;
	}
	if (m_min_distance_y > g_lower_size)
	{
		m_min_distance_y = g_lower_size;
	}
	if (m_min_distance_z > g_lower_size)
	{
		m_min_distance_z = g_lower_size;
	}
	GetInterval(m_min_distance_x, m_xpatchInterval, xInterval);
	GetInterval(m_min_distance_y, m_ypatchInterval, yInterval);
	GetInterval(m_min_distance_z, m_zpatchInterval, zInterval);

	unsigned int xTotalStep = xInterval.size() + 1;
	unsigned int yTotalStep = yInterval.size() + 1;
	unsigned int zTotalStep = zInterval.size() + 1;
	unsigned int yzPlaneNumber = yTotalStep * zTotalStep;

	double xPos = m_patchPos.at(0);
	double yPos = m_patchPos.at(m_xLabel);
	double zPos = m_patchPos.at(m_yLabel);

//	std::cout << " x patches: " << m_xLabel - 1 << std::endl;
//	std::cout << " y patches: " << m_yLabel - 1 << std::endl;
//	std::cout << " z patches: " << m_zLabel - 1 << std::endl;
//	double xPos = originPoint.x;
//	double yPos = originPoint.y;
//	double zPos = originPoint.z;
//	std::cout << "*************************************************" << std::endl;
//	std::cout << " xPos: " << xPos << std::endl;
//	std::cout << " yPos: " << yPos << std::endl;
//	std::cout << " zPos: " << zPos << std::endl;
//	std::cout << "*************************************************" << std::endl;

	for (unsigned int xStep = 0; xStep < xTotalStep; xStep++)
	{
		if (xStep > 0)
		{	xPos += xInterval.at(xStep - 1);
//			if (xInterval.at(xStep - 1) < 0.01)
//			{
//				continue;
//			}
		}
		for (unsigned int yStep = 0; yStep < yTotalStep; yStep++)
		{
			if (yStep == 0) yPos = originPoint.y;
			if (yStep > 0)
			{
				yPos += yInterval.at(yStep - 1);
//				if (yInterval.at(yStep - 1) < 0.01)
//				{
//					continue;
//				}
			}
			for (unsigned int zStep = 0; zStep < zTotalStep; zStep++)
			{
				if (zStep == 0) zPos = originPoint.z;
				if (zStep > 0)
				{
					zPos += zInterval.at(zStep - 1);
//					if (zInterval.at(zStep - 1) < 0.01) continue;
				}
				Vertex v(xPos, yPos, zPos);
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

unsigned int PolycubeMesh::DividePolycube_aniso(const Vertex& originPoint, std::vector<Vertex>& vecVertex, std::vector<Cell>& vecHexahedron,
        const glm::vec3& cellScale, const float delta)
{
    GetPatchMinDistance();

    g_limit = g_upper_size;
    unsigned int prevPointNumber = vecVertex.size();
    unsigned int PointNumber = prevPointNumber;
    std::vector<double> xInterval;
    std::vector<double> yInterval;
    std::vector<double> zInterval;

    if (m_min_distance > g_lower_size)
    {
        swap(m_min_distance, g_lower_size);
    }
    if (m_min_distance_x > g_lower_size)
        m_min_distance_x = g_lower_size;
    if (m_min_distance_y > g_lower_size)
        m_min_distance_y = g_lower_size;
    if (m_min_distance_z > g_lower_size)
        m_min_distance_z = g_lower_size;
    GetInterval_aniso(m_min_distance_x, m_xpatchInterval, xInterval, cellScale.x, delta);
    GetInterval_aniso(m_min_distance_y, m_ypatchInterval, yInterval, cellScale.y, delta);
    GetInterval_aniso(m_min_distance_z, m_zpatchInterval, zInterval, cellScale.z, delta);

    unsigned int xTotalStep = xInterval.size() + 1;
    unsigned int yTotalStep = yInterval.size() + 1;
    unsigned int zTotalStep = zInterval.size() + 1;
    unsigned int yzPlaneNumber = yTotalStep * zTotalStep;

    double xPos = m_patchPos.at(0);
    double yPos = m_patchPos.at(m_xLabel);
    double zPos = m_patchPos.at(m_yLabel);

    for (unsigned int xStep = 0; xStep < xTotalStep; xStep++)
    {
        if (xStep > 0)
        {
            xPos += xInterval.at(xStep - 1);
        }
        for (unsigned int yStep = 0; yStep < yTotalStep; yStep++)
        {
            if (yStep == 0) yPos = originPoint.y;
            if (yStep > 0)
            {
                yPos += yInterval.at(yStep - 1);
            }
            for (unsigned int zStep = 0; zStep < zTotalStep; zStep++)
            {
                if (zStep == 0) zPos = originPoint.z;
                if (zStep > 0)
                {
                    zPos += zInterval.at(zStep - 1);
                }
                Vertex v(xPos, yPos, zPos);
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

unsigned int PolycubeMesh::DividePolycube_round(const Vertex& originPoint, std::vector<Vertex>& vecVertex, std::vector<Cell>& vecHexahedron)
{
    GetPatchMinDistance();

    g_limit = g_upper_size;
    unsigned int prevPointNumber = vecVertex.size();
    unsigned int PointNumber = prevPointNumber;
    std::vector<double> xInterval;
    std::vector<double> yInterval;
    std::vector<double> zInterval;
    if (m_min_distance > g_lower_size)
    {
        swap(m_min_distance, g_lower_size);
    }
    if (m_min_distance_x > g_lower_size)
    {
        m_min_distance_x = g_lower_size;
    }
    if (m_min_distance_y > g_lower_size)
    {
        m_min_distance_y = g_lower_size;
    }
    if (m_min_distance_z > g_lower_size)
    {
        m_min_distance_z = g_lower_size;
    }
    GetInterval(m_min_distance_x, m_xpatchInterval, xInterval);
    GetInterval(m_min_distance_y, m_ypatchInterval, yInterval);
    GetInterval(m_min_distance_z, m_zpatchInterval, zInterval);

    unsigned int xTotalStep = xInterval.size() + 1;
    unsigned int yTotalStep = yInterval.size() + 1;
    unsigned int zTotalStep = zInterval.size() + 1;
    unsigned int yzPlaneNumber = yTotalStep * zTotalStep;

    double xPos = m_patchPos.at(0);
    double yPos = m_patchPos.at(m_xLabel);
    double zPos = m_patchPos.at(m_yLabel);

//  std::cout << " x patches: " << m_xLabel - 1 << std::endl;
//  std::cout << " y patches: " << m_yLabel - 1 << std::endl;
//  std::cout << " z patches: " << m_zLabel - 1 << std::endl;
//  double xPos = originPoint.x;
//  double yPos = originPoint.y;
//  double zPos = originPoint.z;
//  std::cout << "*************************************************" << std::endl;
//  std::cout << " xPos: " << xPos << std::endl;
//  std::cout << " yPos: " << yPos << std::endl;
//  std::cout << " zPos: " << zPos << std::endl;
//  std::cout << "*************************************************" << std::endl;

    for (unsigned int xStep = 0; xStep < xTotalStep; xStep++)
    {
        if (xStep > 0)
        {   xPos += xInterval.at(xStep - 1);
//          if (xInterval.at(xStep - 1) < 0.01)
//          {
//              continue;
//          }
        }
        for (unsigned int yStep = 0; yStep < yTotalStep; yStep++)
        {
            if (yStep == 0) yPos = originPoint.y;
            if (yStep > 0)
            {
                yPos += yInterval.at(yStep - 1);
//              if (yInterval.at(yStep - 1) < 0.01)
//              {
//                  continue;
//              }
            }
            for (unsigned int zStep = 0; zStep < zTotalStep; zStep++)
            {
                if (zStep == 0) zPos = originPoint.z;
                if (zStep > 0)
                {
                    zPos += zInterval.at(zStep - 1);
//                  if (zInterval.at(zStep - 1) < 0.01) continue;
                }
                Vertex v(xPos, yPos, zPos);
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

void PolycubeMesh::GetFacePatches_N()
{
	// std::vector<std::vector<unsigned long> > facePatches;
	// the begin patch's lable is 1, next is 2, and so on. unsigned long is the face index in surface
	facePatches.clear();
	std::vector<unsigned long> a;
	facePatches.resize(500, a);
	int patchIndex = 0;
	m_xLabel = 0;
	m_yLabel = 0;
	m_zLabel = 0;

	//std::vector<unsigned long> faceLabel(surface.size(), 0);
	faceLabel.clear();
	faceLabel.resize(surface.size(), 0);
	const std::multimap<unsigned long, unsigned long>::iterator iterBegin = F_F.begin();
	const std::multimap<unsigned long, unsigned long>::iterator iterEnd = F_F.end();
	std::multimap<unsigned long, unsigned long>::iterator iter = iterBegin;
	for (; iter != iterEnd; ++iter)
	{
		const unsigned long faceId = iter->first;
		if (faceLabel.at(faceId) != 0)
		{
			continue;
		}
		GetFacePatches_N(faceId, ++patchIndex);
		if (faceType.at(faceId) == FACE_X) m_xLabel++;
		else if (faceType.at(faceId) == FACE_Y) m_yLabel++;
		else if (faceType.at(faceId) == FACE_Z) m_zLabel++;
	}
	m_yLabel += m_xLabel;
	m_zLabel += m_yLabel;
	m_patchNumber = patchIndex;
	facePatches.resize(patchIndex);
}

void PolycubeMesh::GetFacePatches_N(const unsigned long faceId, int& patchIndex)
{
	if (faceLabel.at(faceId) != 0)
	{
		return;
	}
	const std::multimap<unsigned long, unsigned long>::iterator iterBegin = F_F.begin();
	const std::multimap<unsigned long, unsigned long>::iterator iterEnd = F_F.end();
	std::multimap<unsigned long, unsigned long>::iterator iter = iterBegin;
//	for (; iter != iterEnd; /*++iter*/)
	{
		std::pair<std::multimap<unsigned long, unsigned long>::iterator, std::multimap<unsigned long, unsigned long>::iterator>
			ret = F_F.equal_range(faceId);
		const std::multimap<unsigned long, unsigned long>::iterator iterLowerBound = ret.first;
		const std::multimap<unsigned long, unsigned long>::iterator iterUpperBound = ret.second;

		faceLabel.at(faceId) = patchIndex;
		std::vector<unsigned long> facePatch;
		//facePatch.push_back(cellId);

		const Face& face = surface.at(faceId);
		Plane facePlane(V.at(face.at(0)), V.at(face.at(1)), V.at(face.at(2)));
		for (iter = iterLowerBound; iter != iterUpperBound; ++iter)
		{
			const unsigned long neighborFaceId = iter->second;
			if (faceLabel.at(neighborFaceId) != 0) continue;
			const Cell& neighborFace = surface.at(neighborFaceId);
			if (neighborFace.size() == 3 || neighborFace.size() == 4)
			{
				const Plane neighborFacePlane(V.at(neighborFace.at(0)), V.at(neighborFace.at(1)), V.at(neighborFace.at(2)));
				const double cosangle = facePlane.IntersectionAngle(neighborFacePlane);
				if (cosangle > 0.866 ) // cos(15) = 0.9659 cos(30) = 0.866
				{
					facePatch.push_back(neighborFaceId);
				}
			}
		}
		for (size_t i = 0; i < facePatch.size(); i++)
		{
			GetFacePatches_N(facePatch.at(i), patchIndex);
		}
	}
}

void PolycubeMesh::GetCellPatches_N()
{
	// std::vector<std::vector<unsigned long> > facePatches;
	// the begin patch's lable is 1, next is 2, and so on. unsigned long is the face index in surface
	facePatches.clear();
	std::vector<unsigned long> a;
	facePatches.resize(500, a);
	int patchIndex = 0;

	//std::vector<unsigned long> faceLabel(surface.size(), 0);
	faceLabel.clear();
	faceLabel.resize(surface.size(), 0);
	const std::multimap<unsigned long, unsigned long>::iterator iterBegin = C_C.begin();
	const std::multimap<unsigned long, unsigned long>::iterator iterEnd = C_C.end();
	std::multimap<unsigned long, unsigned long>::iterator iter = iterBegin;
	for (; iter != iterEnd; ++iter)
	{
		const unsigned long cellId = iter->first;
		if (faceLabel.at(cellId) != 0)
		{
			continue;
		}
		GetFacePatches_N(cellId, ++patchIndex);
	}
	m_patchNumber = patchIndex;
	facePatches.resize(patchIndex);
}

void PolycubeMesh::GetCellPatches_N(const unsigned long cellId, int& patchIndex)
{
	if (faceLabel.at(cellId) != 0)
	{
		return;
	}
	const std::multimap<unsigned long, unsigned long>::iterator iterBegin = C_C.begin();
	const std::multimap<unsigned long, unsigned long>::iterator iterEnd = C_C.end();
	std::multimap<unsigned long, unsigned long>::iterator iter = iterBegin;
//	for (; iter != iterEnd; /*++iter*/)
	{
		std::pair<std::multimap<unsigned long, unsigned long>::iterator, std::multimap<unsigned long, unsigned long>::iterator>
			ret = C_C.equal_range(cellId);
		const std::multimap<unsigned long, unsigned long>::iterator iterLowerBound = ret.first;
		const std::multimap<unsigned long, unsigned long>::iterator iterUpperBound = ret.second;

		faceLabel.at(cellId) = patchIndex;
		std::vector<unsigned long> facePatch;
		//facePatch.push_back(cellId);

		const Cell& cell = C.at(cellId);
		Plane cellPlane(V.at(cell.at(0)), V.at(cell.at(1)), V.at(cell.at(2)));
		for (iter = iterLowerBound; iter != iterUpperBound; ++iter)
		{
			const unsigned long neighborCellId = iter->second;
			if (faceLabel.at(neighborCellId) != 0) continue;
			const Cell& neighborCell = C.at(neighborCellId);
			if (neighborCell.size() == 3)
			{
				const Plane neighborCellPlane(V.at(neighborCell.at(0)), V.at(neighborCell.at(1)), V.at(neighborCell.at(2)));
				const double cosangle = cellPlane.IntersectionAngle(neighborCellPlane);
				if (cosangle > 0.866 ) // cos(15) = 0.9659 cos(30) = 0.866
				{
					facePatch.push_back(neighborCellId);
				}
			}
		}
		for (size_t i = 0; i < facePatch.size(); i++)
		{
			GetFacePatches_N(facePatch.at(i), patchIndex);
		}
	}
}

void PolycubeMesh::GetPatches()
{
	m_patches.clear();
	for (size_t i = 0; i < facePatches.size(); i++)
	{
		Patch patch(*this, facePatches.at(i), edgePatches.at(i), vertexPatches.at(i));
		m_patches.push_back(patch);
	}
}

void PolycubeMesh::ModifyPatchesPosition()
{
	for (size_t i = 0; i < m_patches.size(); i++)
	{
		Patch& patch = m_patches.at(i);
		patch.ModifyPatchesPosition();
	}
}

void PolycubeMesh::ModifyPatchesPosition_ROUND()
{
	for (size_t i = 0; i < m_patches.size(); i++)
	{
		Patch& patch = m_patches.at(i);
		patch.ModifyPatchesPosition_ROUND();
	}
}

void PolycubeMesh::ModifyPatchesPosition_ROUND(float delta)
{
	for (size_t i = 0; i < m_patches.size(); i++)
	{
		Patch& patch = m_patches.at(i);
		patch.ModifyPatchesPosition_ROUND(delta);
	}
}

void PolycubeMesh::ModifyPatchesPosition_ROUND2(float delta)
{
    for (size_t i = 0; i < m_patches.size(); i++)
    {
        Patch& patch = m_patches.at(i);
        patch.ModifyPatchesPosition_ROUND(delta);
    }
}

void PolycubeMesh::AssignPatches()
{
    for (int i = 0; i < m_patches.size(); i++)
    {
        facePatches.at(i) = m_patches.at(i).m_f;
        for (int j = 0; j < facePatches.at(i).size(); j++)
        {
            faceLabel.at(facePatches.at(i).at(j)) = i + 1;
        }
        edgePatches.at(i) = m_patches.at(i).m_e;
        vertexPatches.at(i) = m_patches.at(i).m_v;
    }
}
