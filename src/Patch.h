/*
 * Patch.h
 *
 *  Created on: Jun 19, 2015
 *      Author: cotrik
 */

#ifndef __PATCH_H__
#define __PATCH_H__

//#include "include.h"
#include "GeometricStruct.h"
//#include "PolycubeMesh.h"
class PolycubeMesh;
class Patch
{
public:
	Patch();

public:
	//const bool operator < (const Patch& right) const;
	void GetCorners();

public:
	Patch(const Patch& patch);
	Patch(const PolycubeMesh& mesh);
	Patch(const PolycubeMesh& mesh, const std::vector<unsigned long>& patchVertex);
	Patch(const PolycubeMesh& mesh, const std::vector<unsigned long>& f, std::vector<Edge>& e, std::vector<unsigned long>& v);
	virtual ~Patch();

public:
	const bool operator < (const Patch& patch) const;
	const bool operator == (const Patch& patch) const;
	void GetPatchPosition();
	void ModifyPatchesPosition();
	void ModifyPatchesPosition_ROUND();
	void ModifyPatchesPosition_ROUND(float delta);
	void ModifyPatchesPosition_ROUND2(float delta);
	void GetParametrization(const Patch& triPatch, TriQuadParas& tqp);
	void MapbackToOrigTri(const std::vector<Vertex>& orgV, const Patch& triPatch, const TriQuadParas& p, std::vector<Vertex>& quadV);
	void Smooth();
	void GetInnerVertices(const std::vector<unsigned long>& boundaryVertices, std::vector<unsigned long>& innerVertices);
	void GetBoundaryVertices(std::vector<unsigned long>& boundaryVertices);
public:
	std::vector<unsigned long> m_v;
	std::vector<Edge> m_e;
	std::vector<unsigned long> m_f;
	std::vector<Vertex> m_corners;

	PolycubeMesh* m_pMesh;

	DIRECTION dir;
	double pos;
	glm::vec3 center;
};

#endif /* __PATCH_H__ */
