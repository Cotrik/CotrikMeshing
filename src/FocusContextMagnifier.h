/*
 * FocusContextMagnifier.h
 *
 *  Created on: Feb 10, 2015
 *      Author: cotrik
 */

#ifndef __FOCUSCONTEXTMAGNIFIER_H__
#define __FOCUSCONTEXTMAGNIFIER_H__

#include "Mesh.h"
struct IndexScale
{
	IndexScale()
	{}

	IndexScale(const unsigned long index, const float scale)
	: index(index), scale(scale)
	{}

	unsigned long index;
	float scale;
};
extern const unsigned int enlargeCubeIndices[];
extern const IndexScale indexScale[];
extern std::string strScaleParas;
class FocusContextMagnifier
{
public:
	FocusContextMagnifier(const Mesh& mesh);
	virtual ~FocusContextMagnifier();

public:
	const size_t GetHexVerticesNum() const;
	const size_t GetHexCellsNum() const;
	const size_t GetHexEdgesNum();// const;

	void UpdateVerticesPositions();
	void ConstructLinearEquation();
	void SolveLinearEquation();
	void SolveLinearEquation_TestLLT_PR_LU();
private:
	void Construct_Q();
	void Construct_S();

	void Construct_Vx();
	void Construct_Vy();
	void Construct_Vz();

	void Construct_Rx();
	void Construct_Ry();
	void Construct_Rz();

	void Construct_Hx();
	void Construct_Hy();
	void Construct_Hz();

	void Construct_A();
	void Construct_x();
	void Construct_b();

private:
	Mesh& m_mesh;
	std::vector<Edge> m_vecHexEdge;
};

#endif /* __FOCUSCONTEXTMAGNIFIER_H__ */
