/*
 * Parametrizer.h
 *
 *  Created on: Mar 21, 2016
 *      Author: cotrik
 */

#ifndef __PARAMETRIZER_H__
#define __PARAMETRIZER_H__

#include "include.h"

class Parametrizer
{
public:
	Parametrizer();
	virtual ~Parametrizer();
    // Tao Ju. <<Mean Value Coordinates for Closed Triangular Meshes>>
	void Parametrize3DInnerPoint(const Vertex& innerPoint,
			const std::vector<Vertex>& boundaryVertices,
			const std::vector<Cell>& boundaryTriangles,
			std::vector<glm::vec3>& w,
			const double eps = 1e-8);
};

#endif /* __PARAMETRIZER_H__ */
